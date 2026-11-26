# Peptide consistency weights (vignette Peptide_consistency_weights): weights of the protein x sample abundances from
# how well a protein's peptides agree, the residuals of the median polish that aggregated them. Recomputes the median
# polish of each protein on the prepared peptides (checked against the prepared proteins), derives the weight
# variants, fits the protein-level prolfqua facades lm, limma, deqms and rlm with each, and scores the fits like
# score.R: pAUC_10 and AUC of the signed p-value, all contrasts pooled, and the false discovery proportion and number
# of true discoveries at a nominal FDR of 5%. Writes <out_dir>/scores.tsv (one row per facade and variant) and
# <out_dir>/priors.tsv (the empirical Bayes priors of the dataset), the contrasts under <out_dir>/contrasts.
#   Rscript workflow/scripts/peptide_consistency_weights.R <config> <dataset> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "prepared_dir", "out_dir"))
methods <- c("prolfqua_lm", "prolfqua_limma", "prolfqua_deqms", "prolfqua_rlm")
variants <- c(
  "default",
  "count",
  "count_conf",
  "protein_ivw",
  "sample_ivw",
  "sample_ivw_conf",
  "sample_rel",
  "sample_rel_conf"
)

# Median polish of each protein on its observed peptides, as prolfqua's medpolish aggregation; one row per
# observation with its residual and the protein estimate of its sample
polish <- function(peptides) {
  response <- peptides$response()
  sname <- peptides$sample_name()
  data <- as.data.frame(peptides$data_long())[, c("protein_Id", "peptide_Id", sname, response)]
  data <- data[!is.na(data[[response]]), ]
  fits <- lapply(split(data, data$protein_Id), function(d) {
    features <- unique(d$peptide_Id)
    samples <- unique(d[[sname]])
    i <- match(d$peptide_Id, features)
    j <- match(d[[sname]], samples)
    mat <- matrix(NA_real_, length(features), length(samples))
    mat[cbind(i, j)] <- d[[response]]
    X <- suppressWarnings(stats::medpolish(mat, na.rm = TRUE, trace.iter = FALSE, maxiter = 10, eps = 0.01))
    data.frame(
      protein_Id = d$protein_Id[1],
      sample = d[[sname]],
      residual = X$residuals[cbind(i, j)],
      estimate = (X$col + X$overall)[j],
      n_peptides = length(features),
      n_samples = length(samples)
    )
  })
  dplyr::bind_rows(fits)
}

# Per protein P and sample s: n (observed peptides), RSS (sum of squared residuals) and the residual degrees of
# freedom; the protein variance sigma2_P squeezed over proteins, and the ratio of the sample variance s2_Ps to it
# squeezed over all protein x sample cells (limma::squeezeVar, Smyth 2004). Cells with one peptide get the prior.
consistency <- function(peptides) {
  cells <- dplyr::summarise(
    polish(peptides),
    n = dplyr::n(),
    rss = sum(.data$residual^2),
    estimate = .data$estimate[1],
    .by = c("protein_Id", "sample", "n_peptides", "n_samples")
  )
  protein <- dplyr::summarise(
    cells,
    N = sum(.data$n),
    rss = sum(.data$rss),
    .by = c("protein_Id", "n_peptides", "n_samples")
  )
  # residual df of the additive fit overall + peptide + sample
  protein$d <- pmax(protein$N - protein$n_peptides - protein$n_samples + 1, 0)
  squeezed_protein <- limma::squeezeVar(ifelse(protein$d > 0, protein$rss / protein$d, 0), protein$d)
  protein$sigma2 <- squeezed_protein$var.post
  cells <- dplyr::left_join(cells, protein[, c("protein_Id", "N", "d", "sigma2")], by = "protein_Id")
  # the protein's df split over its samples in proportion to n - 1
  cells$d_cell <- ifelse(cells$N > cells$n_samples, (cells$n - 1) * cells$d / (cells$N - cells$n_samples), 0)
  ratio <- ifelse(cells$d_cell > 0, cells$rss / cells$d_cell / cells$sigma2, 0)
  squeezed_ratio <- limma::squeezeVar(ratio, cells$d_cell)
  cells$ratio <- squeezed_ratio$var.post
  cells$factor <- 1 / cells$ratio
  cells <- dplyr::mutate(cells, factor_rel = .data$factor / mean(.data$factor), .by = "protein_Id")
  list(
    cells = cells,
    protein = protein,
    priors = data.frame(
      proteins = nrow(protein),
      protein_d0 = squeezed_protein$df.prior,
      protein_s02 = squeezed_protein$var.prior,
      cells = nrow(cells),
      cells_with_df = sum(cells$d_cell > 0),
      ratio_d0 = squeezed_ratio$df.prior,
      ratio_s02 = squeezed_ratio$var.prior,
      factor_sd_log2 = stats::median(tapply(log2(cells$factor), cells$protein_Id, stats::sd), na.rm = TRUE)
    )
  )
}

# The prepared proteins with the weight of `variant` in column consistency_weight, named as nr_children so that
# limma takes it as a feature x sample value column; nr_children_feature (the DEqMS count) stays the peptide count.
# The configuration is rebuilt from its values, so prepared files serialized by an older prolfqua take today's fields.
with_weight <- function(proteins, weights, variant) {
  config <- prolfqua::list_to_analysis_configuration(prolfqua::r6_extract_values(proteins$get_config()))
  key <- c("protein_Id", config$sample_name)
  cells <- dplyr::rename(
    weights$cells[, c("protein_Id", "sample", "factor", "factor_rel")],
    !!config$sample_name := "sample"
  )
  data <- proteins$data_long() |>
    dplyr::left_join(cells, by = key) |>
    dplyr::left_join(weights$protein[, c("protein_Id", "sigma2")], by = "protein_Id")
  n <- data[[config$nr_children]]
  confidence <- data[[config$identification_confidence]]
  # rows without observed peptides (no response) keep factor 1
  factor <- ifelse(is.na(data$factor), 1, data$factor)
  factor_rel <- ifelse(is.na(data$factor_rel), 1, data$factor_rel)
  data$consistency_weight <- switch(
    variant,
    count = n,
    count_conf = n * confidence,
    protein_ivw = n / data$sigma2,
    sample_ivw = n * factor / data$sigma2,
    sample_ivw_conf = n * factor / data$sigma2 * confidence,
    sample_rel = n * factor_rel,
    sample_rel_conf = n * factor_rel * confidence
  )
  data <- dplyr::select(data, -dplyr::all_of(c("factor", "factor_rel", "sigma2")))
  config$nr_children <- "consistency_weight"
  config$use_identification_confidence <- FALSE
  prolfqua::LFQData$new(data, config)
}

# FDP and number of true discoveries at the last point with nominal FDR <= 0.05, all contrasts pooled
at_fdr_5 <- function(bench) {
  fdp <- dplyr::filter(bench$get_confusion_FDRvsFDP(), .data$contrast == "all", .data$scorecol <= 0.05)
  if (nrow(fdp) == 0) {
    return(c(FDP_5 = 0, TP_5 = 0))
  }
  last <- fdp[which.max(fdp$scorecol), ]
  c(FDP_5 = last$FDP[1], TP_5 = last$TP_hits[1])
}

peptides <- readRDS(file.path(args$prepared_dir, "peptides.rds"))
proteins <- readRDS(file.path(args$prepared_dir, "proteins.rds"))
weights <- consistency(peptides)

# the recomputed median polish reproduces the prepared proteins, and n their peptide count
sname <- proteins$sample_name()
check <- dplyr::inner_join(
  proteins$data_long(),
  dplyr::rename(weights$cells, !!sname := "sample"),
  by = c("protein_Id", sname)
)
stopifnot(
  nrow(check) == nrow(weights$cells),
  max(abs(check[[proteins$response()]] - check$estimate)) < 1e-8,
  all(check[[proteins$get_config()$nr_children]] == check$n)
)
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.table(
  cbind(dataset = args$dataset, weights$priors),
  file.path(args$out_dir, "priors.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

rows <- list()
for (method in methods) {
  spec <- read_spec(args$config, args$dataset, method)
  for (variant in variants) {
    lfq <- if (variant == "default") proteins else with_weight(proteins, weights, variant)
    weights_arg <- if (variant == "default") list() else list(weights = "consistency_weight")
    analysis <- suppressWarnings(suppressMessages(do.call(
      prolfqua::build_contrast_analysis,
      c(
        list(
          lfqdata = lfq,
          modelstr = paste(spec$dataset$formula, spec$method$formula_terms),
          contrasts = spec$contrasts,
          method = spec$method$facade
        ),
        spec$method$args,
        weights_arg
      )
    )))
    dir <- file.path(args$out_dir, "contrasts", method, variant)
    write_result(
      analysis$get_contrasts(),
      dir,
      spec,
      args$dataset,
      method,
      software_version = paste0("prolfqua ", utils::packageVersion("prolfqua"))
    )
    bench <- prolfquabenchmark::benchmark_from_file(dir)
    scores <- dplyr::filter(
      bench$to_summary_table(dataset = args$dataset),
      .data$contrast == "all",
      .data$score == "scaled.p_value"
    )
    fdr <- at_fdr_5(bench)
    rows[[length(rows) + 1]] <- data.frame(
      dataset = args$dataset,
      method = method,
      variant = variant,
      AUC = scores$AUC,
      pAUC_10 = scores$pAUC_10,
      FDP_5 = fdr[["FDP_5"]],
      TP_5 = fdr[["TP_5"]],
      tool_version = paste0("prolfqua ", utils::packageVersion("prolfqua"))
    )
    message(args$dataset, " ", method, " ", variant)
  }
}
utils::write.table(
  dplyr::bind_rows(rows),
  file.path(args$out_dir, "scores.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
