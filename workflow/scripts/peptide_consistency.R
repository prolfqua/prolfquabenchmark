# Peptide consistency weights (vignette Peptide_consistency_weights), shared by peptide_consistency_weights.R and
# peptide_consistency_illustration.R: the median polish residuals of each protein's peptides, their variances per
# protein and per sample with empirical Bayes shrinkage, and every weight variant as a protein x sample column.

# The weight variants, each a column weight_<variant> of the protein table
weight_variants <- c(
  "count",
  "count_conf",
  "protein_ivw",
  "sample_ivw",
  "sample_ivw_conf",
  "sample_rel",
  "sample_rel_conf"
)

# Median polish of each protein on its observed peptides, as prolfqua's medpolish aggregation; one row per
# observation with its intensity y, residual and the protein estimate of its sample
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
      peptide_Id = d$peptide_Id,
      sample = d[[sname]],
      y = d[[response]],
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
# Returns the observations, the cells (one row per protein and sample), the proteins and the priors.
consistency <- function(peptides) {
  observations <- polish(peptides)
  cells <- dplyr::summarise(
    observations,
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
  protein$s2 <- ifelse(protein$d > 0, protein$rss / protein$d, NA)
  squeezed_protein <- limma::squeezeVar(ifelse(protein$d > 0, protein$s2, 0), protein$d)
  protein$sigma2 <- squeezed_protein$var.post
  cells <- dplyr::left_join(cells, protein[, c("protein_Id", "N", "d", "sigma2")], by = "protein_Id")
  # the protein's df split over its samples in proportion to n - 1
  cells$d_cell <- ifelse(cells$N > cells$n_samples, (cells$n - 1) * cells$d / (cells$N - cells$n_samples), 0)
  cells$ratio_raw <- ifelse(cells$d_cell > 0, cells$rss / cells$d_cell / cells$sigma2, NA)
  squeezed_ratio <- limma::squeezeVar(ifelse(cells$d_cell > 0, cells$ratio_raw, 0), cells$d_cell)
  cells$ratio <- squeezed_ratio$var.post
  cells$factor <- 1 / cells$ratio
  cells <- dplyr::mutate(cells, factor_rel = .data$factor / mean(.data$factor), .by = "protein_Id")
  list(
    observations = observations,
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

# The prepared proteins with the parts of the weights (consistency_sigma2, consistency_factor,
# consistency_factor_rel) and every weight variant (weight_<variant>) as protein x sample columns, as the
# aggregation would report them; rows without observed peptides (no response) have factor 1. The configuration is
# rebuilt from its values, so prepared files serialized by an older prolfqua take today's fields, and fits use the
# identification confidence only through the variants that name it.
with_weight_columns <- function(proteins, weights) {
  config <- prolfqua::list_to_analysis_configuration(prolfqua::r6_extract_values(proteins$get_config()))
  config$use_identification_confidence <- FALSE
  cells <- dplyr::rename(
    weights$cells[, c("protein_Id", "sample", "factor", "factor_rel")],
    !!config$sample_name := "sample",
    consistency_factor = "factor",
    consistency_factor_rel = "factor_rel"
  )
  data <- proteins$data_long() |>
    dplyr::left_join(cells, by = c("protein_Id", config$sample_name)) |>
    dplyr::left_join(
      dplyr::rename(weights$protein[, c("protein_Id", "sigma2")], consistency_sigma2 = "sigma2"),
      by = "protein_Id"
    ) |>
    dplyr::mutate(
      consistency_factor = dplyr::coalesce(.data$consistency_factor, 1),
      consistency_factor_rel = dplyr::coalesce(.data$consistency_factor_rel, 1)
    )
  n <- data[[config$nr_children]]
  confidence <- data[[config$identification_confidence]]
  sample_ivw <- n * data$consistency_factor / data$consistency_sigma2
  sample_rel <- n * data$consistency_factor_rel
  data$weight_count <- n
  data$weight_count_conf <- n * confidence
  data$weight_protein_ivw <- n / data$consistency_sigma2
  data$weight_sample_ivw <- sample_ivw
  data$weight_sample_ivw_conf <- sample_ivw * confidence
  data$weight_sample_rel <- sample_rel
  data$weight_sample_rel_conf <- sample_rel * confidence
  list(data = data, config = config)
}

# LFQData of the proteins fitted with the weight of `variant`: its column is named as nr_children, so that limma
# takes it as a feature x sample value column; nr_children_feature (the DEqMS count) stays the peptide count
lfq_with_weight <- function(weighted, variant) {
  config <- weighted$config$clone(deep = TRUE)
  config$nr_children <- paste0("weight_", variant)
  prolfqua::LFQData$new(weighted$data, config)
}
