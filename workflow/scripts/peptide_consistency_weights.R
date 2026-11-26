# Peptide consistency weights (vignette Peptide_consistency_weights): weights of the protein x sample abundances from
# how well a protein's peptides agree, the residuals of the median polish that aggregated them. Recomputes the median
# polish of each protein on the prepared peptides (checked against the prepared proteins), adds every weight variant as
# a column (peptide_consistency.R), fits the protein-level prolfqua facades lm, limma, deqms and rlm with each, and
# scores the fits like score.R (score_fit()): AUC, pAUC_5 and pAUC_10 of the signed p-value, all contrasts pooled, and
# the false discovery proportion and number of true discoveries at a nominal FDR of 5%. Writes <out_dir>/scores.tsv (one
# row per facade and variant) and <out_dir>/priors.tsv (the empirical Bayes priors of the dataset), the contrasts under
# <out_dir>/contrasts/<method>/<variant>.
#   Rscript workflow/scripts/peptide_consistency_weights.R <config> <dataset> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
source("workflow/scripts/peptide_consistency.R")
args <- script_args(c("config", "dataset", "prepared_dir", "out_dir"))
methods <- c("prolfqua_lm", "prolfqua_limma", "prolfqua_deqms", "prolfqua_rlm")
# default: the prepared proteins as the main benchmark fits them, with the facade's own weights
variants <- c("default", weight_variants)

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
weighted <- with_weight_columns(proteins, weights)

rows <- list()
for (method in methods) {
  spec <- read_spec(args$config, args$dataset, method)
  for (variant in variants) {
    lfq <- if (variant == "default") proteins else lfq_with_weight(weighted, variant)
    weights_arg <- if (variant == "default") list() else list(weights = paste0("weight_", variant))
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
    rows[[length(rows) + 1]] <- score_fit(dir, args$dataset, method, variant)
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
