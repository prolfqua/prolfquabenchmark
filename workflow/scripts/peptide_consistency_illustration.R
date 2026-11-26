# Tables for the illustration in the vignette Peptide_consistency_weights, on one dataset: the peptides, residuals
# and weights of four example proteins, the variances and weight scales of every protein, the raw and shrunk variance
# ratios of 10000 random protein x sample cells, and the lm test of one contrast under four weight variants (read
# from the peptide_consistency results). Example proteins, among those with at least 4 peptides observed in every
# sample, from the lm test of the contrast: the changed protein whose -log10 p-value rises most from count to
# sample_rel, the unchanged one with p < 0.01 under count whose p-value rises most, and the proteins with the smallest
# and largest protein variance. Writes peptide_consistency_{example_peptides,example_weights,proteins,cells,tests}.tsv
# to <out_dir>.
#   Rscript workflow/scripts/peptide_consistency_illustration.R <config> <dataset> <prepared_dir> \
#     <consistency_dir> <contrast> <out_dir>
source("workflow/scripts/common.R")
source("workflow/scripts/peptide_consistency.R")
args <- script_args(c("config", "dataset", "prepared_dir", "consistency_dir", "contrast", "out_dir"))
variants <- c("count", "protein_ivw", "sample_ivw", "sample_rel")
spec <- read_spec(args$config, args$dataset, "prolfqua_lm")
changed_pattern <- spec$dataset$ground_truth$positive$pattern

write_tsv <- function(data, name) {
  utils::write.table(data, file.path(args$out_dir, name), sep = "\t", row.names = FALSE, quote = FALSE)
}

peptides <- readRDS(file.path(args$prepared_dir, "peptides.rds"))
proteins <- readRDS(file.path(args$prepared_dir, "proteins.rds"))
weights <- consistency(peptides)
weighted <- with_weight_columns(proteins, weights)
sname <- weighted$config$sample_name
factor_col <- weighted$config$factor_keys()[1]

# short sample labels: the factor level and the run's rank within it (a1, a2, ...)
samples <- unique(weighted$data[, c(sname, factor_col)])
samples <- samples[order(samples[[factor_col]], samples[[sname]]), ]
samples$sample_label <- paste0(
  samples[[factor_col]],
  stats::ave(seq_len(nrow(samples)), samples[[factor_col]], FUN = seq_along)
)
samples <- dplyr::rename(samples, sample = !!sname, condition = !!factor_col)

# lm test of the contrast under each variant, one row per protein
tests <- dplyr::bind_rows(lapply(variants, function(variant) {
  contrasts <- utils::read.delim(
    file.path(args$consistency_dir, "contrasts", "prolfqua_lm", variant, "contrasts.tsv"),
    stringsAsFactors = FALSE
  )
  contrasts <- contrasts[contrasts$contrast == args$contrast, c("protein_id", "log_fc", "t_statistic", "p_value")]
  cbind(variant = variant, contrasts)
})) |>
  tidyr::pivot_wider(names_from = "variant", values_from = c("log_fc", "t_statistic", "p_value"))
tests$changed <- grepl(changed_pattern, tests$protein_id)

# per protein: variances, consistency spread and the scale of each variant, mean weight over mean peptide count
cell_weights <- dplyr::filter(weighted$data, !is.na(.data[[proteins$response()]]))
scales <- dplyr::summarise(
  cell_weights,
  dplyr::across(
    dplyr::all_of(paste0("weight_", variants)),
    ~ mean(.x) / mean(.data$weight_count),
    .names = "scale_{sub('weight_', '', .col)}"
  ),
  factor_sd_log2 = stats::sd(log2(.data$consistency_factor)),
  .by = "protein_Id"
)
protein_table <- weights$protein |>
  dplyr::select("protein_Id", "n_peptides", "n_samples", "d", "s2", "sigma2") |>
  dplyr::left_join(scales, by = "protein_Id") |>
  dplyr::mutate(changed = grepl(changed_pattern, .data$protein_Id))

# example proteins
max_samples <- max(protein_table$n_samples)
candidates <- dplyr::filter(protein_table, .data$n_peptides >= 4, .data$n_samples == max_samples)
cand_tests <- dplyr::inner_join(tests, candidates, by = c(protein_id = "protein_Id", "changed")) |>
  dplyr::mutate(gain = log10(.data$p_value_count) - log10(.data$p_value_sample_rel))
helped <- dplyr::slice_max(dplyr::filter(cand_tests, .data$changed), .data$gain, n = 1)$protein_id
removed <- dplyr::slice_min(
  dplyr::filter(cand_tests, !.data$changed, .data$p_value_count < 0.01),
  .data$gain,
  n = 1
)$protein_id
others <- dplyr::filter(candidates, !.data$protein_Id %in% c(helped, removed))
consistent <- dplyr::slice_min(others, .data$sigma2, n = 1)$protein_Id
inconsistent <- dplyr::slice_max(others, .data$sigma2, n = 1)$protein_Id
examples <- data.frame(
  protein_Id = c(helped, removed, consistent, inconsistent),
  role = c(
    "changed, p-value falls with sample_rel",
    "unchanged, p-value rises with sample_rel",
    "most consistent peptides",
    "least consistent peptides"
  )
)
stopifnot(nrow(examples) == 4, !anyDuplicated(examples$protein_Id))

example_peptides <- weights$observations |>
  dplyr::filter(.data$protein_Id %in% examples$protein_Id) |>
  dplyr::left_join(examples, by = "protein_Id") |>
  dplyr::left_join(samples, by = "sample") |>
  dplyr::select("protein_Id", "role", "peptide_Id", "sample_label", "condition", "y", "residual", "estimate")
example_weights <- cell_weights |>
  dplyr::filter(.data$protein_Id %in% examples$protein_Id) |>
  dplyr::rename(sample = !!sname) |>
  dplyr::left_join(examples, by = "protein_Id") |>
  dplyr::left_join(samples[, c("sample", "sample_label")], by = "sample") |>
  dplyr::select(
    "protein_Id",
    "role",
    "sample_label",
    "consistency_sigma2",
    "consistency_factor",
    "consistency_factor_rel",
    dplyr::all_of(paste0("weight_", variants))
  )

set.seed(1)
cells <- weights$cells[sample(nrow(weights$cells), min(10000, nrow(weights$cells))), ] |>
  dplyr::select("n", "d_cell", "ratio_raw", "ratio", "factor_rel")

write_tsv(example_peptides, "peptide_consistency_example_peptides.tsv")
write_tsv(example_weights, "peptide_consistency_example_weights.tsv")
write_tsv(protein_table, "peptide_consistency_proteins.tsv")
write_tsv(cells, "peptide_consistency_cells.tsv")
write_tsv(cbind(dataset = args$dataset, contrast = args$contrast, tests), "peptide_consistency_tests.tsv")
