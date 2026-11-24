# Fits proDA to the protein table of the search software (input: protein_table), as prepared for the dataset.
#   Rscript workflow/scripts/fit_proda.R <config> <dataset> <method> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "method", "prepared_dir", "out_dir"))
spec <- read_spec(args$config, args$dataset, args$method)

experiment <- prolfqua::lfqdata_to_summarized_experiment(readRDS(file.path(args$prepared_dir, "protein_table.rds")))
# proDA's fit draws random numbers (two unseeded fits of the same data differ by up to 1.3 in log fold change)
set.seed(1)
fit <- proDA::proDA(
  experiment,
  design = stats::as.formula(paste(spec$dataset$formula, "- 1")),
  data_is_log_transformed = TRUE
)
result <- dplyr::bind_rows(lapply(names(spec$contrasts), function(name) {
  data.frame(contrast = name, proDA::test_diff(fit, contrast = spec$contrasts[[name]]))
}))
result <- dplyr::rename(
  result,
  protein_id = "name",
  log_fc = "diff",
  p_value = "pval",
  p_value_adjusted = "adj_pval",
  avg_intensity = "avg_abundance"
)
write_result(
  result,
  args$out_dir,
  spec,
  args$dataset,
  args$method,
  software_version = paste0("proDA ", utils::packageVersion("proDA"))
)
