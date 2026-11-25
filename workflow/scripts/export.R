# Exports the current benchmark for the display vignettes (vignettes/Benchmark_results.Rmd): scores, curves, the
# scores on the pairs of the Pareto-front methods (results/same_pairs/<dataset>.tsv), a methods table (tool version,
# level, input, fit time), a datasets table (with the labels of the changed and unchanged proteins) and their
# preprocessing steps, as TSV files in <out_dir>. The score
# files (results/scores/<dataset>/<method>.tsv) name the jobs; the Snakefile passes them.
#   Rscript workflow/scripts/export.R <config> <out_dir> <score files...>
source("workflow/scripts/common.R")
args <- script_args(c("config", "out_dir"), rest = "score_files")
config <- yaml::read_yaml(args$config)

jobs <- data.frame(
  dataset = basename(dirname(args$score_files)),
  method = sub("\\.tsv$", "", basename(args$score_files))
)
job_file <- function(dir, extension) file.path("results", dir, jobs$dataset, paste0(jobs$method, extension))
read_tables <- function(paths, colClasses = NA) {
  dplyr::bind_rows(lapply(paths, utils::read.delim, stringsAsFactors = FALSE, colClasses = colClasses))
}
write_tsv <- function(data, name) {
  utils::write.table(data, file.path(args$out_dir, name), sep = "\t", row.names = FALSE, quote = FALSE)
}

write_tsv(read_tables(job_file("scores", ".tsv"), c(tool_version = "character")), "benchmark_scores.tsv")
write_tsv(read_tables(job_file("curves", ".tsv")), "benchmark_curves.tsv")
write_tsv(
  read_tables(file.path("results/same_pairs", paste0(sort(unique(jobs$dataset)), ".tsv"))),
  "benchmark_same_pairs.tsv"
)

methods <- dplyr::bind_rows(lapply(seq_len(nrow(jobs)), function(i) {
  metadata <- yaml::read_yaml(file.path("results/contrasts", jobs$dataset[i], jobs$method[i], "metadata.yaml"))
  fit <- utils::read.delim(file.path("results/benchmarks", jobs$dataset[i], paste0(jobs$method[i], ".tsv")))
  data.frame(
    dataset = jobs$dataset[i],
    method = jobs$method[i],
    tool = config$methods[[jobs$method[i]]]$tool,
    tool_version = metadata$software_version,
    level = metadata$level,
    input = metadata$input,
    description = metadata$method_description,
    fit_seconds = round(fit$s, 1),
    fit_max_rss_mb = round(fit$max_rss)
  )
}))
write_tsv(methods, "benchmark_methods.tsv")

write_tsv(
  data.frame(
    dataset = names(config$datasets),
    input_file = vapply(config$datasets, function(d) d$input_file, ""),
    description = vapply(config$datasets, function(d) d$description, ""),
    positive = vapply(config$datasets, function(d) d$ground_truth$positive$label, ""),
    negative = vapply(config$datasets, function(d) d$ground_truth$negative$label, "")
  ),
  "benchmark_datasets.tsv"
)

steps <- lapply(config$datasets, function(d) unlist(d$preprocessing))
write_tsv(
  data.frame(dataset = rep(names(steps), lengths(steps)), step = unlist(steps, use.names = FALSE)),
  "benchmark_preprocessing.tsv"
)
