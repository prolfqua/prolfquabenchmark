#!/usr/bin/env Rscript
# Appends the scores of the benchmark vignettes and the step timings to inst/extdata/score_history.tsv, one block of
# rows per prolfqua version (vignette `Score_history`). Run from the package root, through the Makefile, which
# renders the benchmark vignettes and times the prolfqua steps first (results in inst/Benchresults):
#   make score-history            # a prolfqua version not yet in the history
#   make score-history REPLACE=1  # rerun a version already in the history
# Only results written after <results-dir>/.score_run (touched at the start of the run) are taken.
# Options for reruns of older prolfqua versions (inst/scripts/prefill_score_history.sh):
#   --results-dir=<dir>  results of that run (default inst/Benchresults)
#   --version=<x.y.z>    prolfqua version of that run (default: the installed prolfqua)
#   --date=<yyyy-mm-dd>  date recorded with the rows (default: today)

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  value <- sub(paste0("^--", name, "="), "", grep(paste0("^--", name, "="), args, value = TRUE))
  if (length(value) == 0) default else value
}
replace <- "--replace" %in% args
results_dir <- option("results-dir", "inst/Benchresults")

devtools::load_all(quiet = TRUE)

since <- file.mtime(file.path(results_dir, ".score_run"))
if (is.na(since)) {
  stop(file.path(results_dir, ".score_run"), " not found; run `make score-history`")
}
rows <- .score_history_rows(
  results_dir,
  since = since,
  prolfqua_version = option("version", as.character(utils::packageVersion("prolfqua"))),
  date = option("date", format(Sys.Date()))
)
history <- .append_score_history("inst/extdata/score_history.tsv", rows, replace = replace)

message(
  "prolfqua ",
  unique(rows$prolfqua_version),
  ": ",
  nrow(rows),
  " rows added; ",
  length(unique(history$prolfqua_version)),
  " versions in inst/extdata/score_history.tsv"
)
print(as.data.frame(table(dataset = rows$dataset, method = rows$method)) |> subset(Freq > 0), row.names = FALSE)
