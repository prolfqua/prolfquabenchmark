# Appends the scores of every dataset and method and the step timings of the installed prolfqua to the score history.
# Refuses a prolfqua version already recorded unless <replace> is TRUE (make score-history REPLACE=1).
#   Rscript workflow/scripts/update_score_history.R <config> <history.tsv> <timing_rows.tsv> <replace>
source("workflow/scripts/common.R")
args <- script_args(c("config", "history", "timing", "replace"))
config <- yaml::read_yaml(args$config)

score_files <- unlist(lapply(names(config$methods), function(method) {
  file.path("results/scores", unlist(config$methods[[method]]$datasets), paste0(method, ".tsv"))
}))
rows <- prolfquabenchmark:::.score_history_rows(
  score_files,
  args$timing,
  prolfqua_version = as.character(utils::packageVersion("prolfqua"))
)
history <- prolfquabenchmark:::.append_score_history(args$history, rows, replace = as.logical(args$replace))
message(
  "prolfqua ",
  unique(rows$prolfqua_version),
  ": ",
  nrow(rows),
  " rows; ",
  length(unique(history$prolfqua_version)),
  " versions in ",
  args$history
)
