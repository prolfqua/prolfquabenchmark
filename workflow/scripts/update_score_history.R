# Appends the scores of every dataset and method and the step timings of one prolfqua version to the score history,
# dated <date> (`today`, or the commit date of a release). Refuses a prolfqua version already recorded unless
# <replace> is TRUE (make score-history REPLACE=1).
#   Rscript workflow/scripts/update_score_history.R <history.tsv> <timing_rows.tsv> <replace> <version> <date> \
#     <score files...>
source("workflow/scripts/common.R")
args <- script_args(c("history", "timing", "replace", "version", "date"), rest = "score_files")

rows <- prolfquabenchmark:::.score_history_rows(
  args$score_files,
  args$timing,
  prolfqua_version = args$version,
  date = if (args$date == "today") format(Sys.Date()) else args$date
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
