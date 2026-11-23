# Score history of the benchmarks across prolfqua releases. The benchmark pipeline (workflow/Snakefile) scores
# every dataset and method; its `score_history` rule appends the scores and the step timings of the installed
# prolfqua to inst/extdata/score_history.tsv (workflow/scripts/update_score_history.R).

.score_history_columns <- c(
  "prolfqua_version",
  "date",
  "dataset",
  "method",
  "level",
  "tool_version",
  "contrast",
  "score",
  "metric",
  "value"
)

.score_history_keys <- c("prolfqua_version", "dataset", "method", "contrast", "score", "metric")

.summary_table_metrics <- c(
  "AUC",
  "pAUC_10",
  "pAUC_20",
  "AP",
  "pAP_50",
  "pAP_80",
  "n_TP",
  "n_TN",
  "n_total",
  "n_missing_contrasts"
)


#' Score history rows of one prolfqua version
#'
#' Reads the score tables of the pipeline (`Benchmark$to_summary_table()` plus `level` and `tool_version`, one per
#' dataset and method) and tables already in the history format (step timings).
#' @param score_files score tables written by workflow/scripts/score.R
#' @param row_files tables with the history columns except `prolfqua_version`, `date` and `level`
#' @param prolfqua_version prolfqua version recorded with the rows
#' @param date date recorded with the rows: the day of the run, or the commit date of a rerun older version
#' @return data.frame with the score history columns; values rounded to 6 significant digits
#' @noRd
.score_history_rows <- function(score_files, row_files, prolfqua_version, date = format(Sys.Date())) {
  read_table <- function(path, character_columns) {
    utils::read.delim(
      path,
      stringsAsFactors = FALSE,
      colClasses = stats::setNames(rep("character", length(character_columns)), character_columns)
    )
  }
  scores <- dplyr::bind_rows(lapply(score_files, read_table, c("tool_version", "level"))) |>
    tidyr::pivot_longer(dplyr::all_of(.summary_table_metrics), names_to = "metric", values_to = "value") |>
    dplyr::transmute(
      dataset = .data$dataset,
      method = .data$model_name,
      level = .data$level,
      tool_version = .data$tool_version,
      contrast = .data$contrast,
      score = .data$score,
      metric = .data$metric,
      value = as.numeric(.data$value)
    )
  rows <- dplyr::bind_rows(c(list(scores), lapply(row_files, read_table, "tool_version")))
  rows$prolfqua_version <- prolfqua_version
  rows$date <- date
  rows$value <- signif(rows$value, 6)
  as.data.frame(rows)[, .score_history_columns]
}


#' Read the score history
#' @param path score history TSV
#' @return data.frame with the score history columns
#' @noRd
.read_score_history <- function(path) {
  utils::read.delim(
    path,
    stringsAsFactors = FALSE,
    colClasses = c(prolfqua_version = "character", date = "character", tool_version = "character", level = "character")
  )
}


#' Append the rows of one prolfqua version to the score history
#' @param path score history TSV; created if it does not exist
#' @param rows data.frame as returned by `.score_history_rows()`
#' @param replace replace the rows of a version already in the history; otherwise such a version aborts
#' @return the new history, invisibly
#' @noRd
.append_score_history <- function(path, rows, replace = FALSE) {
  history <- if (file.exists(path)) .read_score_history(path) else rows[0, .score_history_columns]
  versions <- unique(rows$prolfqua_version)
  recorded <- intersect(versions, history$prolfqua_version)
  if (length(recorded) > 0 && !replace) {
    stop(
      "prolfqua ",
      paste(recorded, collapse = ", "),
      " is already in ",
      path,
      "; rerun with replace = TRUE (make score-history REPLACE=1) to replace its rows"
    )
  }
  rows <- rows[do.call(order, unname(as.list(rows[, .score_history_keys]))), .score_history_columns]
  history <- rbind(history[!history$prolfqua_version %in% versions, ], rows)
  duplicated_keys <- duplicated(history[, .score_history_keys])
  if (any(duplicated_keys)) {
    stop(
      "Rows with the same ",
      paste(.score_history_keys, collapse = ", "),
      ":\n",
      paste(utils::capture.output(print(utils::head(history[duplicated_keys, .score_history_keys]))), collapse = "\n")
    )
  }
  utils::write.table(history, path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(history)
}
