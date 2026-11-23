# Score history of the benchmark vignettes across prolfqua releases.
# `make score-history` renders the benchmark vignettes, which write their results to inst/Benchresults, and
# inst/scripts/update_score_history.R appends the scores to inst/extdata/score_history.tsv.

.score_history_columns <- c(
  "prolfqua_version",
  "date",
  "dataset",
  "method",
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


#' Score one contrast result directory for the score history
#' @param dir directory with contrasts.tsv and metadata.yaml, as written by [write_contrast_results()]
#' @return data.frame, one row per contrast x score x metric
#' @noRd
.score_rows_from_dir <- function(dir) {
  metadata <- yaml::read_yaml(file.path(dir, "metadata.yaml"))
  benchmark_from_file(dir)$to_summary_table(dataset = metadata$dataset) |>
    tidyr::pivot_longer(
      dplyr::all_of(.summary_table_metrics),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::transmute(
      dataset = .data$dataset,
      method = .data$model_name,
      tool_version = metadata$software_version,
      contrast = .data$contrast,
      score = .data$score,
      metric = .data$metric,
      value = as.numeric(.data$value)
    )
}


#' Score rows of one benchmark run
#'
#' Scores every contrast result directory (with metadata.yaml) in `results_dir` and reads every table already in
#' the history format (history_rows.tsv: null-split calibration, step timings). Aborts if a result file is older
#' than `since`: its vignette failed or skipped it, and the file holds scores of an earlier run.
#' @param results_dir directory with one subdirectory per result, as written by the vignettes
#' @param since time the run started
#' @param prolfqua_version prolfqua version recorded with the rows
#' @param date date recorded with the rows: the day of the run, or the commit date of a rerun older version
#' @return data.frame with the score history columns; values rounded to 6 significant digits
#' @noRd
.score_history_rows <- function(
  results_dir,
  since,
  prolfqua_version = as.character(utils::packageVersion("prolfqua")),
  date = format(Sys.Date())
) {
  dirs <- list.dirs(results_dir, recursive = FALSE)
  contrast_dirs <- dirs[file.exists(file.path(dirs, "metadata.yaml"))]
  row_files <- file.path(dirs, "history_rows.tsv")
  row_files <- row_files[file.exists(row_files)]
  result_files <- c(file.path(contrast_dirs, "contrasts.tsv"), row_files)
  if (length(result_files) == 0) {
    stop("No benchmark results in ", results_dir)
  }
  written <- file.mtime(result_files)
  stale <- is.na(written) | written < since
  if (any(stale)) {
    stop(
      "Results not written in this run (render failed or skipped them); ",
      "delete them to record the run without them:\n",
      paste0("  ", result_files[stale], " (", format(written[stale], "%Y-%m-%d %H:%M"), ")", collapse = "\n")
    )
  }
  rows <- dplyr::bind_rows(c(
    lapply(contrast_dirs, .score_rows_from_dir),
    lapply(row_files, utils::read.delim, stringsAsFactors = FALSE)
  ))
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
    colClasses = c(prolfqua_version = "character", date = "character", tool_version = "character")
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
