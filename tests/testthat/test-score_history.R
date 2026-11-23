write_fixture_scores <- function(dir, method) {
  data <- data.frame(
    protein_id = rep(c("P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_ECOLI"), 2),
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 4),
    log_fc = c(0.1, 1.2, -0.3, 0.8, 0.2, 0.9, -0.1, 0.7),
    t_statistic = c(0.5, 3.1, -1.0, 2.5, 0.8, 2.8, -0.3, 2.2),
    p_value = c(0.6, 0.01, 0.3, 0.02, 0.4, 0.02, 0.7, 0.03),
    p_value_adjusted = c(0.8, 0.05, 0.5, 0.08, 0.7, 0.06, 0.9, 0.1)
  )
  metadata <- list(
    dataset = "test_dataset",
    method = method,
    method_description = method,
    input_file = "test.tsv",
    software_version = "tool 1.0",
    date = "2026-11-23",
    ground_truth = list(
      id_column = "protein_id",
      positive = list(label = "ECOLI", pattern = "ECOLI"),
      negative = list(label = "HUMAN", pattern = "HUMAN")
    )
  )
  contrast_dir <- file.path(dir, method)
  suppressMessages(write_contrast_results(data, contrast_dir, metadata))
  scores <- benchmark_from_file(contrast_dir)$to_summary_table(dataset = "test_dataset")
  scores$level <- "protein"
  scores$tool_version <- "tool 1.0"
  path <- file.path(dir, paste0(method, "_scores.tsv"))
  utils::write.table(scores, path, sep = "\t", row.names = FALSE, quote = FALSE)
  path
}

write_fixture_timings <- function(dir) {
  path <- file.path(dir, "history_rows.tsv")
  utils::write.table(
    data.frame(
      dataset = "test_dataset",
      method = c("aggregate", "fit_lm"),
      tool_version = "R 4.6.2",
      contrast = "all",
      score = "elapsed",
      metric = "elapsed_s",
      value = c(5.3, 1.9)
    ),
    path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  path
}


test_that(".score_history_rows reads the score tables and the timing rows", {
  dir <- tempfile("score_history_test")
  dir.create(dir)
  score_files <- c(write_fixture_scores(dir, "method_a"), write_fixture_scores(dir, "method_b"))

  rows <- .score_history_rows(score_files, write_fixture_timings(dir), prolfqua_version = "1.4.0", date = "2025-11-09")

  expect_equal(colnames(rows), .score_history_columns)
  expect_false(anyDuplicated(rows[, .score_history_keys]) > 0)
  expect_equal(unique(rows$prolfqua_version), "1.4.0")
  expect_equal(unique(rows$date), "2025-11-09")
  expect_setequal(unique(rows$method), c("method_a", "method_b", "aggregate", "fit_lm"))
  expect_true(all(c("AUC", "pAUC_10", "AP", "n_missing_contrasts", "elapsed_s") %in% rows$metric))
  expect_equal(unique(rows$level[rows$method == "method_a"]), "protein")
  expect_true(all(is.na(rows$level[rows$metric == "elapsed_s"])))
  expect_type(rows$value, "double")

  unlink(dir, recursive = TRUE)
})


test_that(".append_score_history refuses a recorded version unless replace = TRUE", {
  dir <- tempfile("score_history_test")
  dir.create(dir)
  history_file <- file.path(dir, "history.tsv")
  rows <- .score_history_rows(write_fixture_scores(dir, "method_a"), write_fixture_timings(dir), "1.9.5")

  .append_score_history(history_file, rows)
  expect_equal(nrow(.read_score_history(history_file)), nrow(rows))

  expect_error(.append_score_history(history_file, rows), "prolfqua 1.9.5 is already in")

  rows$value <- rows$value + 1
  history <- .append_score_history(history_file, rows, replace = TRUE)
  expect_equal(nrow(history), nrow(rows))
  expect_equal(sort(history$value), sort(rows$value))

  rows$prolfqua_version <- "1.10.0"
  history <- .append_score_history(history_file, rows)
  expect_equal(nrow(history), 2 * nrow(rows))
  expect_equal(.read_score_history(history_file)$level, history$level)

  unlink(dir, recursive = TRUE)
})
