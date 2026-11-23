write_fixture_result <- function(results_dir, method) {
  data <- data.frame(
    protein_id = rep(c("P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_ECOLI"), 2),
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 4),
    log_fc = c(0.1, 1.2, -0.3, 0.8, 0.2, 0.9, -0.1, 0.7),
    t_statistic = c(0.5, 3.1, -1.0, 2.5, 0.8, 2.8, -0.3, 2.2),
    p_value = c(0.6, 0.01, 0.3, 0.02, 0.4, 0.02, 0.7, 0.03),
    p_value_adjusted = c(0.8, 0.05, 0.5, 0.08, 0.7, 0.06, 0.9, 0.1),
    avg_intensity = c(10, 12, 11, 13, 10, 12, 11, 13)
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
  suppressMessages(write_contrast_results(data, file.path(results_dir, method), metadata))
}

write_fixture_calibration <- function(results_dir) {
  dir.create(file.path(results_dir, "null_split"))
  calibration <- data.frame(
    dataset = "test_null_split",
    method = "none / lm",
    tool_version = "prolfqua 1.9.5",
    contrast = c("random_median", "run_block"),
    score = "p_value",
    metric = "p_0.05",
    value = c(0.05, 0.1)
  )
  utils::write.table(
    calibration,
    file.path(results_dir, "null_split", "history_rows.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

fixture_results <- function() {
  results_dir <- tempfile("score_history_test")
  dir.create(results_dir)
  write_fixture_result(results_dir, "method_a")
  write_fixture_result(results_dir, "method_b")
  write_fixture_calibration(results_dir)
  results_dir
}


test_that(".score_history_rows scores every result directory and reads the history_rows tables", {
  since <- Sys.time() - 60
  results_dir <- fixture_results()

  rows <- .score_history_rows(results_dir, since = since, prolfqua_version = "1.4.0", date = "2025-11-09")

  expect_equal(colnames(rows), .score_history_columns)
  expect_false(anyDuplicated(rows[, .score_history_keys]) > 0)
  expect_equal(unique(rows$prolfqua_version), "1.4.0")
  expect_equal(unique(rows$date), "2025-11-09")
  expect_setequal(unique(rows$method), c("method_a", "method_b", "none / lm"))
  expect_true(all(c("AUC", "pAUC_10", "AP", "n_missing_contrasts", "p_0.05") %in% rows$metric))
  expect_equal(unique(rows$tool_version[rows$method == "method_a"]), "tool 1.0")
  expect_type(rows$value, "double")

  unlink(results_dir, recursive = TRUE)
})


test_that(".score_history_rows aborts on results older than the run, naming them", {
  results_dir <- fixture_results()
  Sys.setFileTime(file.path(results_dir, "method_b", "contrasts.tsv"), as.POSIXct("2020-01-01"))

  expect_error(
    .score_history_rows(results_dir, since = Sys.time() - 60, prolfqua_version = "1.9.5"),
    "method_b/contrasts.tsv \\(2020-01-01"
  )

  unlink(results_dir, recursive = TRUE)
})


test_that(".append_score_history refuses a recorded version unless replace = TRUE", {
  since <- Sys.time() - 60
  results_dir <- fixture_results()
  history_file <- tempfile(fileext = ".tsv")
  rows <- .score_history_rows(results_dir, since = since, prolfqua_version = "1.9.5")

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
  expect_equal(.read_score_history(history_file)$prolfqua_version, history$prolfqua_version)

  unlink(c(results_dir, history_file), recursive = TRUE)
})
