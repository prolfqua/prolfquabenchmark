test_that("benchmark_from_file creates a working Benchmark object", {
  # Create test data in standard format (already snake_case)
  test_data <- data.frame(
    protein_id = c("P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_ECOLI",
                   "P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_ECOLI"),
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 4),
    log_fc = c(0.1, 1.2, -0.3, 0.8, 0.2, 0.9, -0.1, 0.7),
    t_statistic = c(0.5, 3.1, -1.0, 2.5, 0.8, 2.8, -0.3, 2.2),
    p_value = c(0.6, 0.01, 0.3, 0.02, 0.4, 0.02, 0.7, 0.03),
    p_value_adjusted = c(0.8, 0.05, 0.5, 0.08, 0.7, 0.06, 0.9, 0.1),
    avg_intensity = c(10, 12, 11, 13, 10, 12, 11, 13),
    stringsAsFactors = FALSE
  )

  metadata <- list(
    dataset = "test_dataset",
    method = "test_lm",
    method_description = "Test linear model",
    ground_truth = list(
      id_column = "protein_id",
      positive = list(label = "ECOLI", pattern = "ECOLI"),
      negative = list(label = "HUMAN", pattern = "HUMAN")
    )
  )

  tmpdir <- tempfile("benchmark_io_test")
  dir.create(tmpdir)

  # Write files directly in standard format (no column mapping needed)
  utils::write.table(test_data, file.path(tmpdir, "contrasts.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  yaml::write_yaml(metadata, file.path(tmpdir, "metadata.yaml"))

  # Create benchmark from file
  bench <- benchmark_from_file(tmpdir)
  expect_s3_class(bench, "Benchmark")
  expect_equal(bench$model_name, "test_lm")

  # Check pAUC works (includes AP metrics)
  pauc <- bench$pAUC()
  expect_true(is.data.frame(pauc))
  expect_true("AUC" %in% colnames(pauc))
  expect_true(all(c("AP", "pAP_50", "pAP_80") %in% colnames(pauc)))

  # Check to_summary_table works
  summary <- bench$to_summary_table(dataset = "test_dataset")
  expect_true(is.data.frame(summary))
  expect_true(all(c("model_name", "dataset", "contrast", "score",
                     "AUC", "pAUC_10", "pAUC_20",
                     "AP", "pAP_50", "pAP_80") %in% colnames(summary)))
  expect_equal(unique(summary$model_name), "test_lm")
  expect_equal(unique(summary$dataset), "test_dataset")

  # Cleanup
  unlink(tmpdir, recursive = TRUE)
})


test_that("to_summary_table works on existing Benchmark data", {
  dd <- dplyr::filter(
    prolfqua::prolfqua_data("data_benchmarkExample"),
    !is.na(statistic)
  )
  dd <- dd |> dplyr::mutate(avgInt = (c1 + c2) / 2)
  ttd <- ionstar_bench_preprocess(dd)

  bench <- make_benchmark(
    ttd$data,
    benchmark = list(
      list(score = "estimate", desc = TRUE),
      list(score = "statistic", desc = TRUE),
      list(score = "scaled.p.value", desc = TRUE)
    ),
    fcestimate = "estimate",
    model_description = "test model",
    model_name = "test_model"
  )

  summary <- bench$to_summary_table(dataset = "ionstar")
  expect_true(is.data.frame(summary))
  expect_true(nrow(summary) > 0)
  expect_true(all(summary$AUC >= 0 & summary$AUC <= 100))
  expect_equal(unique(summary$dataset), "ionstar")
})


test_that("write and read benchmark_results round-trips", {
  summary <- data.frame(
    model_name = "test",
    model_description = "test model",
    dataset = "ds",
    contrast = "A_vs_B",
    score = "log_fc",
    AUC = 85.3,
    pAUC_10 = 42.1,
    pAUC_20 = 55.7,
    n_TP = 100,
    n_TN = 500,
    n_total = 600,
    n_missing_contrasts = 5,
    stringsAsFactors = FALSE
  )

  tmpfile <- tempfile(fileext = ".tsv")
  write_benchmark_results(summary, tmpfile)
  result <- read_benchmark_results(tmpfile)

  expect_equal(nrow(result), 1)
  expect_equal(result$AUC, 85.3)
  expect_equal(result$model_name, "test")

  unlink(tmpfile)
})


test_that("summary_metrics returns one row per score with valid values", {
  dd <- dplyr::filter(
    prolfqua::prolfqua_data("data_benchmarkExample"),
    !is.na(statistic)
  )
  dd <- dd |> dplyr::mutate(avgInt = (c1 + c2) / 2)
  ttd <- ionstar_bench_preprocess(dd)

  bench <- make_benchmark(
    ttd$data,
    benchmark = list(
      list(score = "estimate", desc = TRUE),
      list(score = "statistic", desc = TRUE),
      list(score = "scaled.p.value", desc = TRUE)
    ),
    fcestimate = "estimate",
    model_description = "test model",
    model_name = "test_model"
  )

  sm <- bench$summary_metrics()
  expect_true(is.data.frame(sm))
  expect_equal(nrow(sm), 3)  # one row per score
  expect_true(all(c("model_name", "score", "AUC", "pAUC_10", "pAUC_20",
                     "AP", "pAP_50", "pAP_80", "n_total") %in% colnames(sm)))
  expect_true(all(sm$AUC >= 0 & sm$AUC <= 100))
  expect_true(all(sm$AP >= 0 & sm$AP <= 100))
  expect_true(all(sm$n_total > 0))
  expect_equal(unique(sm$model_name), "test_model")

  # Geometric mean should be <= arithmetic mean
  pauc <- bench$pAUC() |> dplyr::filter(contrast != "all")
  for (s in sm$score) {
    arith_auc <- mean(pauc$AUC[pauc$what == s])
    expect_true(sm$AUC[sm$score == s] <= arith_auc + 1e-10)
  }
})


test_that("calibration_metrics returns per-contrast and average rows", {
  dd <- dplyr::filter(
    prolfqua::prolfqua_data("data_benchmarkExample"),
    !is.na(statistic)
  )
  dd <- dd |> dplyr::mutate(avgInt = (c1 + c2) / 2)
  ttd <- ionstar_bench_preprocess(dd)

  bench <- make_benchmark(
    ttd$data,
    benchmark = list(
      list(score = "estimate", desc = TRUE),
      list(score = "statistic", desc = TRUE),
      list(score = "scaled.p.value", desc = TRUE)
    ),
    fcestimate = "estimate",
    FDRvsFDP = list(list(score = "p.value", desc = FALSE)),
    model_description = "test model",
    model_name = "test_model"
  )

  cal <- bench$calibration_metrics(fdr_threshold = 0.2)
  expect_true(is.data.frame(cal))
  expect_true(all(c("model_name", "contrast", "score",
                     "FDP_cal", "FDP_bias") %in% colnames(cal)))
  expect_equal(unique(cal$model_name), "test_model")

  # Should have per-contrast rows + "average" row
  expect_true("average" %in% cal$contrast)
  n_contrasts <- length(unique(cal$contrast[cal$contrast != "average"]))
  expect_true(n_contrasts > 0)

  # FDP_cal should be non-negative
  expect_true(all(cal$FDP_cal >= 0))

  # Average should equal arithmetic mean of per-contrast values
  for (s in unique(cal$score)) {
    per_c <- cal |> dplyr::filter(score == s, contrast != "average")
    avg_row <- cal |> dplyr::filter(score == s, contrast == "average")
    expect_equal(avg_row$FDP_cal, mean(per_c$FDP_cal), tolerance = 1e-10)
    expect_equal(avg_row$FDP_bias, mean(per_c$FDP_bias), tolerance = 1e-10)
  }
})


test_that("collect_benchmark_results combines multiple files", {
  s1 <- data.frame(model_name = "m1", contrast = "A", score = "x",
                   AUC = 80, pAUC_10 = 40, pAUC_20 = 50)
  s2 <- data.frame(model_name = "m2", contrast = "A", score = "x",
                   AUC = 90, pAUC_10 = 45, pAUC_20 = 55)

  f1 <- tempfile(fileext = ".tsv")
  f2 <- tempfile(fileext = ".tsv")
  write_benchmark_results(s1, f1)
  write_benchmark_results(s2, f2)

  combined <- collect_benchmark_results(c(f1, f2))
  expect_equal(nrow(combined), 2)
  expect_equal(sort(combined$model_name), c("m1", "m2"))

  unlink(c(f1, f2))
})
