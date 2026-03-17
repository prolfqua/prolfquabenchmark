test_that("write and read contrast results round-trips correctly", {
  # Create test data with prolfqua-native column names
  test_data <- data.frame(
    protein_Id = c("P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_ECOLI",
                   "P1_HUMAN", "P2_ECOLI", "P3_HUMAN", "P4_ECOLI"),
    contrast = rep(c("A_vs_B", "C_vs_D"), each = 4),
    diff = c(0.1, 1.2, -0.3, 0.8, 0.2, 0.9, -0.1, 0.7),
    statistic = c(0.5, 3.1, -1.0, 2.5, 0.8, 2.8, -0.3, 2.2),
    p.value = c(0.6, 0.01, 0.3, 0.02, 0.4, 0.02, 0.7, 0.03),
    FDR = c(0.8, 0.05, 0.5, 0.08, 0.7, 0.06, 0.9, 0.1),
    avgInt = c(10, 12, 11, 13, 10, 12, 11, 13),
    stringsAsFactors = FALSE
  )

  metadata <- list(
    dataset = "test_dataset",
    method = "test_method",
    method_description = "Test method for round-trip",
    input_level = "protein",
    aggregation = "none",
    normalization = "none",
    software_version = "test",
    date = "2026-03-17",
    contrasts = c("A_vs_B", "C_vs_D"),
    ground_truth = list(
      id_column = "protein_id",
      positive = list(label = "ECOLI", pattern = "ECOLI"),
      negative = list(label = "HUMAN", pattern = "HUMAN")
    )
  )

  tmpdir <- tempfile("contrast_io_test")

  # Write
  write_contrast_results(test_data, tmpdir, metadata)

  expect_true(file.exists(file.path(tmpdir, "contrasts.tsv")))
  expect_true(file.exists(file.path(tmpdir, "metadata.yaml")))

  # Read back
  res <- read_contrast_results(tmpdir)

  # Check structure
  expect_true(is.data.frame(res$data))
  expect_true(is.list(res$metadata))

  # Check column renaming happened
  expect_true("log_fc" %in% colnames(res$data))
  expect_true("t_statistic" %in% colnames(res$data))
  expect_true("p_value" %in% colnames(res$data))
  expect_true("p_value_adjusted" %in% colnames(res$data))
  expect_true("avg_intensity" %in% colnames(res$data))
  expect_true("protein_id" %in% colnames(res$data))

  # Check ground truth annotation
  expect_true("species" %in% colnames(res$data))
  expect_true("TP" %in% colnames(res$data))
  expect_true(all(res$data$species %in% c("HUMAN", "ECOLI")))
  expect_equal(sum(res$data$TP), 4)  # 4 ECOLI entries

  # Check metadata preserved
  expect_equal(res$metadata$method, "test_method")
  expect_equal(res$metadata$dataset, "test_dataset")

  # Cleanup
  unlink(tmpdir, recursive = TRUE)
})


test_that("read_contrast_results errors on missing files", {
  expect_error(read_contrast_results(tempdir()), "contrasts.tsv not found")
})


test_that("ms_bench_ap computes correct average precision", {
  # Perfect classifier: precision = 1.0 at all recall levels
  # Note: strict < threshold excludes the boundary point, same as ms_bench_auc
  recall <- c(0, 0.25, 0.5, 0.75, 1.0)
  precision <- c(1.0, 1.0, 1.0, 1.0, 1.0)
  ap <- ms_bench_ap(recall, precision)
  # With < 1.0 threshold, last point excluded: integrates 0..0.75, area = 0.75, scaled = 75
  expect_equal(ap, 75, tolerance = 0.1)

  # Known trapezoid: recall 0-0.5, precision 1.0 to 0.75
  # With < 1.0 threshold, points at recall = {0, 0.5} kept, area = 0.5*0.5*(1.0+0.75) = 0.4375
  recall2 <- c(0, 0.5, 1.0)
  precision2 <- c(1.0, 0.75, 0.5)
  ap2 <- ms_bench_ap(recall2, precision2)
  expect_equal(ap2, 43.75, tolerance = 0.1)

  # Partial AP at recall_threshold = 0.6 (includes points at 0 and 0.5)
  pap <- ms_bench_ap(recall2, precision2, 0.6)
  # Area = 0.5*0.5*(1.0+0.75) = 0.4375, scaled by /0.6*100 = 72.9
  expect_equal(pap, 72.9, tolerance = 0.1)

  # Handles unsorted input
  recall3 <- c(0.5, 0, 1.0)
  precision3 <- c(0.8, 1.0, 0.6)
  ap_unsorted <- ms_bench_ap(recall3, precision3)
  ap_sorted <- ms_bench_ap(sort(recall3), c(1.0, 0.8, 0.6))
  expect_equal(ap_unsorted, ap_sorted)

  # Return type and range
  expect_true(is.numeric(ap))
  expect_true(length(ap) == 1)
  expect_true(ap >= 0 && ap <= 100)
})


test_that("write_contrast_results errors on missing required columns", {
  bad_data <- data.frame(protein_Id = "P1", contrast = "A", x = 1)
  metadata <- list(dataset = "test")
  tmpdir <- tempfile("bad_data_test")
  expect_error(
    write_contrast_results(bad_data, tmpdir, metadata),
    "Missing required columns"
  )
})
