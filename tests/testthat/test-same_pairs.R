# Scoring methods on the same protein-contrast pairs

# One method's results in the contrasts.tsv format, E. coli proteins changed, human proteins unchanged; NA marks the
# pairs the method does not test
same_pairs_result <- function(method, log_fc, p_value) {
  ids <- c("P1_ECOLI", "P2_ECOLI", "P3_HUMAN", "P4_HUMAN", "P5_ECOLI", "P6_HUMAN")
  list(
    data = data.frame(
      protein_id = rep(ids, 2),
      contrast = rep(c("c1", "c2"), each = 6),
      log_fc = log_fc,
      t_statistic = log_fc * 3,
      p_value = p_value,
      p_value_adjusted = pmin(p_value * 2, 1),
      avg_intensity = 10
    ),
    metadata = list(
      dataset = "test",
      method = method,
      method_description = method,
      ground_truth = list(
        id_column = "protein_id",
        positive = list(label = "ECOLI", pattern = "ECOLI"),
        negative = list(label = "HUMAN", pattern = "HUMAN")
      )
    )
  )
}

# selective: P1-P4 only; complete: all six proteins, but P2 in c2 without a p-value; partial: P2-P5;
# equal: the pairs of selective
same_pairs_results <- function() {
  selective_fc <- c(2, 1.5, 0.1, -0.2, NA, NA, 1.8, 1.2, 0.2, 0.1, NA, NA)
  selective_p <- c(0.001, 0.01, 0.8, 0.7, NA, NA, 0.002, 0.02, 0.6, 0.9, NA, NA)
  complete_fc <- c(1.9, 1.4, 0.3, -0.1, 0.5, 1.6, 1.7, 1.1, 0.1, 0.2, 0.4, 1.5)
  complete_p <- c(0.001, 0.02, 0.6, 0.8, 0.3, 0.004, 0.003, NA, 0.7, 0.6, 0.4, 0.01)
  partial_fc <- c(NA, 1.5, 0.1, -0.2, 0.6, NA, NA, 1.2, 0.2, 0.1, 0.5, NA)
  partial_p <- c(NA, 0.01, 0.8, 0.7, 0.2, NA, NA, 0.02, 0.6, 0.9, 0.3, NA)
  list(
    selective = same_pairs_result("selective", selective_fc, selective_p),
    complete = same_pairs_result("complete", complete_fc, complete_p),
    partial = same_pairs_result("partial", partial_fc, partial_p),
    equal = same_pairs_result("equal", selective_fc * 0.9, selective_p)
  )
}

# Score table rows with a known front: selective (pAUC_5 100, 4 proteins) and complete (50, 6); partial and equal are
# dominated
same_pairs_score_table <- function() {
  summary <- data.frame(
    model_name = c("selective", "complete", "partial", "equal"),
    pAUC_5 = c(100, 50, 40, 90),
    n_total = c(4, 6, 4, 4)
  )
  rows <- merge(
    summary,
    expand.grid(
      score = c("log_fc", "t_statistic", "scaled.p_value"),
      contrast = c("all", "c1", "c2"),
      stringsAsFactors = FALSE
    )
  )
  rows
}

test_that("pareto_optimal keeps the points no other point beats in both coordinates", {
  expect_equal(pareto_optimal(c(1, 2, 3), c(3, 2, 1)), c(TRUE, TRUE, TRUE))
  expect_equal(pareto_optimal(c(1, 2, 2), c(1, 2, 1)), c(FALSE, TRUE, FALSE))
  # equal points do not beat each other
  expect_equal(pareto_optimal(c(2, 2), c(2, 2)), c(TRUE, TRUE))
})

test_that("score_on_pairs scores the rows of the given pairs only", {
  result <- same_pairs_results()$complete
  pairs <- data.frame(
    protein_id = c("P1_ECOLI", "P3_HUMAN", "P2_ECOLI", "P4_HUMAN"),
    contrast = c("c1", "c1", "c2", "c2")
  )
  filtered <- result
  filtered$data <- merge(result$data, pairs)
  roc <- benchmark_from_result(filtered)$get_confusion_benchmark()
  roc <- roc[roc$contrast == "all" & roc$what == "log_fc", ]
  expect_equal(
    score_on_pairs(result, pairs, "log_fc"),
    c(pAUC_5 = ms_bench_auc(roc$FPR, roc$TPR, 0.05), pAUC_10 = ms_bench_auc(roc$FPR, roc$TPR, 0.1))
  )
})

test_that("same_pairs_scores scores the strict supersets of a front method on its pairs", {
  results <- same_pairs_results()
  scores <- same_pairs_scores(results, same_pairs_score_table())

  # complete tests every pair of selective with a fold change, so it is scored there for log_fc and t_statistic
  log_fc <- scores[scores$score == "log_fc", ]
  expect_equal(log_fc$front_method, c("selective", "selective"))
  expect_equal(log_fc$method, c("selective", "complete"))
  expect_equal(log_fc$proteins, c(4, 4))
  selective_pairs <- results$selective$data[!is.na(results$selective$data$log_fc), c("protein_id", "contrast")]
  expect_equal(
    unlist(log_fc[2, c("pAUC_5", "pAUC_10")]),
    score_on_pairs(results$complete, selective_pairs, "log_fc")
  )
  expect_equal(scores$method[scores$score == "t_statistic"], c("selective", "complete"))
  # P2 in c2 has no p-value in complete, so complete is no superset for scaled.p_value; partial and equal never are
  expect_false("scaled.p_value" %in% scores$score)
  expect_false(any(c("partial", "equal") %in% scores$method))
})

test_that("same_pairs_scores returns the columns also when no front method has a superset", {
  results <- same_pairs_results()[c("partial", "equal")]
  table <- same_pairs_score_table()
  scores <- same_pairs_scores(results, table[table$model_name %in% c("partial", "equal"), ])
  expect_equal(nrow(scores), 0)
  expect_equal(colnames(scores), c("score", "front_method", "proteins", "method", "pAUC_5", "pAUC_10"))
})
