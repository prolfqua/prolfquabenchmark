# Score methods on the same protein-contrast pairs -----

# Columns that must be defined for a pair to count as tested with a ranking score (ContrastsComparison$pairs() roles):
# the rows each score ranks, the rows n_total counts
.pair_columns <- list(log_fc = "effect", t_statistic = "statistic", scaled.p_value = c("effect", "p.value"))

#' Points no other point beats in both coordinates
#'
#' A point is Pareto optimal when no other point is at least as high in both `x` and `y` and higher in one of them.
#' The benchmark uses it for the methods with the best trade-off between ranking accuracy (pAUC) and the number of
#' proteins tested.
#'
#' @param x numeric, higher is better
#' @param y numeric, the same length as `x`, higher is better
#' @return logical, one per point
#' @export
#' @examples
#' pareto_optimal(c(1, 2, 3), c(3, 2, 1))
#' pareto_optimal(c(1, 2), c(1, 2))
pareto_optimal <- function(x, y) {
  vapply(seq_along(x), function(i) !any(x >= x[i] & y >= y[i] & (x > x[i] | y > y[i])), logical(1))
}

#' Partial AUC of benchmark results restricted to protein and contrast pairs
#'
#' Keeps the rows of `result` in `pairs` and scores them as the benchmark does: [benchmark_from_result()] annotates
#' the ground truth and ranks the proteins by `score`; the partial areas under the ROC curve are those of all contrasts
#' pooled.
#'
#' @param result list with `data`, a contrast table in the format of `contrasts.tsv` (`protein_id`, `contrast`,
#'   `log_fc`, `t_statistic`, `p_value`, `p_value_adjusted`, `avg_intensity`), and `metadata` with the ground truth,
#'   as written by [write_contrast_results()]
#' @param pairs data.frame with the columns `protein_id` and `contrast`
#' @param score ranking score: `"log_fc"`, `"t_statistic"` or `"scaled.p_value"`
#' @return named numeric: `pAUC_5` and `pAUC_10`, the partial areas up to a false positive rate of 0.05 and 0.1
#' @export
#' @family same pairs
score_on_pairs <- function(result, pairs, score) {
  score <- match.arg(score, names(.pair_columns))
  result$data <- dplyr::semi_join(result$data, pairs, by = c("protein_id", "contrast"))
  roc <- benchmark_from_result(result)$get_confusion_benchmark()
  roc <- roc[roc$contrast == "all" & roc$what == score, ]
  c(pAUC_5 = ms_bench_auc(roc$FPR, roc$TPR, 0.05), pAUC_10 = ms_bench_auc(roc$FPR, roc$TPR, 0.1))
}

#' Scores of methods on the protein and contrast pairs of a method on the Pareto front
#'
#' A method on the Pareto front of pAUC_5 against proteins tested may be there because it tests fewer proteins, an
#' easier set. For each ranking score and each front method, the methods that test a strict superset of its pairs are
#' scored on exactly these pairs with [score_on_pairs()]. The pairs a method tests come from
#' [prolfqua::ContrastsComparison]: with `log_fc` the pairs with a fold change, with `t_statistic` those with a t
#' statistic, with `scaled.p_value` those with a fold change and a p-value.
#'
#' @param results named list (names: methods) of the results of one dataset, each as for [score_on_pairs()]
#' @param scores the score table rows of these methods (`model_name`, `contrast`, `score`, `pAUC_5`, `n_total`), as in
#'   `benchmark_scores.tsv`; the front uses pAUC_5 of contrast `all` and `n_total` averaged over the other contrasts
#' @return data.frame with one row per ranking score, front method with a superset and method (the front method on
#'   its own pairs first, then its supersets): `score`, `front_method`, `proteins` (the front method's proteins tested
#'   per contrast), `method`, `pAUC_5`, `pAUC_10`
#' @export
#' @family same pairs
same_pairs_scores <- function(results, scores) {
  comparison <- prolfqua::ContrastsComparison$new(lapply(results, .as_contrasts_table))
  summary <- scores |>
    dplyr::summarize(
      pAUC_5 = .data$pAUC_5[.data$contrast == "all"],
      proteins = mean(.data$n_total[.data$contrast != "all"]),
      .by = c("model_name", "score")
    )
  rows <- list()
  for (score in names(.pair_columns)) {
    score_summary <- summary[summary$score == score, ]
    front <- pareto_optimal(score_summary$pAUC_5, score_summary$proteins)
    for (i in which(front)) {
      front_method <- score_summary$model_name[i]
      rows[[length(rows) + 1]] <- .front_scores(comparison, results, score, front_method, score_summary$proteins[i])
    }
  }
  empty <- data.frame(
    score = character(),
    front_method = character(),
    proteins = numeric(),
    method = character(),
    pAUC_5 = numeric(),
    pAUC_10 = numeric()
  )
  dplyr::bind_rows(empty, rows)
}

#' Scores of one front method and of its strict supersets on its pairs; NULL without a superset
#' @param comparison prolfqua::ContrastsComparison of the methods
#' @param results named list of the methods' results
#' @param score ranking score
#' @param front_method method on the Pareto front
#' @param proteins proteins the front method tests per contrast
#' @noRd
.front_scores <- function(comparison, results, score, front_method, proteins) {
  columns <- .pair_columns[[score]]
  pair_keys <- function(method) do.call(paste, c(comparison$pairs(method, columns), sep = "\r"))
  front_keys <- pair_keys(front_method)
  superset <- vapply(
    names(results),
    function(method) {
      keys <- pair_keys(method)
      length(keys) > length(front_keys) && all(front_keys %in% keys)
    },
    logical(1)
  )
  if (!any(superset)) {
    return(NULL)
  }
  methods <- c(front_method, names(results)[superset])
  front_pairs <- comparison$pairs(front_method, columns)
  pauc <- vapply(methods, function(method) score_on_pairs(results[[method]], front_pairs, score), numeric(2))
  data.frame(
    score = score,
    front_method = front_method,
    proteins = proteins,
    method = methods,
    pAUC_5 = pauc["pAUC_5", ],
    pAUC_10 = pauc["pAUC_10", ],
    row.names = NULL
  )
}

#' A benchmark result as a prolfqua ContrastsTable, its columns named by a ContrastConfiguration
#' @param result list with `data` and `metadata`
#' @noRd
.as_contrasts_table <- function(result) {
  config <- prolfqua::ContrastConfiguration$new(
    subject_id = "protein_id",
    effect_col = "log_fc",
    score_col = "t_statistic",
    pvalue_col = "p_value",
    fdr_col = "p_value_adjusted",
    avg_abundance_col = "avg_intensity"
  )
  prolfqua::ContrastsTable$new(
    result$data,
    subject_id = "protein_id",
    model_name = result$metadata$method,
    config = config
  )
}
