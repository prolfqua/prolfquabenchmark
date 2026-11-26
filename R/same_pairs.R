# Score methods on the same protein-contrast pairs -----

# Columns that must be defined for a pair to count as tested with a ranking score (ContrastsComparison$pairs() roles):
# the rows each score ranks, the rows n_total counts
.pair_columns <- list(log_fc = "effect", t_statistic = "statistic", scaled.p_value = c("effect", "p.value"))

#' Points no other point beats in both coordinates
#'
#' A point is Pareto optimal when no other point is at least as high in both `x` and `y` and higher in one of them.
#' The benchmark uses it for the methods with the best trade-off between ranking accuracy (pAUC) and the number of
#' contrasts computed.
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
  # the ROC of this score only: another score may be missing for all of these pairs
  roc <- benchmark_from_result(result, scores = score)$get_confusion_benchmark()
  roc <- roc[roc$contrast == "all" & roc$what == score, ]
  c(pAUC_5 = ms_bench_auc(roc$FPR, roc$TPR, 0.05), pAUC_10 = ms_bench_auc(roc$FPR, roc$TPR, 0.1))
}

#' Scores of methods on the protein and contrast pairs of a method on the Pareto front
#'
#' A method on the Pareto front of pAUC_5 against contrasts computed may be there because it tests fewer proteins, an
#' easier set. For each ranking score and each front method, the other methods that test every one of its pairs are
#' scored on exactly these pairs with [score_on_pairs()]. The pairs a method tests come from
#' [prolfqua::ContrastsComparison]: with `log_fc` the pairs with a fold change, with `t_statistic` those with a t
#' statistic, with `scaled.p_value` those with a fold change and a p-value.
#'
#' @param results named list (names: methods) of the results of one dataset, each as for [score_on_pairs()]
#' @param scores the score table rows of these methods (`model_name`, `contrast`, `score`, `pAUC_5`, `n_total`), as in
#'   `benchmark_scores.tsv`; the front uses pAUC_5 and `n_total` of contrast `all`, all contrasts pooled
#' @return data.frame with one row per ranking score, front method with a superset and method (the front method on
#'   its own pairs first, then its supersets): `score`, `front_method`, `contrasts` (the front method's contrasts
#'   computed, protein and contrast pairs scored), `method`, `pAUC_5`, `pAUC_10`
#' @export
#' @family same pairs
same_pairs_scores <- function(results, scores) {
  comparison <- prolfqua::ContrastsComparison$new(lapply(results, .as_contrasts_table))
  summary <- scores |>
    dplyr::summarize(
      pAUC_5 = .data$pAUC_5[.data$contrast == "all"],
      contrasts = .data$n_total[.data$contrast == "all"],
      .by = c("model_name", "score")
    )
  rows <- list()
  for (score in names(.pair_columns)) {
    score_summary <- summary[summary$score == score, ]
    front <- pareto_optimal(score_summary$pAUC_5, score_summary$contrasts)
    for (i in which(front)) {
      front_method <- score_summary$model_name[i]
      rows[[length(rows) + 1]] <- .front_scores(comparison, results, score, front_method, score_summary$contrasts[i])
    }
  }
  empty <- data.frame(
    score = character(),
    front_method = character(),
    contrasts = numeric(),
    method = character(),
    pAUC_5 = numeric(),
    pAUC_10 = numeric()
  )
  dplyr::bind_rows(empty, rows)
}

#' Scores of one front method and of the other methods that test all its pairs, on its pairs; NULL without such a method
#' @param comparison prolfqua::ContrastsComparison of the methods
#' @param results named list of the methods' results
#' @param score ranking score
#' @param front_method method on the Pareto front
#' @param contrasts contrasts the front method computes (protein and contrast pairs scored)
#' @noRd
.front_scores <- function(comparison, results, score, front_method, contrasts) {
  columns <- .pair_columns[[score]]
  pair_keys <- function(method) do.call(paste, c(comparison$pairs(method, columns), sep = "\r"))
  front_keys <- pair_keys(front_method)
  # every pair of the front method; the same pairs count, the front method itself not
  superset <- vapply(
    names(results),
    function(method) method != front_method && all(front_keys %in% pair_keys(method)),
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
    contrasts = contrasts,
    method = methods,
    pAUC_5 = pauc["pAUC_5", ],
    pAUC_10 = pauc["pAUC_10", ],
    row.names = NULL
  )
}

#' Scores of every method on the sets of contrasts of the Pareto front and on the contrasts each next set adds
#'
#' The methods on the Pareto front of pAUC_5 against contrasts computed define sets of protein and contrast pairs, their
#' pairs, sorted here by size; front methods with the same pairs form one set, named by the first in alphabetical order.
#' Every method is scored with [score_on_pairs()] on each set if it tests all of its pairs,
#' and on the pairs the next larger set adds if it tests all of those: does a method that computes more contrasts rank
#' the added ones well? The pairs a method tests are those of [same_pairs_scores()].
#'
#' @inheritParams same_pairs_scores
#' @return data.frame with one row per ranking score, front set and method: `score`, `front_method` (the method whose
#'   pairs form the set), `contrasts` (pairs in the set), `method`, `pAUC_5` (on the set; NA if the method does not
#'   test all its pairs), `next_contrasts` (pairs in the next larger set), `added` (pairs of the next set not in this
#'   one), `added_positives` (changed pairs among them, the positives of the ground truth), `pAUC_5_added` (on the added
#'   pairs; NA if the method does not test all of them); the last four are NA for the largest set
#' @export
#' @family same pairs
front_set_scores <- function(results, scores) {
  comparison <- prolfqua::ContrastsComparison$new(lapply(results, .as_contrasts_table))
  summary <- scores |>
    dplyr::summarize(
      pAUC_5 = .data$pAUC_5[.data$contrast == "all"],
      contrasts = .data$n_total[.data$contrast == "all"],
      .by = c("model_name", "score")
    )
  pair_keys <- function(pairs) do.call(paste, c(pairs, sep = "\r"))
  score_if_tested <- function(method_keys, result, pairs, score) {
    if (!all(pair_keys(pairs) %in% method_keys)) {
      return(NA_real_)
    }
    score_on_pairs(result, pairs, score)[["pAUC_5"]]
  }
  rows <- list()
  for (score in names(.pair_columns)) {
    columns <- .pair_columns[[score]]
    score_summary <- summary[summary$score == score, ]
    front <- score_summary$model_name[pareto_optimal(score_summary$pAUC_5, score_summary$contrasts)]
    sets <- lapply(stats::setNames(front, front), function(method) comparison$pairs(method, columns))
    sets <- sets[order(vapply(sets, nrow, integer(1)), names(sets))]
    # front methods with the same pairs form one set
    sets <- sets[!duplicated(lapply(sets, function(pairs) sort(pair_keys(pairs))))]
    method_keys <- lapply(stats::setNames(names(results), names(results)), function(m) {
      pair_keys(comparison$pairs(m, columns))
    })
    for (i in seq_along(sets)) {
      pairs <- sets[[i]]
      row <- data.frame(
        score = score,
        front_method = names(sets)[i],
        contrasts = nrow(pairs),
        method = names(results),
        pAUC_5 = vapply(
          names(results),
          function(m) score_if_tested(method_keys[[m]], results[[m]], pairs, score),
          numeric(1)
        ),
        next_contrasts = NA_integer_,
        added = NA_integer_,
        added_positives = NA_integer_,
        pAUC_5_added = NA_real_,
        row.names = NULL
      )
      if (i < length(sets)) {
        next_method <- names(sets)[i + 1]
        added <- sets[[i + 1]][!pair_keys(sets[[i + 1]]) %in% pair_keys(pairs), ]
        truth <- annotate_ground_truth(
          dplyr::semi_join(results[[next_method]]$data, added, by = c("protein_id", "contrast")),
          results[[next_method]]$metadata$ground_truth
        )
        row$next_contrasts <- nrow(sets[[i + 1]])
        row$added <- nrow(added)
        row$added_positives <- sum(truth$TP)
        row$pAUC_5_added <- vapply(
          names(results),
          function(m) score_if_tested(method_keys[[m]], results[[m]], added, score),
          numeric(1)
        )
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  dplyr::bind_rows(rows)
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
