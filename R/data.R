#' Benchmark example result
#'
#' Contrast results of a prolfqua model on the IonStar E. coli spike-in
#' dilution series (PRIDE PXD003881), used by the examples and tests of the
#' benchmark functions. Moved here from prolfqua 1.9.0.
#'
#' @format A tibble with 16580 rows: `protein_Id`, `contrast`, the compared
#'   groups `c1`, `c2`, and `estimate`, `statistic`, `p.value`,
#'   `p.value.adjusted`, `moderated.p.value`, `moderated.p.value.adjusted`.
#' @source PRIDE PXD003881
#' @family data
#' @docType data
#' @keywords internal
"data_benchmarkExample"

#' Scores for the confusion matrix examples
#'
#' Example scores with known truth, used by the examples of the confusion
#' matrix and ROC functions. Moved here from prolfqua 1.9.0.
#'
#' @format A tibble with 4043 rows: `protein_Id`, `TP` (truth), `estimate`,
#'   `statistic`, `scaled.p.value`.
#' @family data
#' @docType data
#' @keywords internal
"data_test_confusion_matrix_scores"
