#' Create a Benchmark object from a contrast result directory
#'
#' Convenience factory that reads a contrast result directory (contrasts.tsv +
#' metadata.yaml), applies ground truth annotation, and creates a
#' \code{\link{Benchmark}} R6 object.
#'
#' @param path directory containing contrasts.tsv and metadata.yaml
#' @param ... additional arguments passed to \code{\link{make_benchmark}}
#' @return \code{\link{Benchmark}} R6 object
#' @export
#' @family benchmarking
#' @examples
#' \dontrun{
#' bench <- benchmark_from_file("/tmp/ionstar_lm")
#' bench$plot_ROC()
#' bench$pAUC()
#' }
benchmark_from_file <- function(path, ...) {
  res <- read_contrast_results(path)
  make_benchmark(
    prpr = res$data,
    model_name = res$metadata$method,
    model_description = res$metadata$method_description %||% res$metadata$method,
    fcestimate = "log_fc",
    toscale = c("p_value"),
    hierarchy = c("protein_id"),
    avgInt = "avg_intensity",
    benchmark = list(
      list(score = "log_fc", desc = TRUE),
      list(score = "t_statistic", desc = TRUE),
      list(score = "scaled.p_value", desc = TRUE)
    ),
    FDRvsFDP = list(list(score = "p_value_adjusted", desc = FALSE)),
    summarizeNA = "t_statistic",
    ...
  )
}


#' Write benchmark summary results to TSV
#'
#' Writes a benchmark_results.tsv file with AUC metrics per
#' method x contrast x score.
#'
#' @param summary_table data.frame as returned by
#'   \code{Benchmark$to_summary_table()} or \code{collect_benchmark_results()}
#' @param path file path for the output TSV
#' @return Invisibly returns the path written to.
#' @export
#' @family benchmarking
write_benchmark_results <- function(summary_table, path) {
  utils::write.table(summary_table, file = path, sep = "\t",
                     row.names = FALSE, quote = FALSE)
  message("Wrote benchmark results to: ", path)
  invisible(path)
}


#' Read benchmark summary results from TSV
#'
#' @param path file path to benchmark_results.tsv
#' @return data.frame with benchmark summary metrics
#' @export
#' @family benchmarking
read_benchmark_results <- function(path) {
  utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE)
}


#' Collect benchmark results from multiple files
#'
#' Reads multiple benchmark_results.tsv files and combines them into a
#' single data.frame for cross-method comparison.
#'
#' @param paths character vector of file paths to benchmark_results.tsv files
#' @return data.frame with all benchmark results combined
#' @export
#' @family benchmarking
#' @examples
#' \dontrun{
#' files <- list.files("results", pattern = "benchmark_results.tsv",
#'                     recursive = TRUE, full.names = TRUE)
#' combined <- collect_benchmark_results(files)
#' }
collect_benchmark_results <- function(paths) {
  results <- lapply(paths, read_benchmark_results)
  dplyr::bind_rows(results)
}


#' Plot comparison of benchmark results across methods
#'
#' Creates grouped bar charts of AUC metrics by method and contrast.
#'
#' @param combined data.frame as returned by \code{collect_benchmark_results()}
#' @return ggplot object
#' @export
#' @family benchmarking
plot_benchmark_comparison <- function(combined) {
  long <- tidyr::pivot_longer(
    combined,
    cols = c("AUC", "pAUC_10", "pAUC_20"),
    names_to = "metric",
    values_to = "value"
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(x = .data$contrast, y = .data$value, fill = .data$model_name)
  ) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::facet_wrap(~ metric + score, scales = "free_y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::labs(x = "Contrast", y = "Score", fill = "Method")
}
