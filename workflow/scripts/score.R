# Scores one result directory: the summary table of Benchmark$to_summary_table() (AUC, pAUC, AP per contrast and
# ranking score) plus the partial AUC up to a false positive rate of 0.05 (pAUC_5), and the curves of the pooled
# contrast on fixed grids (TPR at FPR 0 to 0.2 by 0.002 per ranking score; false discovery proportion, curve FDP, and
# number of true discoveries, curve TP, at nominal FDR 0 to 0.25 by 0.0025), rounded to 4 digits.
#   Rscript workflow/scripts/score.R <contrast_dir> <scores.tsv> <curves.tsv>
source("workflow/scripts/common.R")
args <- script_args(c("contrast_dir", "scores", "curves"))

metadata <- yaml::read_yaml(file.path(args$contrast_dir, "metadata.yaml"))
bench <- prolfquabenchmark::benchmark_from_file(args$contrast_dir)

roc <- bench$get_confusion_benchmark()
pauc_5 <- roc |>
  dplyr::summarize(pAUC_5 = prolfquabenchmark::ms_bench_auc(.data$FPR, .data$TPR, 0.05), .by = c("contrast", "what")) |>
  dplyr::rename(score = "what")
scores <- bench$to_summary_table(dataset = metadata$dataset) |>
  dplyr::left_join(pauc_5, by = c("contrast", "score")) |>
  dplyr::relocate("pAUC_5", .before = "pAUC_10")
scores$level <- metadata$level
scores$tool_version <- metadata$software_version
utils::write.table(scores, args$scores, sep = "\t", row.names = FALSE, quote = FALSE)

# y at the last point of the curve with x <= each grid value (NA before the first point)
on_grid <- function(x, y, grid) {
  order_x <- order(x)
  x <- x[order_x]
  y <- y[order_x]
  round(vapply(grid, function(g) if (any(x <= g)) y[max(which(x <= g))] else NA_real_, numeric(1)), 4)
}

roc <- dplyr::filter(roc, .data$contrast == "all")
fpr_grid <- seq(0, 0.2, by = 0.002)
roc_curves <- lapply(split(roc, roc$what), function(curve) {
  data.frame(curve = "ROC", score = curve$what[1], x = fpr_grid, y = on_grid(curve$FPR, curve$TPR, fpr_grid))
})
fdp <- dplyr::filter(bench$get_confusion_FDRvsFDP(), .data$contrast == "all")
fdr_grid <- seq(0, 0.25, by = 0.0025)
fdp_curves <- lapply(split(fdp, fdp$what), function(curve) {
  data.frame(curve = "FDP", score = curve$what[1], x = fdr_grid, y = on_grid(curve$scorecol, curve$FDP, fdr_grid))
})
tp_curves <- lapply(split(fdp, fdp$what), function(curve) {
  data.frame(curve = "TP", score = curve$what[1], x = fdr_grid, y = on_grid(curve$scorecol, curve$TP_hits, fdr_grid))
})
curves <- cbind(
  dataset = metadata$dataset,
  method = metadata$method,
  do.call(rbind, c(roc_curves, fdp_curves, tp_curves))
)
utils::write.table(curves, args$curves, sep = "\t", row.names = FALSE, quote = FALSE)
