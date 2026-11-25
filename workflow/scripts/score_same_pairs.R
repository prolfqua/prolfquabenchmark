# Scores the methods of one dataset on the protein and contrast pairs of each Pareto-front method they strictly contain
# (prolfquabenchmark::same_pairs_scores()), rounded to 4 digits.
#   Rscript workflow/scripts/score_same_pairs.R <dataset> <output.tsv> <methods...>
source("workflow/scripts/common.R")
args <- script_args(c("dataset", "output"), rest = "methods")

results <- lapply(stats::setNames(args$methods, args$methods), function(method) {
  dir <- file.path("results/contrasts", args$dataset, method)
  list(
    data = utils::read.delim(file.path(dir, "contrasts.tsv"), stringsAsFactors = FALSE),
    metadata = yaml::read_yaml(file.path(dir, "metadata.yaml"))
  )
})
scores <- dplyr::bind_rows(lapply(args$methods, function(method) {
  path <- file.path("results/scores", args$dataset, paste0(method, ".tsv"))
  utils::read.delim(path, stringsAsFactors = FALSE, colClasses = c(tool_version = "character"))
}))

same_pairs <- prolfquabenchmark::same_pairs_scores(results, scores)
same_pairs <- cbind(dataset = rep(args$dataset, nrow(same_pairs)), same_pairs)
same_pairs$pAUC_5 <- round(same_pairs$pAUC_5, 4)
same_pairs$pAUC_10 <- round(same_pairs$pAUC_10, 4)
utils::write.table(same_pairs, args$output, sep = "\t", row.names = FALSE, quote = FALSE)
