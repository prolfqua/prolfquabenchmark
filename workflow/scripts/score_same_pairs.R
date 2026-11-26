# Scores the methods of one dataset on the protein and contrast pairs of each Pareto-front method whose pairs they all
# test (prolfquabenchmark::same_pairs_scores()), and every method on each front set and on the pairs the next larger
# set adds (prolfquabenchmark::front_set_scores()), rounded to 4 digits.
#   Rscript workflow/scripts/score_same_pairs.R <dataset> <same_pairs.tsv> <front_sets.tsv> <methods...>
source("workflow/scripts/common.R")
args <- script_args(c("dataset", "output", "front_sets"), rest = "methods")

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
utils::write.table(same_pairs, args$output, sep = "\t", row.names = FALSE, quote = FALSE)

front_sets <- prolfquabenchmark::front_set_scores(results, scores)
front_sets <- cbind(dataset = rep(args$dataset, nrow(front_sets)), front_sets)
front_sets$pAUC_5 <- round(front_sets$pAUC_5, 4)
front_sets$pAUC_5_added <- round(front_sets$pAUC_5_added, 4)
utils::write.table(front_sets, args$front_sets, sep = "\t", row.names = FALSE, quote = FALSE)
