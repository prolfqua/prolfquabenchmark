# Exports the peptide consistency weights of every dataset (results/peptide_consistency/<dataset>/scores.tsv and the
# priors.tsv next to it) for the vignette Peptide_consistency_weights, as peptide_consistency_scores.tsv and
# peptide_consistency_priors.tsv in <out_dir>.
#   Rscript workflow/scripts/export_peptide_consistency.R <out_dir> <scores files...>
source("workflow/scripts/common.R")
args <- script_args("out_dir", rest = "score_files")
write_combined <- function(paths, name) {
  data <- dplyr::bind_rows(lapply(paths, utils::read.delim, stringsAsFactors = FALSE))
  utils::write.table(data, file.path(args$out_dir, name), sep = "\t", row.names = FALSE, quote = FALSE)
}
write_combined(args$score_files, "peptide_consistency_scores.tsv")
write_combined(file.path(dirname(args$score_files), "priors.tsv"), "peptide_consistency_priors.tsv")
