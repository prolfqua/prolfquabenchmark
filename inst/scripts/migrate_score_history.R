#!/usr/bin/env Rscript
# One-off (2026-11-23): relabels the score history recorded by the benchmark vignettes to the dataset and method ids
# of the benchmark pipeline (config/config.yaml) and adds the column `level`. Values are unchanged. Drops the
# null-split rows (study removed) and the `limma` rows of the ROPECA vignette (the same model as `prolfqua_limma`).
#   Rscript inst/scripts/migrate_score_history.R inst/extdata/score_history.tsv

path <- commandArgs(trailingOnly = TRUE)[[1]]
history <- utils::read.delim(
  path,
  stringsAsFactors = FALSE,
  colClasses = c(prolfqua_version = "character", date = "character", tool_version = "character")
)

history <- history[history$dataset != "IonStar_MaxQuant_null_split", ]
history <- history[!(history$dataset == "IonStar_MaxQuant" & history$method == "limma"), ]

# proDA ran on two MaxQuant inputs: LFQ intensities (proteinGroups.txt) and median-polished peptides
lfq <- history$dataset == "IonStar_MQ_LFQ"
history$method[lfq] <- "proDA_LFQ"
dataset_ids <- c(
  IonStar_MQ = "IonStar_MaxQuant",
  IonStar_MQ_LFQ = "IonStar_MaxQuant",
  IonStar_MQ_medpolish = "IonStar_MaxQuant",
  IonStar_FragPipe = "IonStar_FragPipe_MSstats",
  IonStar_FragPipe_MSstatsCSV = "IonStar_FragPipe_MSstats",
  IonStar_FragPipe_combined_protein = "IonStar_FragPipe_protein"
)
renamed <- history$dataset %in% names(dataset_ids)
history$dataset[renamed] <- dataset_ids[history$dataset[renamed]]

peptide <- c(
  "prolfqua_mixed",
  "prolfqua_mixedModerated",
  "prolfqua_ropeca",
  "lmer_nested",
  "ropeca_nested",
  "ropecaV2_nested"
)
history$level <- ifelse(
  history$metric == "elapsed_s",
  NA_character_,
  ifelse(
    history$method %in% peptide,
    "peptide",
    ifelse(grepl("^hurdle_|^msqrob2", history$method), "protein + peptide", "protein")
  )
)

keys <- c("prolfqua_version", "dataset", "method", "contrast", "score", "metric")
stopifnot(!anyDuplicated(history[, keys]))
columns <- c(
  "prolfqua_version",
  "date",
  "dataset",
  "method",
  "level",
  "tool_version",
  "contrast",
  "score",
  "metric",
  "value"
)
utils::write.table(history[, columns], path, sep = "\t", row.names = FALSE, quote = FALSE)
print(table(history$prolfqua_version, history$dataset))
