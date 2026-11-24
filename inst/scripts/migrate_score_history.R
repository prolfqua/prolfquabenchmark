#!/usr/bin/env Rscript
# One-off (2026-11-23): relabels the score history recorded by the benchmark vignettes to the dataset and method ids
# of the benchmark pipeline (config/config.yaml) and adds the column `level`. Values are unchanged. Drops the
# null-split rows (study removed) and the `limma` rows of the ROPECA vignette (the same model as `prolfqua_limma`).
# Second step (2026-11-24): drops prolfqua 1.4.0 and the hurdle variants no longer benchmarked, and names
# `hurdle_binomial_firth` after prolfqua's `hurdle_nested` facade, which fits the same model. Rerunning is safe.
#   Rscript inst/scripts/migrate_score_history.R inst/extdata/score_history.tsv

path <- commandArgs(trailingOnly = TRUE)[[1]]
history <- utils::read.delim(
  path,
  stringsAsFactors = FALSE,
  colClasses = c(prolfqua_version = "character", date = "character", tool_version = "character", level = "character")
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
if (!"level" %in% names(history)) {
  history$level <- ifelse(
    history$metric == "elapsed_s",
    NA_character_,
    ifelse(
      history$method %in% peptide,
      "peptide",
      ifelse(grepl("^hurdle_|^msqrob2", history$method), "protein + peptide", "protein")
    )
  )
}

# 2026-11-24: prolfqua 1.4.0 shares too few models and vignettes with today's; the hurdle variants were development
# experiments; hurdle_nested (lm_impute on median-polished proteins, binomial_firth_nested) scores as
# hurdle_binomial_firth did
history <- history[history$prolfqua_version != "1.4.0", ]
dropped <- c(
  "hurdle_binomial",
  "hurdle_binomial_firth_unbound",
  "hurdle_binomial_firth_df4",
  "hurdle_firth",
  "hurdle_firth_df4"
)
history <- history[!history$method %in% dropped, ]
history$method[history$method == "hurdle_binomial_firth"] <- "hurdle_nested"

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
