#!/usr/bin/env Rscript
# One-off (2026-11-23): relabels the score history recorded by the benchmark vignettes to the dataset and method ids
# of the benchmark pipeline (config/config.yaml) and adds the column `level`. Values are unchanged. Drops the
# null-split rows (study removed) and the `limma` rows of the ROPECA vignette (the same model as `prolfqua_limma`).
# Second step (2026-11-24): drops prolfqua 1.4.0 and the hurdle variants no longer benchmarked, and names
# `hurdle_binomial_firth` after prolfqua's `hurdle_nested` facade, which fits the same model. Third step (2026-11-24):
# drops `prolfqua_lm` (lm_missing, now in prolfquatakeoutmissing, which the benchmark no longer uses), names the
# prolfqua models `prolfqua_<facade>`, and keeps proDA and msqrob2_protein only on the search software's protein
# table. Fourth step (2026-11-24): one input file per dataset; proDA and msqrob2_protein on IonStar_MaxQuant move to
# IonStar_MaxQuant_protein, MSstats on MaxQuant's evidence.txt (alone on its file) is dropped. Rerunning is safe.
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

# 2026-11-24, while the prolfqua ids still lack the prefix: lm_missing (prolfqua_lm) left prolfqua for
# prolfquatakeoutmissing; proDA on the MaxQuant LFQ intensities (proDA_LFQ) is proDA, proDA and msqrob2_protein on
# proteins aggregated by prolfqua are dropped; the prolfqua models are named prolfqua_<facade>
if ("lm_impute" %in% history$method) {
  history <- history[history$method != "prolfqua_lm", ]
  aggregated <- c("IonStar_MaxQuant", "IonStar_FragPipe_MSstats", "CPTAC_MaxQuant")
  history <- history[!(history$method %in% c("proDA", "msqrob2_protein") & history$dataset %in% aggregated), ]
  history$method[history$method == "proDA_LFQ"] <- "proDA"
  facades <- c(
    "lm",
    "lm_impute",
    "deqms",
    "deqms_voom",
    "lmer_nested",
    "ropeca_nested",
    "ropecaV2_nested",
    "hurdle_nested"
  )
  renamed <- history$method %in% facades
  history$method[renamed] <- paste0("prolfqua_", history$method[renamed])
}

# 2026-11-24: one input file per dataset; MaxQuant's proteinGroups.txt is a dataset of its own, MSstats on its
# evidence.txt (the only method on that file) is dropped
maxquant <- history$dataset == "IonStar_MaxQuant"
history$dataset[maxquant & history$method %in% c("proDA", "msqrob2_protein")] <- "IonStar_MaxQuant_protein"
history <- history[!(maxquant & history$method == "MSstats") & history$dataset != "IonStar_MaxQuant_evidence", ]

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
