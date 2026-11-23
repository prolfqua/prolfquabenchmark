#!/usr/bin/env Rscript
# Times the core prolfqua steps on the IonStar MaxQuant peptides with the installed prolfqua, for the score history
# (vignette `Score_history`). Median elapsed seconds over 3 repetitions of
#   transform     log2, then robust scaling on the human proteins
#   aggregate     median polish from peptides to proteins
#   fit_lm        build_model() with strategy_lm(), one lm per protein
#   contrasts_lm  Contrasts + ContrastsModerated on that model
#   facade_<m>    build_contrast_analysis() for lm, limma, lm_impute and rlm (prolfqua >= 1.6; 1.4 has no facades)
# Runs with prolfqua 1.4.0 and later; only the setup differs before 1.6 (AnalysisTableAnnotation configuration).
#   Rscript workflow/scripts/time_prolfqua_steps.R <input.rds> <out_dir>
# <input.rds> is written by timing_input.R; the rows go to <out_dir>/history_rows.tsv.

args <- commandArgs(trailingOnly = TRUE)
input <- readRDS(args[[1]])
out_dir <- args[[2]]

has_facades <- exists("build_contrast_analysis", envir = asNamespace("prolfqua"))
contrasts <- c(
  "dilution_(9/7.5)_1.2" = "dilution.e - dilution.d",
  "dilution_(7.5/6)_1.25" = "dilution.d - dilution.c",
  "dilution_(6/4.5)_1.3(3)" = "dilution.c - dilution.b",
  "dilution_(4.5/3)_1.5" = "dilution.b - dilution.a"
)

# prolfqua 1.4 keeps the data in the field `data`; later versions use data_long() and set_data()
data_of <- function(lfq) if (has_facades) lfq$data_long() else lfq$data
replace_data <- function(lfq, data) {
  if (has_facades) {
    lfq$set_data(data)
  } else {
    lfq$data <- data
  }
  invisible(lfq)
}

if (has_facades) {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "raw.file"
  config$hierarchy[["protein_Id"]] <- "leading.razor.protein"
  config$hierarchy[["peptide_Id"]] <- "sequence"
  config$hierarchy_depth <- 1
  config$ident_q_value <- "pep"
  config$set_response("peptide.intensity")
  config$min_peptides_protein <- 2
  if ("isotope_label" %in% names(config)) {
    config$isotope_label <- "isotope"
  }
  config$factors[["dilution."]] <- "sample"
  config$factors[["run_Id"]] <- "run_ID"
  config$factor_depth <- 1
} else {
  config <- prolfqua::create_config_MQ_peptide()
  config$table$factors[["dilution."]] <- "sample"
  config$table$factors[["run_Id"]] <- "run_ID"
  config$table$factorDepth <- 1
}
lfq <- suppressMessages(prolfqua::LFQData$new(prolfqua::setup_analysis(input, config), config))
replace_data(lfq, dplyr::filter(data_of(lfq), !grepl("^REV__|^CON__", .data$protein_Id)))
suppressMessages({
  lfq$filter_proteins_by_peptide_count()
  lfq$remove_small_intensities()
})

time_step <- function(f) {
  elapsed <- numeric(3)
  for (i in seq_along(elapsed)) {
    start <- proc.time()[["elapsed"]]
    value <- suppressMessages(suppressWarnings(f()))
    elapsed[i] <- proc.time()[["elapsed"]] - start
  }
  list(value = value, elapsed = stats::median(elapsed))
}

timings <- list()

step <- time_step(function() {
  human <- lfq$get_copy()
  replace_data(human, dplyr::filter(data_of(human), grepl("HUMAN", .data$protein_Id)))
  human <- human$get_Transformer()$log2()$lfq
  lfq$get_Transformer()$log2()$robscale_subset(human)$lfq
})
timings$transform <- step$elapsed
normalized <- step$value

step <- time_step(function() {
  aggregator <- normalized$get_Aggregator()
  if ("aggregate" %in% names(aggregator)) {
    aggregator$aggregate()
  } else {
    aggregator$medpolish()
  }
  aggregator$lfq_agg
})
timings$aggregate <- step$elapsed
protein <- step$value

step <- time_step(function() {
  prolfqua::build_model(
    protein,
    prolfqua::strategy_lm(paste0(protein$response(), " ~ dilution."), model_name = "Model")
  )
})
timings$fit_lm <- step$elapsed
model <- step$value

timings$contrasts_lm <- time_step(function() {
  prolfqua::ContrastsModerated$new(prolfqua::Contrasts$new(model, contrasts))$get_contrasts()
})$elapsed

if (has_facades) {
  for (method in c("lm", "limma", "lm_impute", "rlm")) {
    timings[[paste0("facade_", method)]] <- time_step(function() {
      prolfqua::build_contrast_analysis(protein, "~ dilution.", contrasts, method = method)$get_contrasts()
    })$elapsed
  }
}

rows <- data.frame(
  dataset = "IonStar_MaxQuant",
  method = names(timings),
  tool_version = paste0("R ", getRversion(), " ", R.version$platform),
  contrast = "all",
  score = "elapsed",
  metric = "elapsed_s",
  value = unlist(timings, use.names = FALSE)
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
utils::write.table(rows, file.path(out_dir, "history_rows.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
print(rows[, c("method", "value")], row.names = FALSE)
