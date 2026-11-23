# Fits MSstats (dataProcess, groupComparison) to its own input of one dataset: MaxQuant evidence.txt and
# proteinGroups.txt for IonStar_MaxQuant, FragPipe MSstats.csv for IonStar_FragPipe_MSstats.
#   Rscript workflow/scripts/fit_msstats.R <config> <dataset> <method> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "method", "prepared_dir", "out_dir"))
spec <- read_spec(args$config, args$dataset, args$method)
annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))

msstats_input <- switch(
  args$dataset,
  IonStar_MaxQuant = {
    zip <- quantdata_file("MAXQuant_IonStar2018_PXD003881.zip")
    protein_groups <- utils::read.table(unz(zip, "proteinGroups.txt"), sep = "\t", header = TRUE)
    evidence <- utils::read.table(unz(zip, "evidence.txt"), sep = "\t", header = TRUE)
    msstats_annotation <- dplyr::select(annotation, "raw.file", Condition = "sample", BioReplicate = "run_ID")
    msstats_annotation$IsotopeLabelType <- "L"
    msstats_annotation <- tidyr::unite(
      msstats_annotation,
      "Experiment",
      "Condition",
      "BioReplicate",
      sep = "_",
      remove = FALSE
    )
    raw_files <- data.frame(Raw.file = unique(evidence$Raw.file), raw.file = tolower(unique(evidence$Raw.file)))
    msstats_annotation <- dplyr::inner_join(raw_files, msstats_annotation, by = "raw.file")
    msstats_annotation$raw.file <- NULL
    raw <- MSstats::MaxQtoMSstatsFormat(
      evidence = evidence,
      annotation = msstats_annotation,
      proteinGroups = protein_groups,
      removeProtein_with1Peptide = TRUE,
      use_log_file = FALSE,
      verbose = FALSE
    )
    raw[!grepl("REV__|CON__", raw$ProteinName), ]
  },
  IonStar_FragPipe_MSstats = {
    fragpipe <- utils::read.csv(
      unz(quantdata_file("MSFragger_IonStar2018_PXD003881.zip"), "IonstarWithMSFragger/MSstats.csv"),
      header = TRUE,
      sep = ",",
      stringsAsFactors = FALSE
    )
    fragpipe$BioReplicate <- NULL
    fragpipe$Condition <- NULL
    fragpipe$Run <- tolower(fragpipe$Run)
    fragpipe_annotation <- dplyr::select(annotation, "raw.file", BioReplicate = "run_ID", Condition = "sample")
    dplyr::inner_join(fragpipe, fragpipe_annotation, by = c(Run = "raw.file"))
  },
  stop("MSstats has no input for dataset ", args$dataset)
)

# One row per contrast, one column per condition (the factor levels a, b, ...), from the contrast expressions
factor_col <- all.vars(stats::as.formula(spec$dataset$formula))
conditions <- sort(unique(as.character(msstats_input$Condition)))
weights <- t(vapply(
  spec$contrasts,
  function(expression) {
    vapply(
      conditions,
      function(level) {
        groups <- stats::setNames(as.list(as.numeric(conditions == level)), paste0(factor_col, conditions))
        eval(parse(text = expression), groups)
      },
      numeric(1)
    )
  },
  numeric(length(conditions))
))
rownames(weights) <- names(spec$contrasts)

invisible(utils::capture.output({
  processed <- MSstats::dataProcess(msstats_input, use_log_file = FALSE, verbose = FALSE)
  comparison <- MSstats::groupComparison(
    contrast.matrix = weights,
    data = processed,
    verbose = FALSE,
    use_log_file = FALSE
  )
}))
result <- dplyr::rename(
  comparison$ComparisonResult,
  protein_id = "Protein",
  contrast = "Label",
  log_fc = "log2FC",
  t_statistic = "Tvalue",
  p_value = "pvalue",
  p_value_adjusted = "adj.pvalue"
)
write_result(
  result,
  args$out_dir,
  spec,
  args$dataset,
  args$method,
  software_version = paste0("MSstats ", utils::packageVersion("MSstats"))
)
