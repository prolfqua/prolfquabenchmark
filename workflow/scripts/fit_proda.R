# Fits proDA to the prepared proteins of one dataset (input: proteins), or to the MaxQuant LFQ protein intensities
# of proteinGroups.txt (input: lfq_proteins), with the same normalization as the prepared datasets.
#   Rscript workflow/scripts/fit_proda.R <config> <dataset> <method> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "method", "prepared_dir", "out_dir"))
spec <- read_spec(args$config, args$dataset, args$method)

# MaxQuant proteinGroups.txt LFQ intensities: contaminants and decoys removed, proteins with at least 2 peptides,
# log2 and robust scaling on the human proteins
mq_lfq_experiment <- function() {
  proteins <- prolfquapp::tidyMQ_ProteinGroups(quantdata_file("MAXQuant_IonStar2018_PXD003881.zip"))
  proteins <- dplyr::filter(proteins, !grepl("^REV__|^CON__", .data$proteinID))
  annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))
  proteins <- dplyr::filter(dplyr::inner_join(proteins, annotation, by = "raw.file"), .data$nr.peptides > 1)
  config <- prolfqua::AnalysisConfiguration$new()
  config$factors[["dilution."]] <- "sample"
  config$file_name <- "raw.file"
  config$set_response("mq.protein.lfq.intensity")
  config$hierarchy[["proteinID"]] <- "proteinID"
  lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(proteins, config), config)
  lfq$remove_small_intensities()
  human <- lfq$get_subset(dplyr::filter(lfq$data_long(), grepl("HUMAN", .data$proteinID)))
  transformer <- lfq$get_Transformer()
  transformer$log2()
  transformer$robscale_subset(human$get_Transformer()$log2()$lfq)
  wide <- transformer$lfq$data_wide(as.matrix = TRUE)
  sample_annotation <- data.frame(wide$annotation)
  rownames(sample_annotation) <- wide$annotation$sampleName
  SummarizedExperiment::SummarizedExperiment(S4Vectors::SimpleList(LFQ = wide$data), colData = sample_annotation)
}

experiment <- switch(
  spec$method$input,
  proteins = prolfqua::lfqdata_to_summarized_experiment(readRDS(file.path(args$prepared_dir, "proteins.rds"))),
  lfq_proteins = mq_lfq_experiment(),
  stop("proDA has no input ", spec$method$input)
)
# proDA's fit draws random numbers (two unseeded fits of the same data differ by up to 1.3 in log fold change)
set.seed(1)
fit <- proDA::proDA(
  experiment,
  design = stats::as.formula(paste(spec$dataset$formula, "- 1")),
  data_is_log_transformed = TRUE
)
result <- dplyr::bind_rows(lapply(names(spec$contrasts), function(name) {
  data.frame(contrast = name, proDA::test_diff(fit, contrast = spec$contrasts[[name]]))
}))
result <- dplyr::rename(
  result,
  protein_id = "name",
  log_fc = "diff",
  p_value = "pval",
  p_value_adjusted = "adj_pval",
  avg_intensity = "avg_abundance"
)
write_result(
  result,
  args$out_dir,
  spec,
  args$dataset,
  args$method,
  software_version = paste0("proDA ", utils::packageVersion("proDA"))
)
