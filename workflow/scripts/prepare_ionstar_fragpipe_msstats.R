# IonStar, FragPipe MSstats.csv: precursor intensities summed to peptides, all proteins (no minimum number of
# peptides), log2 and robust scaling on the human proteins, median polish to proteins.
#   Rscript workflow/scripts/prepare_ionstar_fragpipe_msstats.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))
msstats <- tibble::tibble(utils::read.csv(
  unz(quantdata_file("MSFragger_IonStar2018_PXD003881.zip"), "IonstarWithMSFragger/MSstats.csv"),
  header = TRUE,
  sep = ",",
  stringsAsFactors = FALSE
))
msstats$Run <- tolower(msstats$Run)
peptides <- dplyr::inner_join(annotation, msstats, by = c("raw.file" = "Run")) |>
  dplyr::group_by(.data$raw.file, .data$run_ID, .data$sample, .data$ProteinName, .data$PeptideSequence) |>
  dplyr::summarize(nprecursor = dplyr::n(), pepIntensity = sum(.data$Intensity, na.rm = TRUE), .groups = "drop")

config <- prolfqua::AnalysisConfiguration$new()
config$file_name <- "raw.file"
config$hierarchy[["protein_Id"]] <- "ProteinName"
config$hierarchy[["peptide_Id"]] <- "PeptideSequence"
config$hierarchy_depth <- 1
config$set_response("pepIntensity")
config$factors[["dilution."]] <- "sample"
config$factors[["run"]] <- "run_ID"
config$factor_depth <- 1

lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(peptides, config), config)
lfq$remove_small_intensities()

normalized <- normalize_on_human(lfq)
aggregator <- normalized$get_Aggregator()
aggregator$aggregate()
save_prepared(args$out_dir, proteins = aggregator$lfq_agg, peptides = normalized)
