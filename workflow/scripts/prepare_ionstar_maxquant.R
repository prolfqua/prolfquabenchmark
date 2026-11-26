# IonStar, MaxQuant peptides.txt: contaminants and decoys removed, all proteins (no minimum number of peptides), log2
# and robust scaling on the human proteins, median polish to proteins.
#   Rscript workflow/scripts/prepare_ionstar_maxquant.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

peptides <- prolfquapp::tidyMQ_Peptides(quantdata_file("MAXQuant_IonStar2018_PXD003881.zip"))
annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))
data <- mq_identification_confidence(dplyr::inner_join(peptides, annotation, by = "raw.file"))

config <- create_mq_peptide_config()
config$factors[["dilution."]] <- "sample"
config$factors[["run_Id"]] <- "run_ID"
config$factor_depth <- 1

lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(data, config), config)
lfq$set_data(dplyr::filter(lfq$data_long(), !grepl("^REV__|^CON__", .data$protein_Id)))
lfq$remove_small_intensities()

normalized <- normalize_on_human(lfq)
aggregator <- normalized$get_Aggregator()
aggregator$aggregate()
save_prepared(args$out_dir, proteins = aggregator$lfq_agg, peptides = normalized)
