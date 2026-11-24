# IonStar, MaxQuant peptides.txt: contaminants and decoys removed, proteins with at least 2 peptides, log2 and robust
# scaling on the human proteins, median polish to proteins. Protein table: the LFQ intensities of proteinGroups.txt,
# contaminants and decoys removed, proteins with at least 2 peptides, log2 and robust scaling on the human proteins.
#   Rscript workflow/scripts/prepare_ionstar_maxquant.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

peptides <- prolfquapp::tidyMQ_Peptides(quantdata_file("MAXQuant_IonStar2018_PXD003881.zip"))
annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))
data <- dplyr::inner_join(peptides, annotation, by = "raw.file")

config <- create_mq_peptide_config()
config$factors[["dilution."]] <- "sample"
config$factors[["run_Id"]] <- "run_ID"
config$factor_depth <- 1

lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(data, config), config)
lfq$set_data(dplyr::filter(lfq$data_long(), !grepl("^REV__|^CON__", .data$protein_Id)))
lfq$filter_proteins_by_peptide_count()
lfq$remove_small_intensities()

normalized <- normalize_on_human(lfq)
aggregator <- normalized$get_Aggregator()
aggregator$aggregate()

protein_groups <- prolfquapp::tidyMQ_ProteinGroups(quantdata_file("MAXQuant_IonStar2018_PXD003881.zip"))
protein_groups <- dplyr::filter(protein_groups, !grepl("^REV__|^CON__", .data$proteinID))
protein_groups <- dplyr::inner_join(protein_groups, annotation, by = "raw.file")
protein_groups <- dplyr::filter(protein_groups, .data$nr.peptides > 1)
table_config <- prolfqua::AnalysisConfiguration$new()
table_config$factors[["dilution."]] <- "sample"
table_config$file_name <- "raw.file"
table_config$set_response("mq.protein.lfq.intensity")
table_config$hierarchy[["protein_Id"]] <- "proteinID"
protein_table <- prolfqua::LFQData$new(prolfqua::setup_analysis(protein_groups, table_config), table_config)
protein_table$remove_small_intensities()

save_prepared(
  args$out_dir,
  proteins = aggregator$lfq_agg,
  peptides = normalized,
  protein_table = normalize_on_human(protein_table)
)
