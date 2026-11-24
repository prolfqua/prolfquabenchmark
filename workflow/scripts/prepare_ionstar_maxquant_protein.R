# IonStar, MaxQuant proteinGroups.txt (protein level): LFQ intensities, contaminants and decoys removed, proteins with
# at least 2 peptides, log2 and robust scaling on the human proteins; the number of peptides is the child count (fit
# weights, DEqMS). The proteins are the search software's protein table, so they are saved as both proteins and
# protein_table.
#   Rscript workflow/scripts/prepare_ionstar_maxquant_protein.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))
protein_groups <- prolfquapp::tidyMQ_ProteinGroups(quantdata_file("MAXQuant_IonStar2018_PXD003881.zip"))
protein_groups <- dplyr::filter(protein_groups, !grepl("^REV__|^CON__", .data$proteinID))
protein_groups <- dplyr::inner_join(protein_groups, annotation, by = "raw.file")
protein_groups <- dplyr::filter(protein_groups, .data$nr.peptides > 1)

config <- prolfqua::AnalysisConfiguration$new()
config$factors[["dilution."]] <- "sample"
config$file_name <- "raw.file"
config$set_response("mq.protein.lfq.intensity")
config$hierarchy[["protein_Id"]] <- "proteinID"
config$nr_children <- "nr.peptides"

lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(protein_groups, config), config)
lfq$remove_small_intensities()
proteins <- normalize_on_human(lfq)
save_prepared(args$out_dir, proteins = proteins, protein_table = proteins)
