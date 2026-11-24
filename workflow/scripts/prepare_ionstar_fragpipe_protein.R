# IonStar, FragPipe combined_protein.tsv (protein level): proteins with at least 2 peptides, razor intensities, log2
# and robust scaling on the human proteins; the number of unique peptides is the child count (fit weights, DEqMS).
#   Rscript workflow/scripts/prepare_ionstar_fragpipe_protein.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

annotation <- readxl::read_xlsx(quantdata_file("annotation_Ionstar2018_PXD003881.xlsx"))
protein <- tibble::as_tibble(utils::read.csv(
  unz(quantdata_file("MSFragger_IonStar2018_PXD003881.zip"), "IonstarWithMSFragger/combined_protein.tsv"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
))
protein <- prolfquapp::tidy_FragPipe_combined_protein_deprec(protein)
protein <- dplyr::filter(protein, .data$unique.stripped.peptides > 1)
protein <- dplyr::inner_join(annotation, protein)

config <- prolfqua::AnalysisConfiguration$new()
config$file_name <- "raw.file"
config$hierarchy[["protein_Id"]] <- "protein"
config$hierarchy_depth <- 1
config$set_response("razor.intensity")
config$nr_children <- "unique.stripped.peptides"
config$factors[["dilution."]] <- "sample"
config$factors[["run"]] <- "run_ID"
config$factor_depth <- 1

lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(protein, config), config)
lfq$remove_small_intensities()
save_prepared(args$out_dir, proteins = normalize_on_human(lfq))
