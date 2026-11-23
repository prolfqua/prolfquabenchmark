# CPTAC study 6 (UPS1 spiked into yeast), MaxQuant peptides.txt from Bioconductor MsDataHub: proteins with at least 2
# peptides, vsn, median polish to proteins.
#   Rscript workflow/scripts/prepare_cptac_maxquant.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

peptides <- prolfquapp::tidyMQ_Peptides(utils::read.csv(MsDataHub::cptac_a_b_peptides.txt(), sep = "\t"))
annotation <- data.frame(raw.file = unique(peptides$raw.file))
annotation$group <- gsub("_[0-9]$", "", gsub("^6", "", annotation$raw.file))
peptides <- dplyr::inner_join(annotation, peptides)

config <- create_mq_peptide_config()
config$factors["group."] <- "group"

lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(peptides, config), config)
lfq$filter_proteins_by_peptide_count()
lfq$remove_small_intensities()

normalized <- lfq$get_Transformer()$intensity_matrix(vsn::justvsn)$lfq
aggregator <- normalized$get_Aggregator()
aggregator$aggregate()
save_prepared(args$out_dir, proteins = aggregator$lfq_agg, peptides = normalized)
