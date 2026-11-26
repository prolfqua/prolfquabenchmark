# ProteoBench DIA diaPASEF module (timsTOF; human/yeast/E. coli 65/30/5% in A, 65/15/20% in B), DIA-NN 2.5.0
# submission 53bb223d: see prepare_proteobench() and proteobench_diann_peptides().
#   Rscript workflow/scripts/prepare_diapasef_diann.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

file <- proteobench_submission("53bb223d26d0227ac35d403947eb7de30decf38d", "90d0639bf67ca5f9e30b68339c19680f")
species <- proteobench_species("ProteoBenchFASTA_MixedSpecies_HYE.zip", "8c680a816d98eb31b0c9684662c74ef9")
prepare_proteobench(
  proteobench_diann_peptides(file),
  species,
  conditions = c(A = "Condition_A", B = "Condition_B"),
  out_dir = args$out_dir
)
