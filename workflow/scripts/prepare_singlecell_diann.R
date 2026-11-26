# ProteoBench DIA single-cell module (PXD049412, Orbitrap Astral; human/yeast 240/10 pg in A, 200/50 pg in B), DIA-NN
# 2.5.0 submission 239b93c0: see prepare_proteobench() and proteobench_diann_peptides().
#   Rscript workflow/scripts/prepare_singlecell_diann.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

file <- proteobench_submission("239b93c09a51038aab992c15c0702ff3658fe1c1", "a536ceaa97f65183f71518f9bf07c570")
species <- proteobench_species("ProteoBenchFASTA_MixedSpecies_HY.zip", "ca94c5b62a21b7a3235d7555dd94c559")
prepare_proteobench(
  proteobench_diann_peptides(file),
  species,
  conditions = c(A = "240pg_10pg", B = "200pg_50pg"),
  out_dir = args$out_dir
)
