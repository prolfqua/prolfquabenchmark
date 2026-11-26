# ProteoBench DIA single-cell module (PXD049412, Orbitrap Astral; human/yeast 240/10 pg in A, 200/50 pg in B),
# MaxQuant 2.7.5.0 (DIA) evidence.txt of submission 090d5f63, read with prolfquapp::tidyMQ_Evidence(): proteotypic
# evidence, intensities summed per run and peptide, protein the leading razor protein, confidence 1 - the smallest
# evidence PEP, floored at 0.5; then prepare_proteobench().
#   Rscript workflow/scripts/prepare_singlecell_maxquant.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

file <- proteobench_submission("090d5f63ab4b453cbbc318eb5f8ccfe557ab8621", "fae035ee9685f11e0372d86e4d18b548")
species <- proteobench_species("ProteoBenchFASTA_MixedSpecies_HY.zip", "ca94c5b62a21b7a3235d7555dd94c559")
evidence <- prolfquapp::tidyMQ_Evidence(file)
peptides <- evidence |>
  dplyr::filter(.data$proteotypic, !.data$reverse, !is.na(.data$evidence.intensity)) |>
  dplyr::summarise(
    intensity = sum(.data$evidence.intensity),
    identification_confidence = floor_confidence(1 - min(.data$pep)),
    .by = c("raw.file", "leading.razor.protein", "peptide.id")
  ) |>
  dplyr::transmute(
    .data$raw.file,
    protein = .data$leading.razor.protein,
    peptide = as.character(.data$peptide.id),
    .data$intensity,
    .data$identification_confidence
  )
prepare_proteobench(
  peptides,
  species,
  conditions = c(A = "240pg_10pg", B = "200pg_50pg"),
  out_dir = args$out_dir
)
