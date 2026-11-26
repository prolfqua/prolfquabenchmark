# ProteoBench DIA diaPASEF module (timsTOF; human/yeast/E. coli 65/30/5% in A, 65/15/20% in B), Spectronaut 21.0
# report of submission e0b85747, read with prolfquapp::read_BGS(): precursor quantities summed per run and peptide
# (PEP.GroupingKey), confidence 1 - the smallest precursor q-value floored at 0.5, as prolfquapp's BGS reader; then
# prepare_proteobench().
#   Rscript workflow/scripts/prepare_diapasef_spectronaut.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

file <- proteobench_submission("e0b85747554042cf392ea706aee617abcf70901f", "998810d0f555c6c02056e5be14741871")
species <- proteobench_species("ProteoBenchFASTA_MixedSpecies_HYE.zip", "8c680a816d98eb31b0c9684662c74ef9")
peptides <- prolfquapp::read_BGS(file) |>
  dplyr::filter(!is.na(.data$FG.Quantity)) |>
  dplyr::summarise(
    intensity = sum(.data$FG.Quantity),
    identification_confidence = floor_confidence(1 - min(.data$FG.Qvalue)),
    .by = c("R.FileName", "PG.ProteinGroups", "PEP.GroupingKey")
  ) |>
  dplyr::transmute(
    raw.file = .data$R.FileName,
    protein = .data$PG.ProteinGroups,
    peptide = .data$PEP.GroupingKey,
    .data$intensity,
    .data$identification_confidence
  )
prepare_proteobench(
  peptides,
  species,
  conditions = c(A = "Condition_A", B = "Condition_B"),
  out_dir = args$out_dir
)
