# IonStar reanalysed with MetaMorpheus 1.0.7 and FlashLFQ 1.3 (PXD057758, Solivais et al., J Proteome Res 2025),
# QuantifiedPeaks.tsv from the APB example datasets on Zenodo (record 24418695). Decoys and random-RT controls removed;
# peak intensities summed to base sequences per run; identification confidence 1 for MS/MS peaks and 1 - PIP PEP for
# match-between-runs peaks, the largest of a peptide's peaks; protein groups of one organism (human or E. coli); then
# prepare_peptides() with the dilution of each run (a to e) from its file name.
#   Rscript workflow/scripts/prepare_ionstar_metamorpheus.R <out_dir>
source("workflow/scripts/common.R")
args <- script_args("out_dir")

zip <- downloaded_file(
  "https://zenodo.org/records/24418695/files/metamorpheus_ionstar.zip?download=1",
  "57a3ffb70970de7cfd93713c25f6300b",
  "results/downloads/zenodo"
)
peaks <- readr::read_tsv(
  extracted_file(zip, "QuantifiedPeaks\\.tsv$", "metamorpheus_ionstar"),
  col_select = c(
    "File Name",
    "Base Sequence",
    "Protein Group",
    "Organism",
    "Peak intensity",
    "Peak Detection Type",
    "PIP PEP",
    "Decoy Peptide",
    "Random RT"
  ),
  show_col_types = FALSE
)
organisms <- c("Homo sapiens" = "HUMAN", "Escherichia coli (strain K12)" = "ECOLI")
peptides <- peaks |>
  dplyr::filter(!.data$`Decoy Peptide`, !.data$`Random RT`, .data$Organism %in% names(organisms)) |>
  dplyr::mutate(confidence = ifelse(.data$`Peak Detection Type` == "MSMS", 1, 1 - .data$`PIP PEP`)) |>
  dplyr::summarise(
    intensity = sum(.data$`Peak intensity`),
    identification_confidence = max(.data$confidence),
    .by = c("File Name", "Protein Group", "Organism", "Base Sequence")
  ) |>
  dplyr::transmute(
    raw.file = .data$`File Name`,
    protein_Id = paste0(.data$`Protein Group`, "_", organisms[.data$Organism]),
    peptide = .data$`Base Sequence`,
    .data$intensity,
    .data$identification_confidence,
    dilution = tolower(sub(".*_human_ecoli_([A-E])_.*", "\\1", .data$raw.file))
  )
stopifnot(all(peptides$dilution %in% letters[1:5]))
prepare_peptides(peptides, c(dilution. = "dilution"), args$out_dir)
