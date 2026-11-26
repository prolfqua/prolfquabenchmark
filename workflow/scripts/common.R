# Helpers shared by the pipeline scripts; each script sources this file. Scripts run from the package root.

# Positional command-line arguments, named; with `rest`, the remaining arguments (at least one) as one element
script_args <- function(names, rest = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  n <- length(names)
  complete <- if (is.null(rest)) length(args) == n else length(args) > n
  if (!complete) {
    expected <- paste(c(names, if (!is.null(rest)) paste0(rest, "...")), collapse = " ")
    stop("expected arguments: ", expected, "; got ", length(args))
  }
  named <- stats::setNames(as.list(args[seq_len(n)]), names)
  if (is.null(rest)) named else c(named, stats::setNames(list(args[-seq_len(n)]), rest))
}

# The config entries of one dataset and one method, and the dataset's contrasts as a named character vector
read_spec <- function(config_file, dataset, method) {
  config <- yaml::read_yaml(config_file)
  dataset_spec <- config$datasets[[dataset]]
  list(
    dataset = dataset_spec,
    method = config$methods[[method]],
    contrasts = unlist(config$contrasts[[dataset_spec$contrasts]])
  )
}

# Writes the contrasts of one dataset and method as a result directory (contrasts.tsv, metadata.yaml), the format
# prolfquabenchmark::read_contrast_results() reads; `data` has prolfqua's column names (diff, statistic, p.value, FDR,
# avgInt, protein_Id) or the standard ones. Does not load prolfquabenchmark, which needs today's prolfqua: the fits of
# a release (make score-history-release) run with that release's prolfqua
write_result <- function(data, out_dir, spec, dataset, method, software_version) {
  standard <- c(
    diff = "log_fc",
    statistic = "t_statistic",
    p.value = "p_value",
    FDR = "p_value_adjusted",
    avgInt = "avg_intensity",
    protein_Id = "protein_id"
  )
  renamed <- names(data) %in% names(standard)
  names(data)[renamed] <- standard[names(data)[renamed]]
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.table(data, file.path(out_dir, "contrasts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  yaml::write_yaml(
    list(
      dataset = dataset,
      method = method,
      method_description = spec$method$description,
      input_file = spec$dataset$input_file,
      input = paste(spec$method$input, collapse = " + "),
      level = spec$method$level,
      software_version = software_version,
      date = format(Sys.Date()),
      ground_truth = spec$dataset$ground_truth
    ),
    file.path(out_dir, "metadata.yaml")
  )
}

# MaxQuant peptides.txt configuration shared by the MaxQuant datasets; the data needs identification_confidence,
# see mq_identification_confidence()
create_mq_peptide_config <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "raw.file"
  config$hierarchy[["protein_Id"]] <- "leading.razor.protein"
  config$hierarchy[["peptide_Id"]] <- "sequence"
  config$hierarchy_depth <- 1
  config$set_response("peptide.intensity")
  config$min_peptides_protein <- 2
  config
}

# MaxQuant peptides.txt has one PEP per peptide; its identification confidence is 1 - PEP, as in prolfquapp's reader
mq_identification_confidence <- function(peptides) {
  peptides$identification_confidence <- 1 - peptides$pep
  peptides
}

# log2, then robust scaling on the human proteins
normalize_on_human <- function(lfq) {
  human <- lfq$get_copy()
  human$set_data(dplyr::filter(human$data_long(), grepl("HUMAN", .data$protein_Id)))
  human <- human$get_Transformer()$log2()$lfq
  lfq$get_Transformer()$log2()$robscale_subset(lfqsubset = human)$lfq
}

# Saves the prepared data of a dataset: proteins; peptides unless the input is protein level; the protein table of
# the search software if the dataset has one
save_prepared <- function(out_dir, proteins, peptides = NULL, protein_table = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(proteins, file.path(out_dir, "proteins.rds"))
  if (!is.null(peptides)) {
    saveRDS(peptides, file.path(out_dir, "peptides.rds"))
  }
  if (!is.null(protein_table)) {
    saveRDS(protein_table, file.path(out_dir, "protein_table.rds"))
  }
}

quantdata_file <- function(name) {
  file.path(find.package("prolfquadata"), "quantdata", name)
}

# The file at `url`, downloaded once into `dir` under the last path component of the URL; stops unless its MD5 is
# `md5`
downloaded_file <- function(url, md5, dir) {
  dest <- file.path(dir, basename(sub("\\?.*$", "", url)))
  if (!file.exists(dest)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    old <- options(timeout = 3600)
    on.exit(options(old))
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  found <- unname(tools::md5sum(dest))
  if (!identical(found, md5)) {
    stop(dest, " has MD5 ", found, ", expected ", md5)
  }
  dest
}

# The member of `zip` matching `pattern`, extracted once into the folder `exdir` next to the ZIP
extracted_file <- function(zip, pattern, exdir) {
  member <- grep(pattern, utils::unzip(zip, list = TRUE)$Name, value = TRUE)
  path <- file.path(dirname(zip), exdir, member)
  if (!file.exists(path)) {
    utils::unzip(zip, files = member, exdir = file.path(dirname(zip), exdir))
  }
  path
}

# A file of the ProteoBench server (https://proteobench.cubimed.rub.de/<path>), cached in results/downloads/proteobench
proteobench_file <- function(path, md5) {
  downloaded_file(paste0("https://proteobench.cubimed.rub.de/", path), md5, "results/downloads/proteobench")
}

# The vendor table (input_file.*) of the ProteoBench submission `hash`, extracted next to its ZIP
proteobench_submission <- function(hash, md5) {
  zip <- proteobench_file(paste0("datasets/", hash, "/", hash, "_data.zip"), md5)
  extracted_file(zip, "^input_file\\.", hash)
}

# Species of each accession in a ProteoBench FASTA ZIP, the suffix of the UniProt entry name (sp|P12345|NAME_HUMAN);
# contaminants (accessions Cont_...) are left out
proteobench_species <- function(fasta_zip, md5) {
  zip <- proteobench_file(paste0("fasta/", fasta_zip), md5)
  fasta <- grep("^[^_].*\\.fasta$", utils::unzip(zip, list = TRUE)$Name, value = TRUE)
  headers <- grep("^>", readLines(unz(zip, fasta)), value = TRUE)
  ids <- strsplit(sub("^>(\\S+).*", "\\1", headers), "|", fixed = TRUE)
  species <- tibble::tibble(
    accession = vapply(ids, `[`, "", 2),
    species = sub(".*_", "", vapply(ids, `[`, "", 3))
  )
  dplyr::filter(species, !startsWith(.data$accession, "Cont_"))
}

# Prepares peptides of a spike-in dataset: `peptides` has one row per run and peptide with the columns raw.file,
# protein_Id (ending in the species, e.g. _HUMAN), peptide, intensity, identification_confidence and the factor column
# that `factor` names (c(<factor name> = <column>)). Proteins with at least 2 peptides, intensities below 4 removed,
# log2 and robust scaling on the human proteins, median polish to proteins.
prepare_peptides <- function(peptides, factor, out_dir) {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "raw.file"
  config$hierarchy[["protein_Id"]] <- "protein_Id"
  config$hierarchy[["peptide_Id"]] <- "peptide"
  config$hierarchy_depth <- 1
  config$set_response("intensity")
  config$factors[[names(factor)]] <- factor[[1]]
  config$factor_depth <- 1
  config$min_peptides_protein <- 2

  lfq <- prolfqua::LFQData$new(prolfqua::setup_analysis(peptides, config), config)
  lfq$filter_proteins_by_peptide_count()
  lfq$remove_small_intensities()
  normalized <- normalize_on_human(lfq)
  aggregator <- normalized$get_Aggregator()
  aggregator$aggregate()
  save_prepared(out_dir, proteins = aggregator$lfq_agg, peptides = normalized)
}

# Prepares a ProteoBench dataset: `peptides` has one row per run and peptide with the columns raw.file, protein
# (accessions separated by ";"), peptide, intensity and identification_confidence. A protein group whose accessions
# all have one species keeps it as suffix (protein_Id P12345;Q67890_HUMAN); groups with contaminants, unknown
# accessions or several species are removed. `conditions` names the condition of each run by a pattern of the run
# name (case-insensitive). Then prepare_peptides() with the factor group.
prepare_proteobench <- function(peptides, species, conditions, out_dir) {
  groups <- tibble::tibble(protein = unique(peptides$protein)) |>
    dplyr::mutate(accession = strsplit(.data$protein, ";", fixed = TRUE)) |>
    tidyr::unnest("accession") |>
    dplyr::left_join(species, by = "accession") |>
    dplyr::summarise(
      species = if (anyNA(.data$species) || dplyr::n_distinct(.data$species) > 1) NA else .data$species[1],
      .by = "protein"
    ) |>
    dplyr::filter(!is.na(.data$species))
  peptides <- dplyr::inner_join(peptides, groups, by = "protein")
  peptides$protein_Id <- paste0(peptides$protein, "_", peptides$species)
  peptides$condition <- NA_character_
  for (condition in names(conditions)) {
    peptides$condition[grepl(conditions[[condition]], peptides$raw.file, ignore.case = TRUE)] <- condition
  }
  if (anyNA(peptides$condition)) {
    stop("runs without a condition: ", paste(unique(peptides$raw.file[is.na(peptides$condition)]), collapse = ", "))
  }
  prepare_peptides(peptides, c(group. = "condition"), out_dir)
}

# DIA-NN report of a ProteoBench submission as peptides, read like prolfquapp's DIA-NN reader: protein groups at 1%
# (library and run), precursor quantities summed to stripped sequences, confidence 1 - the smallest precursor PEP
proteobench_diann_peptides <- function(file) {
  report <- prolfquapp::diann_read_output(arrow::read_parquet(file), Lib.PG.Q.Value = 0.01, PG.Q.Value = 0.01)
  diann <- prolfquapp::diann_output_to_peptide(report)
  tibble::tibble(
    raw.file = diann$raw.file,
    protein = diann$Protein.Group,
    peptide = diann$Stripped.Sequence,
    intensity = diann$Peptide.Quantity,
    identification_confidence = 1 - diann$PEP
  )
}
