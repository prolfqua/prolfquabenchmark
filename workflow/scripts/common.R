# Helpers shared by the pipeline scripts; each script sources this file. Scripts run from the package root.

# Positional command-line arguments, named
script_args <- function(names) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != length(names)) {
    stop("expected arguments: ", paste(names, collapse = " "), "; got ", length(args))
  }
  stats::setNames(as.list(args), names)
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

# Writes the contrasts of one dataset and method as a result directory (contrasts.tsv, metadata.yaml); `data` has
# prolfqua's column names (diff, statistic, p.value, FDR, avgInt, protein_Id) or the standard ones
write_result <- function(data, out_dir, spec, dataset, method, software_version) {
  prolfquabenchmark::write_contrast_results(
    data,
    out_dir,
    metadata = list(
      dataset = dataset,
      method = method,
      method_description = spec$method$description,
      input_file = spec$dataset$input_file,
      input = spec$method$input,
      level = spec$method$level,
      software_version = software_version,
      date = format(Sys.Date()),
      ground_truth = spec$dataset$ground_truth
    )
  )
}

# MaxQuant peptides.txt configuration shared by the MaxQuant datasets
create_mq_peptide_config <- function() {
  config <- prolfqua::AnalysisConfiguration$new()
  config$file_name <- "raw.file"
  config$hierarchy[["protein_Id"]] <- "leading.razor.protein"
  config$hierarchy[["peptide_Id"]] <- "sequence"
  config$hierarchy_depth <- 1
  config$ident_q_value <- "pep"
  config$set_response("peptide.intensity")
  config$min_peptides_protein <- 2
  config
}

# log2, then robust scaling on the human proteins
normalize_on_human <- function(lfq) {
  human <- lfq$get_copy()
  human$set_data(dplyr::filter(human$data_long(), grepl("HUMAN", .data$protein_Id)))
  human <- human$get_Transformer()$log2()$lfq
  lfq$get_Transformer()$log2()$robscale_subset(lfqsubset = human)$lfq
}

# Saves the prepared data of a dataset: proteins, and peptides unless the input is protein level
save_prepared <- function(out_dir, proteins, peptides = NULL) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(proteins, file.path(out_dir, "proteins.rds"))
  if (!is.null(peptides)) {
    saveRDS(peptides, file.path(out_dir, "peptides.rds"))
  }
}

quantdata_file <- function(name) {
  file.path(find.package("prolfquadata"), "quantdata", name)
}
