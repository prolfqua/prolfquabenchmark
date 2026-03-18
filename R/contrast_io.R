#' @importFrom dplyr mutate filter case_when
#' @importFrom rlang sym .data
#' @importFrom utils read.delim write.table
NULL

# Column mapping: prolfqua-native → standardized snake_case
.col_map_to_standard <- c(
  "diff" = "log_fc",
  "statistic" = "t_statistic",
  "p.value" = "p_value",
  "FDR" = "p_value_adjusted",
  "avgInt" = "avg_intensity",
  "protein_Id" = "protein_id"
)

.col_map_from_standard <- stats::setNames(
  names(.col_map_to_standard),
  .col_map_to_standard
)

.required_standard_cols <- c("protein_id", "contrast", "log_fc", "t_statistic",
                             "p_value", "p_value_adjusted")

.required_metadata_fields <- c("dataset", "method", "method_description",
                               "input_file", "software_version", "date",
                               "ground_truth")


#' Write contrast results to a standardized file format
#'
#' Writes a contrasts.tsv file and a metadata.yaml sidecar file to the
#' specified directory. Accepts prolfqua-native column names and renames
#' them to the standardized snake_case format.
#'
#' @param data data.frame with contrast results (one row per protein x contrast)
#' @param path directory to write to (created if it does not exist)
#' @param metadata list with required fields: \code{dataset}, \code{method},
#'   \code{method_description}, \code{input_file}, \code{software_version},
#'   \code{date}, and \code{ground_truth}
#' @param col_map named character vector mapping input column names to standard
#'   names. Defaults to the prolfqua-native mapping. Names are input columns,
#'   values are output columns.
#' @return Invisibly returns the path written to.
#' @export
#' @family benchmarking
#' @examples
#' \dontrun{
#' metadata <- list(
#'   dataset = "ionstar_maxquant",
#'   method = "prolfqua_lm_mod",
#'   method_description = "Linear model with variance moderation",
#'   input_file = "MAXQuant_IonStar2018_PXD003881.zip/evidence.txt",
#'   input_level = "peptide",
#'   aggregation = "median_polish",
#'   normalization = "robscale",
#'   software_version = paste0("prolfqua ", packageVersion("prolfqua")),
#'   date = Sys.Date(),
#'   contrasts = c("dilution_(4.5/3)_1.5", "dilution_(6/4.5)_1.3(3)"),
#'   ground_truth = list(
#'     id_column = "protein_id",
#'     positive = list(label = "ECOLI", pattern = "ECOLI"),
#'     negative = list(label = "HUMAN", pattern = "HUMAN")
#'   )
#' )
#' write_contrast_results(contrast_data, "/tmp/ionstar_lm", metadata)
#' }
write_contrast_results <- function(data, path, metadata,
                                   col_map = .col_map_to_standard) {
  # Validate required metadata fields
  missing_meta <- setdiff(.required_metadata_fields, names(metadata))
  if (length(missing_meta) > 0) {
    stop("Missing required metadata fields: ",
         paste(missing_meta, collapse = ", "))
  }

  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  # Rename columns that exist in the mapping
  data_out <- data
  present <- intersect(names(col_map), colnames(data_out))
  for (old_name in present) {
    colnames(data_out)[colnames(data_out) == old_name] <- col_map[[old_name]]
  }

  # Validate required columns
  missing_cols <- setdiff(.required_standard_cols, colnames(data_out))
  if (length(missing_cols) > 0) {
    stop("Missing required columns after mapping: ",
         paste(missing_cols, collapse = ", "),
         "\nPresent columns: ", paste(colnames(data_out), collapse = ", "))
  }

  # Write TSV
  tsv_path <- file.path(path, "contrasts.tsv")
  utils::write.table(data_out, file = tsv_path, sep = "\t",
                     row.names = FALSE, quote = FALSE)

  # Write metadata YAML

  yaml_path <- file.path(path, "metadata.yaml")
  yaml::write_yaml(metadata, file = yaml_path)

  message("Wrote contrast results to: ", path)
  invisible(path)
}


#' Read contrast results from standardized file format
#'
#' Reads contrasts.tsv and metadata.yaml from the specified directory.
#' Applies ground truth annotation (species and TP columns) based on the
#' ground_truth specification in metadata.yaml.
#'
#' @param path directory containing contrasts.tsv and metadata.yaml
#' @return list with elements:
#'   \describe{
#'     \item{data}{data.frame with contrast results including species and TP columns}
#'     \item{metadata}{list parsed from metadata.yaml}
#'   }
#' @export
#' @family benchmarking
#' @examples
#' \dontrun{
#' res <- read_contrast_results("/tmp/ionstar_lm")
#' head(res$data)
#' res$metadata$method
#' }
read_contrast_results <- function(path) {
  tsv_path <- file.path(path, "contrasts.tsv")
  yaml_path <- file.path(path, "metadata.yaml")

  if (!file.exists(tsv_path)) {
    stop("contrasts.tsv not found in: ", path)
  }
  if (!file.exists(yaml_path)) {
    stop("metadata.yaml not found in: ", path)
  }

  data <- utils::read.delim(tsv_path, sep = "\t", stringsAsFactors = FALSE)
  meta <- yaml::read_yaml(yaml_path)

  # Validate required columns
  missing_cols <- setdiff(.required_standard_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("contrasts.tsv missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Apply ground truth annotation from metadata
  gt <- meta$ground_truth
  if (!is.null(gt)) {
    id_col <- gt$id_column
    if (!id_col %in% colnames(data)) {
      stop("ground_truth$id_column '", id_col, "' not found in data")
    }
    data <- data |>
      dplyr::mutate(
        species = dplyr::case_when(
          grepl(gt$positive$pattern, !!rlang::sym(id_col)) ~ gt$positive$label,
          grepl(gt$negative$pattern, !!rlang::sym(id_col)) ~ gt$negative$label,
          TRUE ~ "OTHER"
        )
      ) |>
      dplyr::filter(.data$species != "OTHER") |>
      dplyr::mutate(TP = (.data$species == gt$positive$label))
  }

  list(data = data, metadata = meta)
}
