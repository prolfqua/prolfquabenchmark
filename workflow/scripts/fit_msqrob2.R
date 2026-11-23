# Fits msqrob2's hurdle model to the prepared peptides of one dataset: peptides median-polished to proteins
# (QFeatures::aggregateFeatures), msqrobHurdle, hypothesisTestHurdle. Two ways to report it (method `combine`):
#   intensity_first  the intensity estimate where it exists, else the detection (count) estimate; BH per contrast
#   fisher           msqrob2's Fisher combination of both p-values (fisherPval, fisherAdjPval)
#   Rscript workflow/scripts/fit_msqrob2.R <config> <dataset> <method> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "method", "prepared_dir", "out_dir"))
spec <- read_spec(args$config, args$dataset, args$method)
suppressMessages({
  library(QFeatures)
  library(msqrob2)
})

peptides <- readRDS(file.path(args$prepared_dir, "peptides.rds"))
proteins <- readRDS(file.path(args$prepared_dir, "proteins.rds"))
factor_col <- all.vars(stats::as.formula(spec$dataset$formula))
levels <- sort(unique(as.character(peptides$data_long()[[factor_col]])))

# msqrob2 tests model parameters; with treatment coding a difference to the first level is that level's parameter
reference <- paste0(factor_col, levels[1])
hypotheses <- paste0(gsub(" ", "", sub(paste0(" - ", reference), "", spec$contrasts, fixed = TRUE)), "=0")
L <- makeContrast(hypotheses, parameterNames = paste0(factor_col, levels[-1]))
contrast_names <- stats::setNames(names(spec$contrasts), paste0("hurdle_", colnames(L)))

median_polish <- function(x, verbose = FALSE, ...) {
  fit <- stats::medpolish(x, na.rm = TRUE, trace.iter = verbose, maxiter = 10)
  fit$overall + fit$col
}
experiment <- prolfqua::lfqdata_to_summarized_experiment(lfqdata = peptides)
pe <- QFeatures::QFeatures(list(peptide = experiment), colData = colData(experiment))
pe <- QFeatures::aggregateFeatures(pe, i = "peptide", fcol = "protein_Id", name = "protein", fun = median_polish)
pe <- msqrobHurdle(pe, i = "protein", formula = stats::as.formula(spec$dataset$formula), overwrite = TRUE)
pe <- hypothesisTestHurdle(pe, i = "protein", L, overwrite = TRUE)

row_data <- rowData(pe[["protein"]])
tests <- dplyr::bind_rows(lapply(names(contrast_names), function(column) {
  tibble::as_tibble(as.data.frame(row_data[[column]]), rownames = "protein_Id") |>
    dplyr::mutate(contrast = contrast_names[[column]])
}))

result <- switch(
  spec$method$combine,
  intensity_first = {
    intensity <- !is.na(tests$logFCt)
    combined <- dplyr::transmute(
      tests,
      protein_Id = .data$protein_Id,
      contrast = .data$contrast,
      diff = ifelse(intensity, .data$logFC, .data$logOR),
      statistic = ifelse(intensity, .data$logFCt, .data$logORt),
      p.value = ifelse(intensity, .data$logFCpval, .data$logORpval)
    )
    dplyr::ungroup(prolfqua::adjust_p_values(combined, column = "p.value", group_by_col = "contrast"))
  },
  fisher = dplyr::transmute(
    tests,
    protein_Id = .data$protein_Id,
    contrast = .data$contrast,
    diff = dplyr::coalesce(.data$logFC, .data$logOR),
    statistic = dplyr::coalesce(.data$logFCt, .data$logORt),
    p.value = .data$fisherPval,
    FDR = .data$fisherAdjPval
  ),
  stop("unknown combine ", spec$method$combine)
)

# Average intensity per protein and contrast: the mean of the protein means of the groups in the contrast
group_means <- proteins$get_Stats()$stats()
group_means$group <- paste0(factor_col, group_means[[factor_col]])
average <- dplyr::bind_rows(lapply(names(spec$contrasts), function(name) {
  groups <- all.vars(parse(text = spec$contrasts[[name]]))
  dplyr::filter(group_means, .data$group %in% groups) |>
    dplyr::summarize(avgInt = mean(.data$meanAbundance, na.rm = TRUE), .by = "protein_Id") |>
    dplyr::mutate(contrast = name)
}))
# Rows by contrast and protein: the scores break ties (detection estimates share values) by row order
result <- dplyr::left_join(result, average, by = c("protein_Id", "contrast")) |>
  dplyr::arrange(.data$contrast, .data$protein_Id)

write_result(
  result,
  args$out_dir,
  spec,
  args$dataset,
  args$method,
  software_version = paste0("msqrob2 ", utils::packageVersion("msqrob2"))
)
