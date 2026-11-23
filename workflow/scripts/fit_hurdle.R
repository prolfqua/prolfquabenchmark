# Fits a hurdle method to one prepared dataset: an intensity facade and a detection facade, combined by
# ContrastsHurdleFacade (Fisher's method where both p-values exist). With `combine: fisher_always` both p-values are
# combined whenever they exist, the msqrob2 rule.
#   Rscript workflow/scripts/fit_hurdle.R <config> <dataset> <method> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "method", "prepared_dir", "out_dir"))
spec <- read_spec(args$config, args$dataset, args$method)

fit_component <- function(component) {
  lfq <- readRDS(file.path(args$prepared_dir, paste0(component$input, ".rds")))
  do.call(
    prolfqua::build_contrast_analysis,
    c(
      list(lfqdata = lfq, modelstr = spec$dataset$formula, contrasts = spec$contrasts, method = component$facade),
      component$args
    )
  )
}

# Fisher's method on the intensity and detection p-values whenever at least one exists
fisher_always <- function(hurdle) {
  log_p <- cbind(log(hurdle$p.value.intensity), log(hurdle$p.value.detection))
  n_components <- rowSums(!is.na(log_p))
  chisq <- -2 * rowSums(log_p, na.rm = TRUE)
  hurdle$p.value <- ifelse(n_components > 0, stats::pchisq(chisq, 2 * n_components, lower.tail = FALSE), NA_real_)
  dplyr::ungroup(prolfqua::adjust_p_values(hurdle))
}

hurdle <- prolfqua::ContrastsHurdleFacade$new(
  fit_component(spec$method$intensity),
  fit_component(spec$method$detection)
)$get_contrasts()
if (identical(spec$method$combine, "fisher_always")) {
  hurdle <- fisher_always(hurdle)
}
write_result(
  hurdle,
  args$out_dir,
  spec,
  args$dataset,
  args$method,
  software_version = paste0("prolfqua ", utils::packageVersion("prolfqua"))
)
