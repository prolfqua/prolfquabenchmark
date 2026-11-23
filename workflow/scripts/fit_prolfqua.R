# Fits one prolfqua method (a facade of build_contrast_analysis()) to one prepared dataset.
#   Rscript workflow/scripts/fit_prolfqua.R <config> <dataset> <method> <prepared_dir> <out_dir>
source("workflow/scripts/common.R")
args <- script_args(c("config", "dataset", "method", "prepared_dir", "out_dir"))
spec <- read_spec(args$config, args$dataset, args$method)

for (package in spec$method$packages) {
  library(package, character.only = TRUE)
}
lfq <- readRDS(file.path(args$prepared_dir, paste0(spec$method$input, ".rds")))
analysis <- do.call(
  prolfqua::build_contrast_analysis,
  c(
    list(
      lfqdata = lfq,
      modelstr = paste(spec$dataset$formula, spec$method$formula_terms),
      contrasts = spec$contrasts,
      method = spec$method$facade
    ),
    spec$method$args
  )
)
write_result(
  analysis$get_contrasts(),
  args$out_dir,
  spec,
  args$dataset,
  args$method,
  software_version = paste0("prolfqua ", utils::packageVersion("prolfqua"))
)
