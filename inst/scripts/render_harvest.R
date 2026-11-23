#!/usr/bin/env Rscript
# Renders a benchmark vignette from before 2026-03-17, which keeps its results only as Benchmark objects, and
# writes every prolfqua Benchmark object left in the render environment as a contrast result directory
# (contrasts.tsv + metadata.yaml) that the score history collector reads. Used by prefill_score_history.sh.
#   Rscript render_harvest.R <vignette.Rmd> <dataset> <results_dir>
# Run in the vignette's directory, with the library of the old prolfqua first on the library path (R_LIBS).

args <- commandArgs(trailingOnly = TRUE)
vignette <- args[[1]]
dataset <- args[[2]]
results_dir <- args[[3]]

# The same analysis as a method benchmarked today: recorded under today's name, so its rows line up across releases
same_as_today <- c(
  benchmark_merged = "prolfqua_lm", # med. polish, lm moderated, merged with the imputed contrasts
  benchmark_prolfqua = "prolfqua_lm", # the same, FragPipe proteins
  benchmark_Limma = "prolfqua_limma"
)

env <- new.env()
rmarkdown::render(vignette, envir = env, quiet = TRUE)

objects <- ls(env)
benchmarks <- objects[vapply(objects, function(x) inherits(get(x, envir = env), "Benchmark"), logical(1))]
benchmarks <- benchmarks[!grepl("proDA|msstats|msqrob", benchmarks, ignore.case = TRUE)]
for (name in benchmarks) {
  bench <- get(name, envir = env)
  data <- as.data.frame(bench$data())
  # p-value and FDR columns as the Benchmark scaled and calibrated them (ROPECA: beta.based.significance); the average
  # intensity only feeds plot_scatter(), kept when the data has it
  columns <- c(
    protein_id = bench$hierarchy[[1]],
    contrast = bench$contrast,
    log_fc = bench$fcestimate,
    t_statistic = "statistic",
    p_value = bench$toscale[[1]],
    p_value_adjusted = bench$FDRvsFDP[[1]]$score
  )
  if (bench$avgInt %in% colnames(data)) {
    columns <- c(columns, avg_intensity = bench$avgInt)
  }
  absent <- setdiff(columns, colnames(data))
  if (length(absent) > 0) {
    message(name, ": not recorded, no column ", paste(absent, collapse = ", "))
    next
  }
  method <- if (name %in% names(same_as_today)) {
    same_as_today[[name]]
  } else {
    paste0("prolfqua_", sub("^benchmark_", "", name))
  }
  path <- file.path(results_dir, paste0(dataset, "_", method))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    stats::setNames(data[, columns], names(columns)),
    file.path(path, "contrasts.tsv"),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  yaml::write_yaml(
    list(
      dataset = dataset,
      method = method,
      method_description = bench$model_description,
      input_file = vignette,
      software_version = paste0("prolfqua ", utils::packageVersion("prolfqua")),
      date = format(Sys.Date()),
      ground_truth = list(
        id_column = "protein_id",
        positive = list(label = "ECOLI", pattern = "ECOLI"),
        negative = list(label = "HUMAN", pattern = "HUMAN")
      )
    ),
    file.path(path, "metadata.yaml")
  )
  message(name, ": recorded as ", method)
}
