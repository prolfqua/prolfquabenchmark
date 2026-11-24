# CLAUDE.md — prolfquabenchmark

## Changelog

Record every change to this package's code or user-visible behaviour in `NEWS.md` at the repository root. Add a bullet
under the heading for the current `DESCRIPTION` version, creating a new version heading when a change opens one. Log the
user-visible effect, not the implementation detail, and match the existing format. Do not log changes that only touch
agent-instruction or repository meta files (`AGENTS.md`, `CLAUDE.md`, and similar) unless explicitly asked. `make
new-version` bumps the version, commits, and tags but does not write `NEWS.md`, so add the entry as part of the change
itself.

## Purpose

Benchmark of the prolfqua models against MSstats, proDA and msqrob2 on the IonStar and CPTAC spike-in datasets. A
Snakemake pipeline (`workflow/Snakefile`, `config/config.yaml`, scripts in `workflow/scripts`) runs the benchmark;
the vignettes `Benchmark_results` and `Score_history` only read the tables it exports to `inst/extdata`.

## Build & Test

```bash
make check-fast          # R CMD check without vignettes (quick)
make install             # the pipeline loads the installed prolfquabenchmark
make benchmark J=4       # run the pipeline, export inst/extdata/benchmark_*.tsv
make benchmark-dry       # list the jobs
make score-history       # append scores and step timings of the installed prolfqua (REPLACE=1 to rerun a version)
make score-history-release V=1.9.0 PROLFQUA=f4dfb7f3 PROLFQUAPP=18914707   # an earlier release, own library
make vignette V=Benchmark_results
```

- Add a method: an entry under `methods` in `config/config.yaml` with the `input` it needs; it runs on every dataset whose `inputs` contain all of them; a prolfqua facade needs no new script; prolfqua methods are named `prolfqua_<facade>`
- External tools take only their own input (feature table, peptides, the search software's `protein_table`), never the `proteins` prolfqua aggregated
- Add a dataset: one input file; an entry under `datasets` with its `inputs` and a `workflow/scripts/prepare_<dataset>.R` writing `proteins.rds` (and `peptides.rds` if it provides peptides, `protein_table.rds` if the search software reports proteins)
- Score history: starts at prolfqua 1.9.0; after each prolfqua tag, install it and run `make score-history`; an earlier release runs from `results/releases/<version>/lib` (fit and prepare scripts must not load prolfquabenchmark, which needs today's prolfqua)
- One name per dataset and method, used in the results, the exported tables and the score history; never a second label for the same input
- Python dependencies (snakemake) in the uv environment (`pyproject.toml`, `uv.lock`); R packages are the installed ones, declared in `DESCRIPTION`

## Dependency Management

DESCRIPTION is the single source of truth for all dependencies. If a package is missing, add it to `DESCRIPTION` in the appropriate field instead of relying on ad hoc local installs.

- **Imports:** packages used in vignettes and R code (prolfqua, tidyverse, QFeatures, etc.)
- **Suggests:** dev tooling (devtools, roxygen2, covr, lintr, etc.)
- **Remotes:** non-CRAN packages (prolfqua from GitHub, prolfquadata from GitLab)

Use the normal user / system R libraries for this workspace; `renv` autoload is disabled.

## Key Notes

- The analyses live in `workflow/scripts`; vignettes contain no analysis and render without data downloads
- prolfqua is on the `main` branch (specified in Remotes)
- The old `AnalysisTableAnnotation` class was merged into `AnalysisConfiguration` in prolfqua
