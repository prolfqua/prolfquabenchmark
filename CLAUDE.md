# CLAUDE.md — prolfquabenchmark

## Purpose

Benchmarking vignettes comparing MSstats, proDA, msqrob2, and prolfqua using IonStar and CPTAC datasets.

## Build & Test

```bash
make check-fast    # R CMD check without vignettes (quick)
make check         # Full R CMD check including vignettes (~15 min)
make vignette V=BenchmarkingIonstarData  # Build a single vignette
```

## Dependency Management

DESCRIPTION is the single source of truth for all dependencies. If a package is missing, add it to `DESCRIPTION` in the appropriate field instead of relying on ad hoc local installs.

- **Imports:** packages used in vignettes and R code (prolfqua, tidyverse, QFeatures, etc.)
- **Suggests:** dev tooling (devtools, roxygen2, covr, lintr, etc.)
- **Remotes:** non-CRAN packages (prolfqua from GitHub, prolfquadata from GitLab)

Use the normal user / system R libraries for this workspace; `renv` autoload is disabled.

## Key Notes

- Vignettes are the main content — there is minimal R code outside of them
- prolfqua is on the `Modelling2R6` branch (specified in Remotes)
- The old `AnalysisTableAnnotation` class was merged into `AnalysisConfiguration` in prolfqua
