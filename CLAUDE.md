# CLAUDE.md — prolfquabenchmark

## Purpose

Benchmarking vignettes comparing MSstats, proDA, msqrob2, and prolfqua using IonStar and CPTAC datasets.

## Build & Test

```bash
make check-fast    # R CMD check without vignettes (quick)
make check         # Full R CMD check including vignettes (~15 min)
make vignette V=BenchmarkingIonstarData  # Build a single vignette
make renv-init     # First-time renv setup (installs all deps from DESCRIPTION)
make renv-restore  # Restore from lockfile
```

## Dependency Management

DESCRIPTION is the single source of truth for all dependencies. Never manually run `renv::install("pkg")` to fix missing packages — add them to DESCRIPTION (Imports or Suggests) instead, then `make renv-init`.

- **Imports:** packages used in vignettes and R code (prolfqua, tidyverse, QFeatures, etc.)
- **Suggests:** dev tooling (devtools, roxygen2, covr, lintr, etc.)
- **Remotes:** non-CRAN packages (prolfqua from GitHub, prolfquadata from GitLab)

The `renv/settings.json` includes `Suggests` in `package.dependency.fields` so dev tools are installed by renv.

## Key Notes

- Vignettes are the main content — there is minimal R code outside of them
- prolfqua is on the `Modelling2R6` branch (specified in Remotes)
- The old `AnalysisTableAnnotation` class was merged into `AnalysisConfiguration` in prolfqua
