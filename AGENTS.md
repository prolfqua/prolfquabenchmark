# CLAUDE.md — prolfquabenchmark

## Changelog

Record every change to this package's code or user-visible behaviour in `NEWS.md` at the repository root. Add a bullet
under the heading for the current `DESCRIPTION` version, creating a new version heading when a change opens one. Log the
user-visible effect, not the implementation detail, and match the existing format. Do not log changes that only touch
agent-instruction or repository meta files (`AGENTS.md`, `CLAUDE.md`, and similar) unless explicitly asked. `make
new-version` bumps the version, commits, and tags but does not write `NEWS.md`, so add the entry as part of the change
itself.

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
- prolfqua is on the `main` branch (specified in Remotes)
- The old `AnalysisTableAnnotation` class was merged into `AnalysisConfiguration` in prolfqua
