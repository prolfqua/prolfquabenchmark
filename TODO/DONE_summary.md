# Changes Summary — 2026-03-18 to 2026-03-30

## Overview

Completed two major refactoring efforts:
1. **Vignette consolidation** (11 → 6) by organizing vignettes per tool instead of per dataset
2. **Pipeline interface** (`write_contrast_results` → `benchmark_from_result` workflow)

Both plans from `TODO/` are now fully addressed:
- `RefactorVignettes.md` — vignette merging (achieved 6 instead of the planned 7)
- `snuggly-drifting-mccarthy.md` — standardized contrast file format & pipeline architecture

## Vignette Consolidation (11 → 6)

The strategy shifted from the original merge plan (by dataset) to organizing **by tool**, which produced fewer, more focused vignettes.

### Final 6 vignettes

| Vignette | Tool(s) | Datasets |
|----------|---------|----------|
| `Benchmark_prolfqua.Rmd` | prolfqua (missing + limma) | IonStar/MQ, FragPipe MSstats, FragPipe protein, CPTAC |
| `Benchmark_msqrob2.Rmd` | msqrob2 hurdle | IonStar/MQ, IonStar/FragPipe, CPTAC |
| `Benchmark_proDA.Rmd` | proDA | IonStar MQ LFQ, IonStar MQ medpolish |
| `Benchmark_MSStats.Rmd` | MSstats | IonStar/MQ, IonStar/FragPipe |
| `Benchmark_pipeline_demo.Rmd` | prolfqua (demo) | IonStar/MQ |
| `Benchmark_Model_IonStar_With2Factors.Rmd` | prolfqua (2-factor) | IonStar/MQ |

### Deleted vignettes (absorbed into the above)

- `BenchmarkingIonstarData.Rmd` → `Benchmark_prolfqua.Rmd`
- `Benchmark_rlm.Rmd` → `Benchmark_prolfqua.Rmd`
- `BenchmarkMSqRob2.Rmd` → `Benchmark_msqrob2.Rmd`
- `Benchmark_proDA_medpolish.Rmd` → `Benchmark_proDA.Rmd`
- `Benchmark_proDA_fromMQlfq.Rmd` → `Benchmark_proDA.Rmd`
- `BenchmarkFragPipeMSStats.Rmd` / `BenchmarkFragPipe_MSstatsFormat_prolfqua_msqrob2.Rmd` → split across tool vignettes
- `BenchmarkFragPipeProteinIonStar.Rmd` / `BenchmarkFragPipe_combinedProtein_prolfqua.Rmd` → `Benchmark_prolfqua.Rmd`
- `Benchmark_cptac.Rmd` / `Benchmark_cptac_prolfqua_msqrob2.Rmd` → split across tool vignettes

## R Code Changes

### `R/benchmark_io.R`
- **`benchmark_from_result()`** (new) — convenience factory that takes the list returned by `write_contrast_results()`, applies ground truth annotation, and creates a `Benchmark` object. Eliminates repeated manual `make_benchmark()` calls in vignettes.

### `R/contrast_io.R`
- **`annotate_ground_truth()`** (new, internal) — extracted from `read_contrast_results()` for reuse by both `read_contrast_results()` and `benchmark_from_result()`.
- **`write_contrast_results()`** — now returns `list(data = data_out, metadata = metadata)` instead of just the path. Enables direct piping into `benchmark_from_result()`.

### `R/Benchmark.R`
- Fixed roxygen: moved `.geomean` helper above `ms_bench_ap` docs block so documentation attaches to the correct function.

### Vignette improvements
- msqrob2 vignette: extracted shared `run_msqrob2_hurdle()` helper, uses `benchmark_from_result()` workflow
- Added `avgAbd` → `avg_intensity` to msqrob2 column map
- Added missing-package notice blocks to MSStats, proDA, and prolfqua vignettes (user-friendly message when optional deps are missing)
- Fixed `Benchmark_pipeline_demo.Rmd`: added missing `input_file` metadata field

## Pipeline Architecture (fully implemented)

The full pipeline from `snuggly-drifting-mccarthy.md` is now in place:

```
fit → write_contrast_results() → contrasts.tsv + metadata.yaml
                ↓ (returns list)
        benchmark_from_result()  ← in-memory path
                or
        benchmark_from_file()    ← file-based path
                ↓
        write_benchmark_results() → benchmark_results.tsv
                ↓
        collect_benchmark_results() + plot_benchmark_comparison()
```

## TODO documents archived

The following planning documents are now complete and can be archived:
- `RefactorVignettes.md`
- `VignetteStructure.md`
- `snuggly-drifting-mccarthy.md`
