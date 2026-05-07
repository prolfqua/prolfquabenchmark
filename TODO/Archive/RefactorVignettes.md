# Plan: Reduce prolfquabenchmark Vignette Runtime via Merging

## Context

The 11 vignettes in prolfquabenchmark take too long to build during `R CMD check`. The root cause: **6 vignettes** independently load, configure, filter, normalize, and aggregate the same IonStar/MaxQuant peptide dataset. By merging vignettes that use the same dataset and overlapping tools, we eliminate ~4 redundant full preprocessing cycles.

## Strategy: Merge vignettes (11 → 7)

No new helper functions. We merge 4 vignettes into existing ones that already share the same dataset and preprocessing. Each merge target already does the same data loading — we just add the extra model/benchmark sections from the absorbed vignette.

## Merges

### Merge 1: `Benchmark_rlm.Rmd` → `BenchmarkingIonstarData.Rmd`

**Rationale:** Both use IonStar/MQ peptides with identical preprocessing. `BenchmarkingIonstarData` already benchmarks 7 prolfqua model variants (lm, lm+mod, lm+DEqMS, missing, mixed effects, ROPECA, limma). RLM is just another prolfqua variant.

**What to move:**
- The `df.residual.rlm` and `sigma.rlm` helper definitions (lines 145-152 of Benchmark_rlm)
- RLM model fitting via `prolfqua::strategy_rlm()` (line 159)
- RLM contrast + benchmark creation (lines 223-246)
- RLM + moderation variant (lines 254-277)

**Where to insert:** After the existing model variants section in BenchmarkingIonstarData (around line 500+), add a new "## Robust Linear Model" section.

**Delete:** `vignettes/Benchmark_rlm.Rmd`

---

### Merge 2: `Benchmark_proDA_medpolish.Rmd` → `Benchmark_proDA_fromMQlfq.Rmd` (rename to `Benchmark_proDA.Rmd`)

**Rationale:** Both are proDA-only benchmarks on IonStar. One uses MQ peptides + median polish, the other uses MQ protein-level LFQ intensities. Combining them into one "proDA benchmark" vignette avoids loading the IonStar dataset twice.

**What to move:**
- The peptide-based preprocessing block from `_medpolish` (loading, config, filter, normalize, medpolish — lines 33-78)
- The proDA fit + contrasts on medpolished data (lines 97-113)
- The benchmark creation for medpolished proDA (lines 121-143)

**Structure of merged vignette:**
1. Section 1: "proDA on MQ protein-level LFQ" (existing `_fromMQlfq` content)
2. Section 2: "proDA on MQ peptides + median polish" (moved from `_medpolish`)
3. Section 3: Comparison of both approaches (combine benchmark plots)

**Rename:** `Benchmark_proDA_fromMQlfq.Rmd` → `Benchmark_proDA.Rmd`
**Delete:** `vignettes/Benchmark_proDA_medpolish.Rmd`
**Update:** VignetteIndexEntry in YAML header

---

### Merge 3: `Benchmark_MSStats.Rmd` → `BenchmarkFragPipeMSStats.Rmd`

**Rationale:** `BenchmarkFragPipeMSStats` already runs MSstats on FragPipe data. `Benchmark_MSStats` runs MSstats on MaxQuant evidence.txt. Both produce MSstats benchmarks. Combining avoids a separate vignette for a single-tool benchmark.

**What to move:**
- The MaxQuant evidence.txt loading + MSstats annotation setup (lines 38-63)
- `MaxQtoMSstatsFormat()` call (lines 68-73)
- Contaminant filtering (lines 80-82)
- `dataProcess()` + `groupComparison()` (lines 86-118)
- Benchmark creation (lines 126-168)

**Where to insert:** Add a new top-level section "# MSstats on MaxQuant evidence.txt" in BenchmarkFragPipeMSStats, after the existing FragPipe-based MSstats section. The final comparison section can then include both MSstats runs.

**Delete:** `vignettes/Benchmark_MSStats.Rmd`

---

### Merge 4: `BenchmarkMSqRob2.Rmd` → `BenchmarkFragPipeMSStats.Rmd`

**Rationale:** `BenchmarkFragPipeMSStats` already runs msqrob2 on FragPipe data. `BenchmarkMSqRob2` runs msqrob2 on MaxQuant peptides. But both do IonStar + msqrob2 hurdle. The MQ-based msqrob2 run is the more standard one (peptide.txt), so it's valuable to keep but as a section within the multi-tool comparison vignette.

**What to move:**
- MQ peptide loading + preprocessing (lines 51-98 — same as other IonStar vignettes)
- msqrob2 hurdle fitting + extraction (lines 121-227)
- Benchmark creation (lines 234-269)

**Complication:** BenchmarkFragPipeMSStats uses FragPipe data, not MQ data. Adding the MQ-based msqrob2 run means this vignette would load BOTH FragPipe AND MQ data. This adds some runtime back.

**Alternative:** Instead of merging into BenchmarkFragPipeMSStats, merge BenchmarkMSqRob2 into `BenchmarkingIonstarData.Rmd` which already loads MQ peptides. Add msqrob2 as an external-tool comparison alongside prolfqua's own variants.

**Decision:** Merge into `BenchmarkingIonstarData.Rmd` — it already has MQ peptide data loaded, so msqrob2 can run without extra preprocessing. Add a "## msqrob2 comparison" section.

**Delete:** `vignettes/BenchmarkMSqRob2.Rmd`

---

## Summary of changes

| Vignette | Action |
|----------|--------|
| `BenchmarkingIonstarData.Rmd` | **Absorbs** Benchmark_rlm + BenchmarkMSqRob2 |
| `Benchmark_proDA.Rmd` (renamed) | **Absorbs** Benchmark_proDA_medpolish |
| `BenchmarkFragPipeMSStats.Rmd` | **Absorbs** Benchmark_MSStats |
| `Benchmark_rlm.Rmd` | **Deleted** |
| `BenchmarkMSqRob2.Rmd` | **Deleted** |
| `Benchmark_proDA_medpolish.Rmd` | **Deleted** |
| `Benchmark_MSStats.Rmd` | **Deleted** |
| `Benchmark_cptac.Rmd` | Unchanged |
| `Benchmark_pipeline_demo.Rmd` | Unchanged |
| `Benchmark_Model_IonStar_With2Factors.Rmd` | Unchanged |
| `BenchmarkFragPipeProteinIonStar.Rmd` | Unchanged |

**Result: 7 vignettes** (down from 11), eliminating 4 redundant IonStar preprocessing cycles.

## Files to modify

1. `vignettes/BenchmarkingIonstarData.Rmd` — add RLM + msqrob2 sections
2. `vignettes/Benchmark_proDA_fromMQlfq.Rmd` — add medpolish section, rename to `Benchmark_proDA.Rmd`
3. `vignettes/BenchmarkFragPipeMSStats.Rmd` — add MQ-based MSstats section
4. Delete: `vignettes/Benchmark_rlm.Rmd`, `vignettes/BenchmarkMSqRob2.Rmd`, `vignettes/Benchmark_proDA_medpolish.Rmd`, `vignettes/Benchmark_MSStats.Rmd`

## Verification

1. `make vignette V=BenchmarkingIonstarData` — verify merged vignette renders
2. `make vignette V=Benchmark_proDA` — verify renamed/merged proDA vignette renders
3. `make vignette V=BenchmarkFragPipeMSStats` — verify merged MSstats vignette renders
4. `make check-fast` — quick check (no vignettes, confirms package structure is valid)
5. `make check` — full check, verify runtime is reduced and all vignettes pass
