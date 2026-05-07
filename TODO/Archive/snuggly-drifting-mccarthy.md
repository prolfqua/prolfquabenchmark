# prolfquabenchmark: Inventory & Interface Design Plan

## Context

The goal is to decouple **model fitting** from **benchmark scoring** in prolfquabenchmark by defining a clear interchange file format. Today, each vignette loads raw data, fits models, creates `Benchmark` objects, and produces plots — all in one monolithic script. The vision is a pipeline architecture where:

1. A fitting step produces a **standardized contrast result file** (per dataset × method)
2. A scoring step loads that file into the `Benchmark` class and computes all metrics
3. A reporting step collects benchmark results across methods and produces summary tables/plots

This enables a future Snakemake pipeline: `for each dataset × method → fit → score → aggregate`.

---

## 1. Inventory: Datasets

| Dataset | Source | Organism | Design | Data Levels in Vignettes |
|---------|--------|----------|--------|--------------------------|
| **IonStar (MaxQuant)** | PXD003881 via `prolfquadata` | HUMAN background + ECOLI spike-in | 5 dilutions (3, 4.5, 6, 7.5, 9), replicates per dilution | Peptide→Protein (median polish), Protein (MQ LFQ) |
| **IonStar (FragPipe)** | PXD003881, FragPipe v14 | Same | Same | Protein (`combined_protein.tsv`), Precursor/Peptide (`MSstats.csv`) |
| **CPTAC Study 6** | `msdata` R package | YEAST background + UPS1 spike-in | 2 concentrations (0.25, 0.74 fmol), 3 reps each | Peptide→Protein (median polish + VSN) |

**Ground truth assignment:**
- IonStar: TP = `ECOLI` in protein_Id, TN = `HUMAN`
- CPTAC: TP = `UPS` in protein_Id, TN = `YEAST`

---

## 2. Inventory: Contrasts

### IonStar — 10 possible, 4 used in benchmarks (smallest FCs)

| Contrast | Comparison | Expected log2 FC |
|----------|------------|-------------------|
| `dilution_(4.5/3)_1.5` | b − a | log2(1.5) ≈ 0.58 |
| `dilution_(6/4.5)_1.3(3)` | c − b | log2(1.33) ≈ 0.41 |
| `dilution_(7.5/6)_1.25` | d − c | log2(1.25) ≈ 0.32 |
| `dilution_(9/7.5)_1.2` | e − d | log2(1.2) ≈ 0.26 |

### CPTAC — 1 contrast

| Contrast | Comparison |
|----------|------------|
| `b_vs_a` | 0.74 fmol vs 0.25 fmol UPS1 |

---

## 3. Inventory: Methods Benchmarked (across vignettes)

| Vignette | Method(s) | Input Level |
|----------|-----------|-------------|
| `BenchmarkingIonstarData` | prolfqua: missing, lm, lm_mod, lm_DEqMS, mix_eff, mix_eff_mod, ropeca | Peptide→Protein |
| `Benchmark_rlm` | prolfqua: rlm, rlm_mod | Peptide→Protein |
| `Benchmark_MSStats` | MSstats | Peptide |
| `Benchmark_proDA_fromMQlfq` | proDA | Protein (MQ LFQ) |
| `Benchmark_proDA_medpolish` | proDA | Peptide→Protein |
| `BenchmarkMSqRob2` | msqrob2 (msqrobHurdle) | Peptide→Protein (QFeatures) |
| `BenchmarkFragPipeProteinIonStar` | prolfqua lm | Protein (FragPipe) |
| `BenchmarkFragPipeMSStats` | prolfqua, proDA | Precursor/Peptide (FragPipe) |
| `Benchmark_cptac` | prolfqua merged, proDA | Peptide→Protein (CPTAC) |
| `Benchmark_Model_IonStar_With2Factors` | prolfqua (2-factor) | Peptide/Protein (WIP) |

---

## 4. Inventory: Benchmark Scores

The `Benchmark` class ranks proteins by **3 score columns** (configurable), then computes confusion-matrix metrics cumulatively:

### Score columns used for ranking

| Score | Direction | Meaning |
|-------|-----------|---------|
| `diff` (or `logFC`, `log2FC`) | desc | Fold-change estimate |
| `statistic` (or `Tvalue`, `t_statistic`) | desc | Test statistic |
| `scaled.p.value` (or `scaled.pvalue`, `scaled.pval`) | desc | P-value × sign(FC), computed internally from `toscale` columns |

### FDR calibration score
| Score | Direction |
|-------|-----------|
| `FDR` (or `adj.pvalue`, `adj_pval`) | asc |

### Derived metrics (per score, per contrast)

| Metric | Formula |
|--------|---------|
| TPR (Recall) | cumsum(TP) / total_TP |
| FPR | cumsum(!TP) / total_TN |
| Precision | cumsum(TP) / (cumsum(TP) + cumsum(!TP)) |
| FDP | cummean(!TP) |
| ACC | (TP_hits + TN_hits) / N |

### Summary metrics

| Metric | Meaning |
|--------|---------|
| AUC | Full area under ROC |
| pAUC_10 | Partial AUC at FPR ≤ 0.1 |
| pAUC_20 | Partial AUC at FPR ≤ 0.2 |

### Plots

ROC curves, Precision-Recall, FDR vs FDP, score distributions (ridgeline), intensity vs score scatter — all faceted by contrast.

---

## 5. Current Interface (as-is)

The current `Benchmark$new()` constructor expects a **data.frame** with:

**Required columns:**
- `protein_Id` (or custom `hierarchy`)
- `contrast` — contrast name
- `TP` — boolean, pre-annotated by `ionstar_bench_preprocess()` or `cptac_bench_preprocess()`
- `species` — "HUMAN"/"ECOLI" or "YEAST"/"UPS"
- `diff` — fold-change estimate
- `statistic` — test statistic
- `p.value` — raw p-value (gets internally scaled)
- `FDR` — adjusted p-value
- `avgInt` — average intensity

**Current serialization:** ad-hoc `saveRDS()` of Benchmark objects or contrast lists into `inst/Benchresults/`.

---

## 6. Proposed Design: Standardized Contrast Result File Format

### 6.1 Fitting Output File (TSV + metadata YAML)

Each fitting method produces **one directory per dataset × method** containing:

**`contrasts.tsv`** — long format, one row per protein × contrast, **tool-agnostic snake_case columns**:

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `protein_id` | char | yes | Protein identifier |
| `contrast` | char | yes | Contrast name |
| `log_fc` | numeric | yes | log2 fold-change estimate |
| `t_statistic` | numeric | yes | Test statistic (t, z, etc.) |
| `p_value` | numeric | yes | Raw p-value |
| `p_value_adjusted` | numeric | yes | BH-adjusted p-value |
| `avg_intensity` | numeric | recommended | Mean log2 intensity |
| `df` | numeric | optional | Degrees of freedom |
| `sigma` | numeric | optional | Residual standard error |
| `n_obs` | integer | optional | Number of observations used |
| `n_missing` | integer | optional | Number of missing values |

**Column mapping from current code:** The `write_contrast_results()` function will accept prolfqua-native names and rename:
- `diff` → `log_fc`, `statistic` → `t_statistic`, `p.value` → `p_value`, `FDR` → `p_value_adjusted`, `avgInt` → `avg_intensity`

**`metadata.yaml`** — fitting provenance + ground truth specification:

```yaml
# Fitting provenance
dataset: ionstar_maxquant
method: prolfqua_lm_mod
method_description: "Linear model with variance moderation on median-polished proteins"
input_level: peptide          # peptide, precursor, protein
aggregation: median_polish    # none, median_polish, topN, sum
normalization: robscale       # robscale, vsn, quantile, none
software_version: "prolfqua 4.5.0"
date: "2026-03-17"
contrasts:
  - dilution_(4.5/3)_1.5
  - dilution_(6/4.5)_1.3(3)
  - dilution_(7.5/6)_1.25
  - dilution_(9/7.5)_1.2

# Ground truth specification (enables automated TP annotation)
ground_truth:
  id_column: protein_id
  positive:
    label: ECOLI             # species label for true positives
    pattern: "ECOLI"         # regex pattern to match in protein_id
  negative:
    label: HUMAN             # species label for true negatives
    pattern: "HUMAN"         # regex pattern to match in protein_id
  # proteins matching neither pattern are filtered out as "OTHER"
```

For CPTAC, the ground_truth block would be:
```yaml
ground_truth:
  id_column: protein_id
  positive:
    label: UPS
    pattern: "UPS"
  negative:
    label: YEAST
    pattern: "YEAST"
```

### 6.2 Loading into Benchmark

The ground truth annotation is now **data-driven from metadata.yaml** — no need to pick a preprocessing function:

```r
# Read contrasts + metadata, annotate TP/species from ground_truth spec
load_contrast_results <- function(path) {
  data <- read_tsv(file.path(path, "contrasts.tsv"))
  meta <- yaml::read_yaml(file.path(path, "metadata.yaml"))

  # Automated ground truth annotation from metadata
  gt <- meta$ground_truth
  data <- data |>
    mutate(species = case_when(
      grepl(gt$positive$pattern, !!sym(gt$id_column)) ~ gt$positive$label,
      grepl(gt$negative$pattern, !!sym(gt$id_column)) ~ gt$negative$label,
      TRUE ~ "OTHER"
    )) |>
    filter(species != "OTHER") |>
    mutate(TP = (species == gt$positive$label))

  list(data = data, metadata = meta)
}
```

Convenience constructor:

```r
benchmark_from_file <- function(path, ...) {
  res <- load_contrast_results(path)
  # Map tool-agnostic column names to Benchmark expectations
  make_benchmark(
    prpr = res$data,
    model_name = res$metadata$method,
    model_description = res$metadata$method_description,
    fcestimate = "log_fc",
    toscale = c("p_value"),
    hierarchy = c("protein_id"),
    benchmark = list(
      list(score = "log_fc", desc = TRUE),
      list(score = "t_statistic", desc = TRUE),
      list(score = "scaled.p_value", desc = TRUE)  # generated by Benchmark internally
    ),
    FDRvsFDP = list(list(score = "p_value_adjusted", desc = FALSE)),
    ...
  )
}
```

### 6.3 Benchmark Results Output (serialized summary)

After scoring, produce a **results file** that captures metrics without the full confusion curve data:

**`benchmark_results.tsv`** — one row per method × contrast × score:

| Column | Type | Description |
|--------|------|-------------|
| `model_name` | char | Method identifier |
| `model_description` | char | Human-readable description |
| `dataset` | char | Dataset identifier |
| `contrast` | char | Contrast name |
| `score` | char | Score column used (diff, statistic, scaled.p.value) |
| `AUC` | numeric | Full AUC |
| `pAUC_10` | numeric | Partial AUC at FPR ≤ 0.1 |
| `pAUC_20` | numeric | Partial AUC at FPR ≤ 0.2 |
| `n_TP` | integer | Number of true positives |
| `n_TN` | integer | Number of true negatives |
| `n_total` | integer | Total proteins scored |
| `n_missing_contrasts` | integer | Proteins with missing contrasts |

This table is the key input for cross-method summary visualization.

Add to Benchmark class:

```r
# Method to export summary
Benchmark$set("public", "to_summary_table", function(dataset = "unknown") { ... })

# Standalone function to collect across methods
collect_benchmark_results <- function(result_files) {
  # Read all benchmark_results.tsv files
  # Bind into single table
  # Return for visualization
}
```

---

## 7. Implementation Steps

### Step 1: Define and document the contrast result file format
- Create `R/contrast_io.R` with `write_contrast_results()` and `read_contrast_results()`
- TSV for data, YAML sidecar for metadata
- `write_contrast_results()` accepts prolfqua-native column names and renames to snake_case
- `read_contrast_results()` validates required columns, reads metadata, applies ground truth annotation from `metadata.yaml$ground_truth`
- Validate required columns on read

### Step 2: Add `benchmark_from_file()` factory
- In `R/Benchmark.R` or new `R/benchmark_io.R`
- Wraps `read_contrast_results()` (which now handles ground truth) + `make_benchmark()`
- Maps snake_case column names to Benchmark constructor params

### Step 3: Add `to_summary_table()` method to Benchmark
- Extracts pAUC metrics + counts into a flat data.frame
- Add `write_benchmark_results()` / `read_benchmark_results()` for TSV serialization

### Step 4: Add `collect_benchmark_results()` aggregator
- Reads multiple result files
- Produces combined summary table
- Add basic comparison plots (grouped bar charts of AUC by method × contrast)

### Step 5: Refactor one vignette as proof of concept
- Split `BenchmarkingIonstarData.Rmd` into fitting script + scoring script
- Fitting script produces `contrasts.tsv` + `metadata.yaml`
- Scoring script loads and benchmarks

### Step 6: (Future) Snakemake pipeline skeleton
- Not in this PR, but the file format enables it

---

## 8. Files to Create/Modify

| File | Action | Purpose |
|------|--------|---------|
| `R/contrast_io.R` | **create** | `write_contrast_results()`, `read_contrast_results()` |
| `R/benchmark_io.R` | **create** | `benchmark_from_file()`, `write_benchmark_results()`, `read_benchmark_results()`, `collect_benchmark_results()` |
| `R/Benchmark.R` | **modify** | Add `to_summary_table()` method to Benchmark class |
| `DESCRIPTION` | **modify** | Add `yaml` to Imports |
| One vignette | **modify** | Refactor as proof of concept |

---

## 9. Verification

- `make document` succeeds after changes
- `make check-fast` passes
- Refactored vignette produces identical benchmark results (compare pAUC values)
- Round-trip test: write contrast results → read → benchmark → same AUC as direct path
