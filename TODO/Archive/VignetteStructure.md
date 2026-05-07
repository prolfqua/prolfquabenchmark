# Vignette Structure Summary

## Overview

11 vignettes, ~3,880 lines total. Organized below by dataset and purpose.

---

## IonStar / MaxQuant peptides (6 vignettes, shared preprocessing)

These 6 vignettes all start from `MAXQuant_IonStar2018_PXD003881.zip` → `tidyMQ_Peptides()` → identical config/filter/normalize/aggregate pipeline (~100-150 lines each). This is the main source of redundancy.

### 1. `BenchmarkingIonstarData.Rmd` (840 lines) — **Comprehensive prolfqua benchmark**

- **Dataset:** IonStar / MQ peptides.txt
- **Preprocessing:** load → config (dilution factor) → filter (REV/CON, 2 peptides, small intensities) → log2 + robscale (HUMAN subset) → medpolish aggregation to protein level
- **Models benchmarked (all prolfqua):**
  - Linear model on proteins (`strategy_lm`)
  - Linear model + moderation (`ContrastsModerated`)
  - Linear model + DEqMS moderation (`ContrastsDEqMS`)
  - Missing data imputation (`ContrastsMissing`)
  - Mixed effects on peptides (`strategy_lmer` with peptide + sample random effects)
  - ROPECA (peptide-level lm, then aggregation)
  - Limma (`strategy_limma` + `build_model_limma`)
- **Contrasts:** 4 standard IonStar dilution contrasts
- **Output:** ROC/pAUC comparison of all 7 prolfqua model variants. Saves benchmark RDS files to `inst/Benchresults/`.

### 2. `Benchmark_rlm.Rmd` (278 lines) — **Robust linear model benchmark**

- **Dataset:** IonStar / MQ peptides.txt (identical preprocessing to #1)
- **Models benchmarked:**
  - Robust linear model (`strategy_rlm`) on protein level
  - RLM + moderation (`ContrastsModerated`)
- **Contrasts:** Defines all 10 possible IonStar contrasts but uses only the 4 standard ones
- **Output:** ROC, FDR vs FDP plots for RLM variants
- **Note:** Defines custom `df.residual.rlm()` and `sigma.rlm()` helpers. Has `SUBSET` flag for speedup.

### 3. `BenchmarkMSqRob2.Rmd` (297 lines) — **msqrob2 standalone benchmark**

- **Dataset:** IonStar / MQ peptides.txt (identical preprocessing to #1)
- **Models benchmarked:**
  - msqrob2 hurdle model (`msqrobHurdle` + `hypothesisTestHurdle`)
- **Preprocessing detail:** Creates QFeatures object, aggregates with custom `my_medianPolish()` function, runs hurdle model which combines intensity model (logFC) and count model (logOR) via anti_join preference logic.
- **Contrasts:** 4 standard IonStar dilution contrasts (encoded as msqrob2's `makeContrast` format)
- **Output:** ROC, FDR vs FDP, precision-recall. Saves benchmark RDS.
- **Dependencies:** QFeatures, msqrob2, limma, plotly (gated by `evalAll`)

### 4. `Benchmark_proDA_medpolish.Rmd` (186 lines) — **proDA on median-polished peptides**

- **Dataset:** IonStar / MQ peptides.txt (identical preprocessing to #1)
- **Models benchmarked:**
  - proDA on protein-level (after peptide → medpolish aggregation)
- **Contrasts:** 4 standard IonStar dilution contrasts via `proDA::test_diff()`
- **Output:** ROC, FDR vs FDP, precision-recall. Saves benchmark RDS.

### 5. `Benchmark_pipeline_demo.Rmd` (411 lines) — **File-based pipeline architecture demo**

- **Dataset:** IonStar / MQ peptides.txt (identical preprocessing to #1)
- **Models benchmarked:**
  - 4 prolfqua variants: lm, lm+moderated, missing, merged
- **Key difference:** Demonstrates the **file-based pipeline** workflow:
  - Fit → `write_contrast_results()` (writes contrasts.tsv + metadata.yaml)
  - Score → `benchmark_from_file()` (reads back and creates Benchmark)
  - Aggregate → `collect_benchmark_results()` + `plot_benchmark_comparison()`
- **Contrasts:** 4 standard IonStar dilution contrasts
- **Output:** Demonstrates decoupled fit/score/aggregate architecture. Comparison bar plots.

### 6. `Benchmark_Model_IonStar_With2Factors.Rmd` (260 lines) — **Two-factor interaction model**

- **Dataset:** IonStar / MQ peptides.txt, BUT with different factor design:
  - Drops sample "e" (only uses a, b, c, d)
  - Recodes samples into F1 (L1/L2) × F2 (L1/L2) factorial design
- **Preprocessing:** Same load/filter/normalize/aggregate pipeline but with 2 factors instead of dilution
- **Models benchmarked:**
  - `strategy_lm("abundance ~ F1. * F2.")` with interaction term
  - Moderated contrasts
  - Missing data imputation (`ContrastsMissing`)
  - Merged contrasts
- **Contrasts:** 5 custom contrasts including main effects, conditional effects, and interaction
- **Output:** ANOVA, volcano, MA plots. No benchmark/ROC (no ground truth evaluation).
- **Note:** More of a prolfqua tutorial than a benchmark. Has TODO comment for adding description.

---

## IonStar / MaxQuant protein-level LFQ (1 vignette)

### 7. `Benchmark_proDA_fromMQlfq.Rmd` (207 lines) — **proDA on MQ LFQ intensities**

- **Dataset:** IonStar / MQ proteinGroups.txt (protein-level LFQ intensities, NOT peptides)
- **Preprocessing:** `tidyMQ_ProteinGroups()` → manual config (AnalysisConfiguration) → filter REV/CON → filter nr.peptides > 1 → log2 + robscale (HUMAN subset). Different from peptide-based vignettes.
- **Models benchmarked:**
  - proDA on protein LFQ intensities
- **Contrasts:** 4 standard IonStar dilution contrasts via `proDA::test_diff()`
- **Output:** ROC, FDR vs FDP, precision-recall. Saves benchmark RDS.

---

## IonStar / MaxQuant evidence.txt (1 vignette)

### 8. `Benchmark_MSStats.Rmd` (188 lines) — **MSstats on MaxQuant evidence.txt**

- **Dataset:** IonStar / MQ evidence.txt + proteinGroups.txt (MSstats-native input format)
- **Preprocessing:** `MaxQtoMSstatsFormat()` → filter REV/CON → `MSstats::dataProcess()`
- **Models benchmarked:**
  - MSstats (`groupComparison()`)
- **Contrasts:** 4 standard IonStar dilution contrasts (via contrast matrix)
- **Output:** ROC, FDR vs FDP, precision-recall. Saves benchmark RDS.
- **Dependencies:** MSstats (gated by `evalAll`)

---

## IonStar / FragPipe (2 vignettes)

### 9. `BenchmarkFragPipeMSStats.Rmd` (618 lines) — **Multi-tool comparison on FragPipe data**

- **Dataset:** IonStar / FragPipe MSstats.csv (peptide-level, not MaxQuant)
- **Preprocessing:** Read FragPipe MSstats.csv → sum precursors to peptides → manual config → filter → log2 + robscale (HUMAN subset) → medpolish aggregation
- **Models benchmarked:**
  - prolfqua (lm + merged + moderated)
  - proDA (on aggregated proteins)
  - msqrob2 hurdle (on peptides via QFeatures)
  - MSstats (on FragPipe MSstats.csv directly)
- **Contrasts:** 4 standard IonStar dilution contrasts
- **Output:** Side-by-side benchmark comparison of all 4 tools. Bar plots, ROC tests via pROC::roc.test(). Saves benchmark RDS.

### 10. `BenchmarkFragPipeProteinIonStar.Rmd` (247 lines) — **prolfqua + proDA on FragPipe protein-level**

- **Dataset:** IonStar / FragPipe combined_protein.tsv (protein-level, not peptides)
- **Preprocessing:** `tidy_FragPipe_combined_protein_deprec()` → manual config → filter unique.stripped.peptides > 1 → log2 + robscale (HUMAN subset). No aggregation needed (already protein-level).
- **Models benchmarked:**
  - prolfqua (lm + merged + moderated)
  - proDA (on protein intensities)
- **Contrasts:** 4 standard IonStar dilution contrasts
- **Output:** Side-by-side prolfqua vs proDA comparison. Saves RDS.
- **Note:** Uses deprecated function `tidy_FragPipe_combined_protein_deprec()`

---

## CPTAC dataset (1 vignette)

### 11. `Benchmark_cptac.Rmd` (351 lines) — **Multi-tool comparison on CPTAC spike-in**

- **Dataset:** CPTAC (UPS spike-in in yeast background) from `msdata::quant()`
- **Preprocessing:** Read via `tidyMQ_Peptides()` → infer group from filename → filter → **VSN normalization** (not log2+robscale) → medpolish aggregation
- **Models benchmarked:**
  - prolfqua (lm + merged + moderated)
  - proDA (on aggregated proteins)
  - msqrob2 hurdle (on peptides via QFeatures)
- **Contrasts:** 1 contrast: `b_vs_a` (group.b - group.a)
- **Ground truth:** YEAST (background) vs UPS (spike-in), via `cptac_bench_preprocess()`
- **Output:** Side-by-side comparison of 3 tools. Benchmark plots.

---

## Redundancy Map

| Preprocessing step | Vignettes that repeat it |
|---|---|
| IonStar MQ peptide load + config + filter | #1, #2, #3, #4, #5, #6 |
| log2 + robscale (HUMAN subset) | #1, #2, #3, #4, #5, #6, #7, #9, #10 |
| medpolish aggregation | #1, #2, #3, #4, #5, #9, #11 |
| proDA fit + 4 IonStar contrasts | #4, #7, #9, #10 |
| msqrob2 hurdle + extraction logic | #3, #9, #11 |
| IonStar 4-contrast definition | #1, #2, #3, #4, #5, #7, #8, #9, #10 |
| Benchmark tail (pAUC + ROC + FDR + PR plots) | all 11 |

## Tool × Dataset Coverage Matrix

| Tool \ Dataset | IonStar/MQ pep | IonStar/MQ prot | IonStar/MQ evidence | IonStar/FragPipe pep | IonStar/FragPipe prot | CPTAC |
|---|---|---|---|---|---|---|
| **prolfqua lm** | #1, #5 | — | — | #9 | #10 | #11 |
| **prolfqua lm+mod** | #1, #5 | — | — | #9 | #10 | #11 |
| **prolfqua DEqMS** | #1 | — | — | — | — | — |
| **prolfqua missing** | #1, #5 | — | — | — | — | — |
| **prolfqua merged** | #5 | — | — | #9 | #10 | #11 |
| **prolfqua lmer** | #1 | — | — | — | — | — |
| **prolfqua ROPECA** | #1 | — | — | — | — | — |
| **prolfqua limma** | #1 | — | — | — | — | — |
| **prolfqua rlm** | #2 | — | — | — | — | — |
| **prolfqua 2-factor** | #6 | — | — | — | — | — |
| **proDA** | #4 | #7 | — | #9 | #10 | #11 |
| **msqrob2** | #3 | — | — | #9 | — | #11 |
| **MSstats** | — | — | #8 | #9 | — | — |
