# prolfquabenchmark

Benchmarking vignettes comparing differential expression analysis methods — **prolfqua**, **MSstats**, **proDA**, **msqrob2**, **limma**, and **DEqMS** — using the IonStar and CPTAC gold-standard datasets.

Full documentation and rendered vignettes: **https://prolfqua.github.io/prolfquabenchmark/**

## Vignettes

| Dataset | Vignette |
|---------|----------|
| IonStar / MaxQuant | [Benchmarking normalization, aggregation and models](https://prolfqua.github.io/prolfquabenchmark/articles/BenchmarkingIonstarData.html) |
| IonStar / MaxQuant | [MSstats](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_MSStats.html) |
| IonStar / MaxQuant | [proDA (LFQ intensities)](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_proDA_fromMQlfq.html) |
| IonStar / MaxQuant | [proDA (peptides)](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_proDA_medpolish.html) |
| IonStar / MaxQuant | [msqrob2](https://prolfqua.github.io/prolfquabenchmark/articles/BenchmarkMSqRob2.html) |
| IonStar / MaxQuant | [Robust linear model](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_rlm.html) |
| IonStar / MaxQuant | [Two-factor model](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_Model_IonStar_With2Factors.html) |
| IonStar / FragPipe | [combined_protein.tsv](https://prolfqua.github.io/prolfquabenchmark/articles/BenchmarkFragPipeProteinIonStar.html) |
| IonStar / FragPipe | [MSstats.tsv](https://prolfqua.github.io/prolfquabenchmark/articles/BenchmarkFragPipeMSStats.html) |
| CPTAC / MaxQuant | [peptide.txt](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_cptac.html) |

The `inst/MyArticle/` directory contains the source for the prolfqua manuscript (`paper.Rmd`) and its supplement (`prolfqua_supplement.Rmd`), buildable via `make` from that directory.

## Installation

```r
remotes::install_github("wolski/prolfquabenchmark", dependencies = TRUE)
```

## Building vignettes

```bash
make vignette V=BenchmarkingIonstarData   # single vignette
make check                                 # all vignettes + R CMD check
```
