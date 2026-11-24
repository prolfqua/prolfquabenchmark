# prolfquabenchmark

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18775961.svg)](https://doi.org/10.5281/zenodo.18775961)

Benchmark of differential expression methods for label-free proteomics — the **prolfqua** models and **MSstats**, **proDA** and **msqrob2** — on the IonStar (E. coli in human) and CPTAC (UPS1 in yeast) spike-in datasets, with the scores of every prolfqua release since 1.9.0.

Rendered results: **https://prolfqua.github.io/prolfquabenchmark/**

- [Benchmark results](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_results.html): current scores, ROC and FDP-vs-FDR curves per dataset and method
- [Benchmark scores across prolfqua releases](https://prolfqua.github.io/prolfquabenchmark/articles/Score_history.html): AUC, pAUC and step timings per release

## The benchmark pipeline

A Snakemake pipeline in this repository runs the benchmark; the vignettes only show its results.

- `config/config.yaml`: the datasets (input file, the inputs they provide, preprocessing, contrasts, ground truth) and the methods (tool, the inputs they need, level); every method runs on every dataset that provides its inputs
- `workflow/Snakefile`, `workflow/scripts/`: prepare each dataset once, fit every method, score, export to `inst/extdata`
- Uses the installed R packages (prolfqua, prolfquapp, prolfquadata and the Suggests of `DESCRIPTION`); Snakemake comes from the uv environment of `pyproject.toml`

```bash
make install            # the pipeline loads the installed prolfquabenchmark
make benchmark J=4      # run the pipeline, export inst/extdata/benchmark_*.tsv
make benchmark-dry      # list the jobs
make score-history      # append the scores and step timings of the installed prolfqua to the history
make score-history-release V=1.9.0 PROLFQUA=f4dfb7f3 PROLFQUAPP=18914707  # the same for an earlier release
```

A method is added in `config/config.yaml` (a prolfqua facade needs no new code); a dataset needs a `prepare_<dataset>.R` script in `workflow/scripts`.

The `inst/MyArticle/` directory contains the source for the prolfqua manuscript (`paper.Rmd`) and its supplement (`prolfqua_supplement.Rmd`), buildable via `make` from that directory.

## Installation

```r
remotes::install_github("prolfqua/prolfquabenchmark", dependencies = TRUE)
```
