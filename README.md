
This github project has an accompanying website created with pkgdown.
https://prolfqua.github.io/prolfquabenchmark/

# prolfquaBenchmark - Benchmarking LFQ data analysis using the prolfqua package.

The R package prolfquabenchmark contains vignettes which show how we used the IonStar dataset to asses the performance of the methods implemented in the [prolfqua](https://github.com/fgcz/prolfqua) R packages and other packages (MSstats, proDA).

- [Benchmarking FragPipe output using Ionstar Dataset](https://prolfqua.github.io/prolfquabenchmark/articles/BenchmarkFragPipeProteinIonStar.html)
- [Benchmarking MSstats using the Ionstar Dataset](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_MSStats.html)
- [Benchmarking the proDA package using the Ionstar Dataset MQ LFQ intensities](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_proDA_fromMQlfq.html)
- [Benchmarking the proDA package using the Ionstar Dataset starting from peptides](https://prolfqua.github.io/prolfquabenchmark/articles/Benchmark_proDA_medpolish.html)
- [Benchmarking normalization, aggregation and models using the Ionstar Dataset](https://prolfqua.github.io/prolfquabenchmark/articles/BenchmarkingIonstarData.html)



## How to install and generate the vignettes

Download the latest _prolfqua_ release from from https://github.com/fgcz/prolfqua/releases
And then install it by running on the command:

```
R CMD INSTALL prolfqua_X.Y.Z.tar.gz
```

Or in the R session:

```r
install.packages("prolfqua_X.Y.Z.tar.gz",repos = NULL, type="source")
```

Afterwards you also will need to install the `prolfquadata` package

```r
install.packages('remotes')
remotes::install_gitlab("wolski/prolfquadata", host="gitlab.bfabric.org")

```

Finally you can clone the repository and open it in RStudio and build the vignettes with

```r
devtools::build_vignettes()
```

Or build them on the command line by running

```r
R CMD build prolfquaBenchmark
```
