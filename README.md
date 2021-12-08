# proLFQuaBenchmark - Benchmarking LFQ data analysis using the prolfqua package.

This package contains vignettes only. The vignettes show how we used benchmark datasets to asses the performance of the methods implemented in prolfqua and other packages.
The results best be studied here by browsing the vignettes online : https://wolski.github.io/proLFQuaBenchmark/.


## How to install and generate the vignettes

Download the latest prolfqua release from from https://github.com/wolski/prolfqua/releases

And then install it by running on the command:

```
R CMD INSTALL prolfqua_X.Y.Z.tar.gz
```

Or in the R session:
```
install.packages("prolfqua_X.Y.Z.tar.gz",repos = NULL, type="source")
```

Afterwards you also will need to install the `prolfquaData` package

```
install.packages('remotes')
remotes::install_gitlab("wolski/prolfquaData", host="gitlab.bfabric.org")

```

Finally you can clone the repository and open it in RStudio and build the vignettes with

```{r}
devtools::build_vignettes()
```

Or build them on the command line by running

```{r}
R CMD build prolfquaBenchmark
```
