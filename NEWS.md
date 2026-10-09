# prolfquabenchmark 0.3.4

- The example datasets `data_benchmarkExample` and `data_test_confusion_matrix_scores` now ship with prolfquabenchmark (they were removed from prolfqua 1.9.0); load them with `prolfqua::prolfqua_data(<name>, package = "prolfquabenchmark")`.
- `SqueezeVar_comparison` no longer lists `lm_missing` among the prolfqua facades that feed variance moderation.
- The vignettes `Benchmark_prolfqua`, `Benchmark_pipeline_demo` and `Benchmark_Model_IonStar_With2Factors` take `ContrastsMissing` and `ContrastsLMMissingFacade` from the new package prolfquatakeoutmissing (Suggests), since prolfqua 1.9.0 no longer has them.

- New vignette `SqueezeVar_comparison` comparing variance moderation across `limma::squeezeVar()`, `msqrob2`, and
  `prolfqua::squeezeVarRob()`. It shows msqrob2 delegates to limma, and that `prolfqua::squeezeVarRob()` is limma's
  classic `fitFDist` plus a `min_df` filter — explaining why old MSqRob forked it (classic limma collapses to
  `df.prior = Inf` on fractional/low-df features). It documents that `prolfqua` reaches this regime through its `rlm` /
  `rfit`, `lmer_nested`, and imputation facades (fractional / low residual df; `firth` is not a consumer), and includes
  a benchmark (synthetic + real IonStar spike-in data) showing the former fork and modern `limma::squeezeVar()` give
  near-identical results — the evidence behind prolfqua migrating its moderation to `limma::squeezeVar()`.
- Began tracking user-visible changes in `NEWS.md`. For changes before this version, see the git history.
