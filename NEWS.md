# prolfquabenchmark 0.3.4

- New vignette `Benchmark_null_split`: splits the runs of each IonStar dilution at random into two pseudo-groups and tests them against each other with `lm`, `limma` and `lm_impute` after no normalization, `robscale()` or human-subset `robscale_subset()`. It reports the share of p < 0.05 and p < 0.01, FDR hits, KS distance and p-value histograms over 20 splits, plus a split along the two run blocks. Both normalizations are calibrated; without normalization single splits are far off; the run-block split gives about 100 FDR hits per 1000 proteins.

- New vignette `Benchmark_ropeca` comparing prolfqua's `ropeca_nested` with the new `ropecaV2_nested`, `lmer_nested` and `limma`: on simulated proteins with correlated peptides and on the unchanged human proteins of IonStar it tabulates the share of p < 0.05 and p < 0.01 by number of peptides, and on IonStar it scores the E. coli vs human ranking (AUC, partial AUC, FDR vs FDP). Requires prolfqua >= 1.9.4.
- `Benchmark_hurdle` drops `firth_nested`, which prolfqua 1.9.3 removed, and its hurdle variants; the msqrob2-rule variant now combines `lm_impute` with `binomial_firth_nested` (`hurdle_binomial_firth_df4`). The calibration text keeps the `firth_nested` result: the share of unchanged proteins with detection p < 0.05 rose from 1.3% to 7.7% with the number of peptides.
- New vignette `Benchmark_hurdle` benchmarking prolfqua's hurdle contrasts (`ContrastsHurdleFacade`: `lm_impute` combined with `firth_nested`, `binomial_nested` or `binomial_firth_nested`) against `lm_impute` alone and msqrob2's hurdle model on the IonStar MaxQuant peptides. It also scores msqrob2's rule of combining the detection p-value for fully detected proteins too, and tabulates how often each detection model calls the unchanged human proteins significant, by number of peptides.
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
