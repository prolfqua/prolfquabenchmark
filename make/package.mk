INSTALL_BUILD_TARGET = build-install

.PHONY: build-install score-history

build-install: document
	Rscript -e "devtools::build(vignettes = FALSE)"

SCORE_VIGNETTES = Benchmark_prolfqua Benchmark_MSStats Benchmark_proDA Benchmark_msqrob2 Benchmark_hurdle Benchmark_ropeca Benchmark_null_split

score-history:
	touch inst/Benchresults/.score_run
	for v in $(SCORE_VIGNETTES); do $(MAKE) vignette V=$$v || exit 1; done
	Rscript inst/scripts/timing_input.R inst/Benchresults/timing_input.rds
	Rscript inst/scripts/time_prolfqua_steps.R inst/Benchresults/timing_input.rds inst/Benchresults/timing
	Rscript inst/scripts/update_score_history.R $(if $(REPLACE),--replace)

help-package:
	@echo ""
	@echo "Package-specific:"
	@echo "  make install builds without the long-running benchmark vignettes"
	@echo "  make score-history          - render the benchmark vignettes (~18 min), time the core prolfqua steps"
	@echo "                                and append scores and timings for the installed prolfqua to"
	@echo "                                inst/extdata/score_history.tsv"
	@echo "  make score-history REPLACE=1 - same, replacing the rows of a version already in the history"
