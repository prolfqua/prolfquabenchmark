INSTALL_BUILD_TARGET = build-install

.PHONY: build-install benchmark benchmark-dry score-history

build-install: document
	Rscript -e "devtools::build(vignettes = FALSE)"

# The benchmark pipeline (workflow/Snakefile, config/config.yaml) uses the installed packages: run `make install`
# first. Python dependencies (snakemake) come from the uv environment of pyproject.toml.
SNAKEMAKE = uv run snakemake -s workflow/Snakefile --configfile config/config.yaml
J ?= 4

benchmark:
	$(SNAKEMAKE) -j $(J)

benchmark-dry:
	$(SNAKEMAKE) -n

score-history:
	$(SNAKEMAKE) -j $(J) score_history $(if $(REPLACE),--forcerun score_history --config replace=true)

help-package:
	@echo ""
	@echo "Package-specific:"
	@echo "  make install builds without vignettes"
	@echo "  make benchmark [J=4]         - run the benchmark pipeline and export inst/extdata/benchmark_*.tsv"
	@echo "  make benchmark-dry           - list the jobs the pipeline would run"
	@echo "  make score-history           - append scores and step timings of the installed prolfqua to"
	@echo "                                 inst/extdata/score_history.tsv"
	@echo "  make score-history REPLACE=1 - same, replacing the rows of a version already in the history"
