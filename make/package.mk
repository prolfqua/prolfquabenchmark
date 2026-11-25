# Build the documentation site with altdoc instead of pkgdown.
SITE_CMD := Rscript -e "altdoc::render_docs(freeze = FALSE)"

INSTALL_BUILD_TARGET = build-install

.PHONY: build-install benchmark benchmark-dry score-history score-history-release

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

# An earlier prolfqua release: installs it with the prolfquapp of its time into results/releases/$(V)/lib and appends
# its scores and step timings, dated with its prolfqua commit. PROLFQUA defaults to the tag $(V).
#   make score-history-release V=1.9.0 PROLFQUA=f4dfb7f3 PROLFQUAPP=18914707 [REPLACE=1]
PROLFQUA ?= $(V)
score-history-release:
	inst/scripts/install_release.sh $(V) $(PROLFQUA) $(PROLFQUAPP)
	$(SNAKEMAKE) -j $(J) score_history $(if $(REPLACE),--forcerun score_history) --config release=$(V) \
	  release_date=$$(git -C ../prolfqua log -1 --format=%ad --date=short $(PROLFQUA)) $(if $(REPLACE),replace=true)

help-package:
	@echo ""
	@echo "Package-specific:"
	@echo "  make install builds without vignettes"
	@echo "  make benchmark [J=4]         - run the benchmark pipeline and export inst/extdata/benchmark_*.tsv"
	@echo "  make benchmark-dry           - list the jobs the pipeline would run"
	@echo "  make score-history           - append scores and step timings of the installed prolfqua to"
	@echo "                                 inst/extdata/score_history.tsv"
	@echo "  make score-history REPLACE=1 - same, replacing the rows of a version already in the history"
	@echo "  make score-history-release V=<version> [PROLFQUA=<commit>] PROLFQUAPP=<commit> [REPLACE=1]"
	@echo "                               - same for an earlier release, installed into results/releases/<version>"
