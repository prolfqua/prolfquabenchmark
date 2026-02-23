.PHONY: all check check-fast test build build-vignettes document coverage install lint format clean help site deploy renv-init renv-restore renv-snapshot renv-reset

all: check

help:
	@echo "prolfquabenchmark development targets:"
	@echo ""
	@echo "  make all             - full pipeline: document -> build -> check (default)"
	@echo "  make check           - R CMD check (runs document, build first)"
	@echo "  make check-fast      - R CMD check without rebuilding vignettes"
	@echo "  make build           - build tarball (runs document first)"
	@echo "  make build-vignettes - build all vignettes into inst/doc"
	@echo "  make vignette V=Name - build a single vignette (without .Rmd extension)"
	@echo "  make document        - generate roxygen2 docs"
	@echo "  make test            - run testthat tests (runs document first)"
	@echo "  make coverage        - code coverage report"
	@echo "  make install         - install package locally"
	@echo "  make lint            - run lintr"
	@echo "  make format          - format with air"
	@echo "  make site            - build pkgdown site locally"
	@echo "  make deploy          - build pkgdown site and push to gh-pages"
	@echo "  make clean           - remove build artifacts"
	@echo ""
	@echo "  Environment (renv):"
	@echo "  make renv-init       - initialize renv and install all deps (first time)"
	@echo "  make renv-restore    - restore environment from renv.lock"
	@echo "  make renv-snapshot   - update renv.lock after installing new packages"
	@echo "  make renv-local      - reinstall prolfqua + prolfquapp from local sibling dirs"
	@echo "  make renv-reset      - nuke renv library + lockfile, reinit from DESCRIPTION"

document:
	Rscript -e "devtools::document()"

build: document
	Rscript -e "devtools::build()"

check: build
	Rscript -e "devtools::check()"

check-fast: document
	Rscript -e "devtools::check(build_args = '--no-build-vignettes', args = '--no-vignettes')"

build-vignettes: document
	Rscript -e "devtools::build_vignettes()"
	mkdir -p inst/doc
	cp doc/*.html doc/*.Rmd doc/*.R inst/doc/ 2>/dev/null || true

# Build a single vignette by name: make vignette V=BenchmarkingIonstarData
vignette:
ifndef V
	$(error Usage: make vignette V=<vignette_name>, e.g. make vignette V=BenchmarkingIonstarData)
endif
	Rscript -e "rmarkdown::render('vignettes/$(V).Rmd')"

test: document
	Rscript -e "devtools::test()"

coverage: document
	Rscript -e "covr::package_coverage() |> print()"

install: document
	Rscript -e "devtools::install()"

lint:
	Rscript -e "lintr::lint_package()"

format:
	air format .

site: document
	Rscript -e "pkgdown::build_site()"

deploy: document
	Rscript -e "pkgdown::deploy_to_branch()"

renv-init:
	Rscript -e "renv::init(bioconductor = TRUE)"

renv-restore:
	Rscript -e "renv::restore()"

renv-snapshot:
	Rscript -e "renv::snapshot()"

renv-local:
	Rscript -e "devtools::install('../prolfqua', upgrade = 'never')"
	Rscript -e "devtools::install('../prolfquapp', upgrade = 'never')"

renv-reset:
	rm -rf renv/library renv.lock
	Rscript -e "renv::init(bioconductor = TRUE)"

clean:
	rm -rf *.Rcheck
	rm -f Rplots.pdf
	rm -rf inst/doc doc Meta
	rm -f vignettes/*.html vignettes/*.R
