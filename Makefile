.PHONY: all check check-fast test build build-vignettes document coverage install lint format clean help site deploy

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

# Build a single vignette by name: make vignette V=Benchmark_prolfqua
vignette:
ifndef V
	$(error Usage: make vignette V=<vignette_name>, e.g. make vignette V=Benchmark_prolfqua)
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

clean:
	rm -rf *.Rcheck
	rm -f Rplots.pdf
	rm -rf inst/doc doc Meta
	rm -f vignettes/*.html vignettes/*.R
