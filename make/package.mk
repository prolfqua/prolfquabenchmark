INSTALL_BUILD_TARGET = build-install

.PHONY: build-install

build-install: document
	Rscript -e "devtools::build(vignettes = FALSE)"

help-package:
	@echo ""
	@echo "Package-specific:"
	@echo "  make install builds without the long-running benchmark vignettes"
