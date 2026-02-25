.PHONY: all paper supplement clean help

# Project root (where renv .Rprofile lives)
PKG_ROOT := $(realpath ../..)
ARTICLE_DIR := $(realpath .)

# Auto-detect latest dated benchmark directory, or override with: make paper BENCH_DATE=20260221
BENCH_DATE ?= $(shell ls -1d $(PKG_ROOT)/inst/Benchresults/[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9] 2>/dev/null | sort -r | head -1 | xargs basename)

export BENCH_DATE

all: paper supplement

help:
	@echo "MyArticle build targets:"
	@echo ""
	@echo "  make all              - build paper + supplement PDFs (default)"
	@echo "  make paper            - render paper.Rmd to PDF"
	@echo "  make supplement       - render prolfqua_supplement.Rmd to PDF"
	@echo "  make clean            - remove generated outputs"
	@echo ""
	@echo "  Override benchmark date:  make paper BENCH_DATE=20260221"
	@echo "  Current BENCH_DATE: $(BENCH_DATE)"

paper:
	@echo "Building paper.Rmd with BENCH_DATE=$(BENCH_DATE)"
	cd $(PKG_ROOT) && Rscript -e "rmarkdown::render('$(ARTICLE_DIR)/paper.Rmd')"

supplement:
	@echo "Building prolfqua_supplement.Rmd with BENCH_DATE=$(BENCH_DATE)"
	cd $(PKG_ROOT) && Rscript -e "rmarkdown::render('$(ARTICLE_DIR)/prolfqua_supplement.Rmd')"

clean:
	rm -f paper.pdf paper.tex paper.log paper.html
	rm -f prolfqua_supplement.pdf prolfqua_supplement.tex prolfqua_supplement.log prolfqua_supplement.html
	rm -rf paper_files prolfqua_supplement_files
	rm -rf _bookdown_files
