.PHONY: all datasets

all: _output/prospect-traits-manuscript.html

_output/prospect-traits-manuscript.html: prospect-traits-manuscript.qmd
	quarto render

datasets: \
	data/ecosis-processed/wisc-leaf-trait-vine/spectra.arrow \
	data/ecosis-processed/lopex/spectra.arrow

data/ecosis-processed/%/spectra.arrow: \
	scripts/R/process-%.R
	Rscript $<
