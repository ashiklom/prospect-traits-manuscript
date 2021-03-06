.PHONY: all build deploy pkg

all: build deploy

pkg:
	Rscript -e 'devtools::document(".")'
	Rscript -e 'devtools::install(".")'

build: extdata/results.db data/siteinfo.RData
	./_build.sh

deploy: build
	./_deploy.sh

## Results SQLite table
extdata/results.db:
	mkdir -p extdata
	rsync -avz --progress geo:~/dietzelab/prospectinversion/scripts/results.db $@

## Site climate data
data/site_info.RData: extdata/temperature_matrix.rds extdata/precipitation_matrix.rds
	Rscript scripts/make_site_info.R

extdata/temperature_matrix.rds extdata/precipitation_matrix.rds:
	Rscript scripts/get_climate_matrices.R

## Results data ready for analysis
data/results_all.RData: scripts/setup.R extdata/results.db
	Rscript scripts/setup.R
