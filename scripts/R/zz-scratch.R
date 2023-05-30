library(arrow)
library(dplyr)

metadata <- read_feather("data/ecosis-processed/wisc-leaf-trait-vine/metadata.arrow")
spectra <- arrow_table("data/ecosis-processed/wisc-leaf-trait-vine/spectra.arrow")
