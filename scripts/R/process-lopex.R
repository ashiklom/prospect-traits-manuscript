stopifnot(
  requireNamespace("yaml"),
  requireNamespace("tidyr")
)

library(dplyr)
library(readr)
library(arrow)
library(purrr)

data_root <- "~/projects/spectra_db/lopex/"

dataset_id <- "lopex"
header <- list(
  dataset_id = dataset_id,
  dataset_description = paste(
    "Leaf optical properties experiment (LOPEX 93)"
  ),
  ecosis_id = "leaf-optical-properties-experiment-database--lopex93-"
)

metadata <- read_csv(file.path(data_root, "metadata.csvy"), comment = "#") %>%
  mutate(entity_id = observation_id)

read_csvy_header <- function(filename) {
  fo <- file(filename, "r")
  on.exit(close(fo), add = TRUE)
  header_lines <- character(0)
  repeat {
    nextline <- readLines(fo, 1)
    if (grepl("^ *#", nextline)) {
      header_lines <- c(header_lines, nextline)
    } else {
      break
    }
  }
  header <- paste(gsub("^#", "", header_lines), collapse = "\n")
  yaml::yaml.load(header)
}

parse_spectra <- function(filename) {
  meta <- read_csvy_header(filename)
  columns <- map_chr(meta$resources$fields, "name") %>%
    grep(pattern = "wavelengths", invert = TRUE, value = TRUE)
  spec_types <- tibble(
    spectra_id = columns,
    spectral_measurement = c(
      "R" = "reflectance",
      "T" = "transmittance"
    )[meta$spectra_types]
  )
  dat <- read_csv(filename, comment = "#", show_col_types = FALSE)
  dat_long <- dat %>%
    tidyr::pivot_longer(
      !wavelengths,
      names_to = "spectra_id",
      values_to = "value"
    ) %>%
    rename(wavelength_nm = wavelengths) %>%
    mutate(observation_id = gsub("(lopex_[[:digit:]]+)\\..*", "\\1", spectra_id)) %>%
    left_join(spec_types, by = "spectra_id")
  dat_long
}

spectra_files <- list.files(file.path(data_root, "spectra"), full.names = TRUE)
spectra <- map(spectra_files, parse_spectra) %>%
  list_rbind()

spectra_arrow <- arrow_table(spectra)
metadata_arrow <- arrow_table(metadata)
metadata_arrow$metadata <- header

outdir <- file.path("data", "ecosis-processed", dataset_id)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write_feather(metadata_arrow, file.path(outdir, "metadata.arrow"))
write_feather(spectra_arrow, file.path(outdir, "spectra.arrow"))
