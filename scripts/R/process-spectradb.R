stopifnot(
  requireNamespace("yaml"),
  requireNamespace("tidyr")
)

source("scripts/R/common.R")

library(dplyr)
library(readr)
library(arrow)
library(purrr)
library(fs)

dataset_ids <- path_file(dir_ls("~/projects/spectra_db", type = "d"))
overwrite <- FALSE

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

for (dataset_id in dataset_ids) {
  print(paste0("Dataset ID: ", dataset_id))
  outdir <- file.path("data", "ecosis-processed", dataset_id)
  metadata_file <- file.path(outdir, "metadata.arrow")
  spectra_file <- file.path(outdir, "spectra.arrow")
  if (!overwrite && file.exists(metadata_file)) {
    print("Skipping")
    next
  }
  data_root <- paste0("~/projects/spectra_db/", dataset_id, "/")

  # TODO: Fill in other metadata?
  header <- list(dataset_id = dataset_id)

  metadata <- read_csv(file.path(data_root, "metadata.csvy"), comment = "#", show_col_types = FALSE) %>%
    mutate(entity_id = observation_id)

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
      mutate(observation_id = gsub("(.*)\\.[[:digit:]]+$", "\\1", spectra_id)) %>%
      left_join(spec_types, by = "spectra_id")
    dat_long
  }

  spectra_files <- list.files(file.path(data_root, "spectra"), full.names = TRUE)
  spectra <- map(spectra_files, parse_spectra) %>%
    list_rbind()

  stopifnot(all(unique(spectra[["observation_id"]]) %in% metadata[["observation_id"]]))

  spectra_arrow <- arrow_table(spectra)
  metadata_arrow <- arrow_table(metadata)
  metadata_arrow$metadata <- header

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write_feather(metadata_arrow, metadata_file)
  write_feather(spectra_arrow, spectra_file)
}
