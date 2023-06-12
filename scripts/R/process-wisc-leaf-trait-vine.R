library(dplyr)
library(tidyr)
library(arrow)

dataset_id <- "wisc-leaf-trait-vine"

header <- list(
  dataset_id = dataset_id,
  dataset_description = paste0(
    "Leaf level spectra and LMA for a set of trees, forbs, vines,",
    "and grasses collected in Madison, WI"
  ),
  ecosis_id = "leaf-level-spectra-and-lma-for-a-set-of-trees--forbs--vines-and-grasses-collected-in-madison--wi",
  leaf_mass_per_area = "g cm-2",
  functional_type = "Plant functional type. One of: tree, forb, grass, vine."
)

data_path <- "data/ecosis-download/leaf-level-spectra-and-lma-for-a-set-of-trees--forbs--vines-and-grasses-collected-in-madison--wi.csv"
message("Reading dataset...")
dat <- read_csv_arrow(data_path)

message("Updating metadata...")
dat2 <- dat %>%
  arrange(id) %>%
  mutate(
    entity_id = paste(
      dataset_id,
      sprintf("%05d", id),
      sep = "|"
    ),
    id = NULL,
    observation_id = paste(
      entity_id,
      sprintf("%05d", row_number()),
      sep = "|"
    ),
    spectra_id = paste(observation_id, "R", sep = "|"),
    leaf_mass_per_area = lma_g_m2,
    lma_g_m2 = NULL,
    spectral_measurement = "reflectance",
    .before = 1
  )

message("Pivoting spectra...")
spectra <- dat2 %>%
  select(observation_id, spectra_id, `350`:`2500`) %>%
  pivot_longer(
    !c(spectra_id, observation_id),
    names_to = "wavelength_nm",
    names_transform = as.numeric,
    values_to = "value"
  ) %>%
  arrow_table()

message("Post-processing metadata...")
metadata <- dat2 %>%
  select(-(`350`:`2500`), -spectra_id) %>%
  arrow_table()
metadata$metadata <- header

message("Writing...")
outdir <- file.path("data", "ecosis-processed", dataset_id)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
message("...metadata...")
write_feather(metadata, file.path(outdir, "metadata.arrow"))
message("...spectra...")
write_feather(spectra, file.path(outdir, "spectra.arrow"))

message("Done!")
