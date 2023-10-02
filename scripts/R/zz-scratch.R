library(arrow)
library(dplyr)

metadata <- read_feather("data/ecosis-processed/wisc-leaf-trait-vine/metadata.arrow")
spectra <- arrow_table("data/ecosis-processed/wisc-leaf-trait-vine/spectra.arrow")

################################################################################
library(tidyverse)
library(fs)
library(arrow)

params <- read_csv("data/results/params.csv")

read_metadata <- function(fname) {
  # fname <- metadata_files[1]
  draw <- read_feather(fname, as_data_frame = FALSE)
  dataset_id <- draw$metadata[["dataset_id"]]
  dat <- as_tibble(draw)
  dat[["dataset_id"]] <- dataset_id
  # Force certain cols to be character
  for (ccol in c("collection_date", "variety")) {
    if (ccol %in% colnames(dat)) {
      dat[[ccol]] <- as.character(dat[[ccol]])
    }
  }
  dat
}

metadata_files <- dir_ls("data/ecosis-processed/", glob = "**/metadata.arrow", recurse = TRUE)
metadata <- metadata_files %>%
  map(read_metadata) %>%
  list_rbind()

p4wide <- params %>%
  filter(version == "4", parameters %in% c("N" ,"Ccab", "Cw", "Cm")) %>%
  select(observation_id, dataset_id, parameters, mean) %>%
  pivot_wider(names_from = "parameters", values_from = "mean")

p4results <- p4wide %>%
  inner_join(metadata)

metadata %>%
  filter(dataset_id == "ecosis_cedarcreek_biodiversity") %>%
  select(where(~any(!is.na(.)))) %>%
  ggplot() +
  aes(x = dataset_id, y = leaf_chltot_per_area) +
  geom_violin() +
  geom_jitter()

p4results %>%
  filter(!is.na(leaf_chltot_per_area)) %>%
  ggplot() +
  aes(x = leaf_chltot_per_area, y = Ccab) +
  geom_point() +
  geom_abline(color = "red") +
  facet_wrap(vars(dataset_id)) +
  theme_bw() +
  labs(x = "Observed", y = "Predicted")

# Quick summary pairs plot
params %>%
  filter(version == "4", parameters %in% c("N", "Ccab", "Cw", "Cm")) %>%
  select(observation_id, parameters, mean) %>%
  pivot_wider(names_from = "parameters", values_from = "mean") %>%
  select(-observation_id) %>%
  as.matrix() %>%
  pairs(pch = 19, alpha = 0.2)

# Join with metadata
metadata <- 

ggplot(params) +
  aes(x = mean, fill = version) +
  geom_density() +
  facet_wrap(~parameters, scales = "free")
