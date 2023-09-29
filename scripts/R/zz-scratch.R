library(arrow)
library(dplyr)

metadata <- read_feather("data/ecosis-processed/wisc-leaf-trait-vine/metadata.arrow")
spectra <- arrow_table("data/ecosis-processed/wisc-leaf-trait-vine/spectra.arrow")

################################################################################
library(tidyverse)

params <- read_csv("data/params.csv")

params %>%
  filter(version == "4", parameters %in% c("N", "Ccab", "Cw", "Cm")) %>%
  select(parameters, mean) %>%
  mutate(rn = (row_number()-1) %/% 4) %>%
  pivot_wider(names_from = "parameters", values_from = "mean") %>%
  select(-rn) %>%
  as.matrix() %>%
  pairs()

ggplot(params) +
  aes(x = mean, fill = version) +
  geom_density() +
  facet_wrap(~parameters, scales = "free")
