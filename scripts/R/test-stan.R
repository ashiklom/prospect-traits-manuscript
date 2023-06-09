invisible(lapply(list.files("./R", pattern = ".R$", full.names = TRUE), source))
dat <- read.csv("~/projects/spectra_db/lopex/spectra/lopex_01.csvy", comment.char = "#")
waves <- dat[, "wavelengths"]
refl <- as.matrix(dat[, c(2, 4, 6)])
system.time(
  test4 <- fit_prospect(waves, refl, "4", chains = 1)
)
test4
# test4_2 <- fit_prospect(waves, refl, "4", chains = 1)
# test5 <- fit_prospect(waves, refl, "5", chains = 1)
# test5b <- fit_prospect(waves, refl, "5b", chains = 1)
# testd <- fit_prospect(waves, refl, "d", chains = 1)

# testpro <- fit_prospect(waves, refl, "pro", chains = 1)

library(ggplot2)
library(posterior)
library(tidyverse)
test4_df <- as_draws_df(test4)

test4_long <- test4_df %>%
  select(-lp__) %>%
  pivot_longer(
    (N:rsd),
    names_to = "variable",
    values_to = "value"
  )

ggplot(test4_long) +
  aes(x = value) +
  geom_histogram() +
  facet_wrap(vars(variable), scales = "free_x")

ggplot(test4_df) +
  aes(x = Ccab) +
  geom_histogram()
