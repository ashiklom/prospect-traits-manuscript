invisible(lapply(list.files("./R", pattern = ".R$", full.names = TRUE), source))
dat <- read.csv("~/projects/spectra_db/lopex/spectra/lopex_01.csvy", comment.char = "#")
waves <- dat[, "wavelengths"]
refl <- as.matrix(dat[, c(2, 4, 6)])
test4 <- fit_prospect(waves, refl, "4", chains = 1)
test4
# test4_2 <- fit_prospect(waves, refl, "4", chains = 1)
# test5 <- fit_prospect(waves, refl, "5", chains = 1)
# test5b <- fit_prospect(waves, refl, "5b", chains = 1)
# testd <- fit_prospect(waves, refl, "d", chains = 1)

# testpro <- fit_prospect(waves, refl, "pro", chains = 1)

library(ggplot2)
library(dply)
test4_df <- posterior::as_draws_df(test4)

ggplot(test4_df) +
  aes(x = Ccab) +
  geom_histogram()
