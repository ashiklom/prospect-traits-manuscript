stopifnot(
  requireNamespace("cmdstanr"),
  requireNamespace("posterior"),
  requireNamespace("rrtm")
)

resample_spec_avg <- function(x, waves) {
  stopifnot(
    min(waves) >= 400,
    max(waves) <= 2500
  )
  x <- as.matrix(x)
  waves_ref <- seq(400, 2500, 1)
  wstart <- head(waves, -1)
  wend <- tail(waves, -1)
  xout <- numeric(waves)
  for (i in seq_along(waves)) {
    ii <- which(waves_ref >= wstart, waves_ref <= wend)
    xout[i] <- colMeans(x[ii, ])
  }
  drop(xout)
}

resample_spec_interp <- function(x, waves) {
  waves_ref <- seq(400, 2500, 1)
  x <- as.matrix(x)
  xout <- matrix(NA_real_, length(waves), NCOL(x))
  for (i in seq_len(NCOL(x))) {
    xout[, i] <- approx(waves_ref, x[, i], waves)[["y"]]
  }
  drop(xout)
}

make_data_list <- function(waves, reflectance,
                           prospect_version = "pro",
                           resample_func = resample_spec_interp) {
  stopifnot(
    length(waves) == NROW(reflectance),
    prospect_version %in% c("4", "5", "5b", "d", "pro")
  )
  if (prospect_version %in% c("4", "5", "5b")) {
    refr <- rrtm:::refractive_p45
    if (prospect_version == "4") {
      kmat <- rrtm:::dataspec_p4
    } else if (prospect_version == "5") {
      kmat <- rrtm:::dataspec_p5[, c("Cab", "Car", "Cw", "Cm")]
    } else if (prospect_version == "5b") {
      kmat <- rrtm::dataspec_p5
    }
    talf_angle <- 40
  } else if (prospect_version == "d") {
    refr <- rrtm:::refractive_pd
    kmat <- rrtm:::dataspec_pd
    talf_angle <- 40
  } else if (prospect_version == "pro") {
    refr <- rrtm:::refractive_ppro
    kmat <- rrtm:::dataspec_ppro
    talf_angle <- 40
  }

  # Pre-calculate tav -- this is a relatively expensive operation
  talf <- rrtm:::tav_abs(talf_angle, refr)
  t12 <- rrtm:::tav_abs(90, refr)
  t21 <- t12 / (refr ^ 2)

  spec_quantities <- list(
    talf = talf,
    t12 = t12,
    t21 = t21,
    kmat = kmat
  )
  obs <- as.matrix(reflectance)
  c(
    list(
      nwl = length(waves),
      nobs = ncol(obs),
      obs = obs
    ),
    lapply(spec_quantities, resample_func, waves = waves)
  )
}

fit_prospect <- function(waves,
                         reflectance,
                         prospect_version, ...) {
  data_list <- make_data_list(waves, reflectance, prospect_version = prospect_version)
  modstring <- prospect_stan(prospect_version)
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(modstring))
  mod$sample(data_list, ...)
}
