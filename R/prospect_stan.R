stopifnot(
  requireNamespace("glue")
)
library(distributional)

prospect_priors <- c(
  N = dist_truncated(dist_normal(1.4, 0.2), lower = 1),
  Ccab = dist_truncated(dist_normal(40, 20), lower = 0),
  Ccar = dist_truncated(dist_normal(40, 25), lower = 0),
  Canth = dist_truncated(dist_normal(40, 25), lower = 0),
  Cbrown = dist_exponential(2.0),
  Cw = dist_truncated(dist_normal(0.01, 0.01), lower = 0),
  Cm = dist_truncated(dist_normal(0.01, 0.01), lower = 0),
  Cprot = dist_truncated(dist_normal(0.01, 0.01), lower = 0),
  Ccbc = dist_truncated(dist_normal(0.01, 0.01), lower = 0),
  rsd = dist_truncated(dist_cauchy(0, 1.0), lower = 0)
)

dist2stan <- function(x) {
  xraw <- unclass(x)[[1]]
  if ("dist" %in% names(xraw)) {
    # Truncated distribution
    xdist <- xraw[["dist"]]
  } else {
    xdist <- xraw
  }
  params <- unlist(xdist)
  paramstring <- paste(params, collapse = ", ")
  dist_type <- attr(xdist, "class")[1]
  stanstring <- sprintf(
    "%s ~ %s(%s)",
    names(x),
    gsub("dist_", "", dist_type),
    paramstring
  )
  stanstring
}

prospect_stan <- function(prospect_version) {
  expint_func <- "
    real expint(real x) {
      real A;
      real B;
      A = log((0.56146 / x + 0.65) * (1 + x));
      B = x^4 * exp(7.7 * x) * (2 + x)^3.7;
      return (A^-7.7 + B)^-0.13;
    }
  "
  transfun_func <- "
    real transfun(real k) {
      real trans;
      if (k <= 0) {
        trans = 1.0;
      } else {
        trans = (1 - k) * exp(-k) + k^2 * expint(k);
      }
      if (trans < 0) {
        trans = 0.0;
      }
      return trans;
    }
  "

  gpm_func <- "
    real gpm(real N, real k,
             data real talf, data real ralf,
             data real t12, data real r12,
             data real t21, data real r21) {
      real trans = transfun(k);
      real denom = 1 - r21 ^ 2 * trans ^ 2;
      real Ta = talf * trans * t21 / denom;
      real Ra = ralf + r21 * trans * Ta;
      real t = t12 * trans * t21 / denom;
      real r = r12 + r21 * trans * t;
      real Tsub;
      real Rsub;
      if (r + t >= 1) {
        Tsub = t / (t + (1 - t) * (N - 1));
        Rsub = 1 - Tsub;
      } else {
        real D = sqrt((1 + r + t) * (1 + r - t) * (1 - r + t) * (1 - r - t));
        real r2 = r ^ 2;
        real t2 = t ^ 2;
        real va = (1 + r2 - t2 + D) / (2 * r);
        real vb = (1 - r2 + t2 + D) / (2 * t);
        real vbNN = vb ^ (N - 1);
        real vbNN2 = vbNN ^ 2;
        real va2 = va ^ 2;
        real denomx = va2 * vbNN2 - 1;
        Rsub = va * (vbNN2 - 1) / denomx;
        Tsub = vbNN * (va2 - 1) / denomx;
      }
      real denomy = 1 - Rsub * r;
      // Reflectance
      real RN = Ra + Ta * Rsub * t / denomy;
      // Transmittance
      // TN <- Ta * Tsub / denomy
      return RN;
    }
  "

  prospect_Cparams <- list(
    "4" = c("Ccab", "Cw", "Cm"),
    "5" = c("Ccab", "Ccar", "Cw", "Cm"),
    "5b" = c("Ccab", "Ccar", "Cbrown", "Cw", "Cm"),
    "d" = c("Ccab", "Ccar", "Canth", "Cbrown", "Cw", "Cm"),
    "pro" = c("Ccab", "Ccar", "Canth", "Cbrown", "Cw", "Cprot", "Ccbc")
  )[[prospect_version]]
  prospect_args <- paste(lapply(prospect_Cparams, sprintf, fmt = "real %s"), collapse = ", ")
  n_col <- length(prospect_Cparams)
  ccmat <- paste(
    glue::glue("cc[{seq_along(prospect_Cparams)},1] = {prospect_Cparams} / N;"),
    collapse = "\n"
  )

  prospect_func <- sprintf(
    "
    vector prospect(
      real N, %s, data matrix kmat,
      data vector talf, data vector ralf,
      data vector t12, data vector r12,
      data vector t21, data vector r21
    ) {
      matrix[%d, 1] cc;
      int nwl = size(talf);
      %s
      vector[nwl] k = to_vector(kmat * cc);
      vector[nwl] result;
      for (i in 1:nwl) {
        result[i] = gpm(N, k[i], talf[i], ralf[i], t12[i], r12[i], t21[i], r21[i]);
      }
      return result;
    }
    ", prospect_args, n_col, ccmat
  )

  funcblock <- paste(
    "functions {",
    expint_func, transfun_func, gpm_func, prospect_func,
    "}", collapse = "\n"
  )

  datablock <- sprintf(
    "
    data {
      int<lower=0> nwl;
      int<lower=0> nobs;
      // array[nobs] row_vector[nwl] obs;
      array[nwl, nobs] real obs;
      vector[nwl] talf;
      vector[nwl] t12;
      vector[nwl] t21;
      matrix[nwl, %d] kmat;
    }

    transformed data {
      vector[nwl] ralf = 1 - talf;
      vector[nwl] r12 = 1 - t12;
      vector[nwl] r21 = 1 - t21;
    }
    ", n_col
  )
  paramlist <- paste(glue::glue("real<lower=0> {prospect_Cparams};"), collapse = "\n")
  paramblock <- sprintf(
    "
    parameters {
      real<lower=1> N;
      %s
      real<lower=0> rsd;
    }
    ", paramlist
  )

  priorblock <- paste(
    glue::glue("{lapply(prospect_priors[c('N', prospect_Cparams, 'rsd')], dist2stan)};"),
    collapse = "\n"
  )

  modelblock <- sprintf(
    "
    model {
      // Priors
      %s
      // Likelihood
      vector[nwl] mod = prospect(N, %s, kmat,
                                 talf, ralf, t12, r12, t21, r21);
      for (i in 1:nobs) {
        obs[,i] ~ normal(mod, rsd);
      }
    }
    ", priorblock, paste(prospect_Cparams, collapse = ", ")
  )

  paste(
    funcblock,
    datablock,
    paramblock,
    modelblock,
    collapse = "\n"
  )

}
