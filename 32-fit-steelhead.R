library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 4L) # One core per chain

source(here::here("helper-functions.R"))
source(here::here("data-functions.R"))

base_dir <- here::here("results", "steelhead")

## Read preprocessed data
sth_prep <- read_rds(here::here(base_dir, "steelhead-prep.rds"))

## Model specification --------------------------------------------------------
avail <- list(
  formula = ~ 0 +
    s(tt_lower, k = 4, bs = "gp", m = c(3, 7), by = lowersite) +
    s(tt_yellow, k = 8, bs = "gp", m = c(3, 7), by = yellowgate) +
    s(tt_jacks, k = 8, bs = "gp", m = c(3, 7), by = jackscreek),
  priors = fce_priors(
    fixed = list(
      `tt_lower_null` = normal(0, 1e-2),
      `tt_yellow_null` = normal(0, 1e-2),
      `tt_jacks_null` = normal(0, 1e-2)
    ),
    random = list(
      `tt_lower` = normal(0, 4),
      `tt_yellow` = normal(0, 4),
      `tt_jacks` = normal(0, 4)
    )
  )
)
## pcap does not appear to change much over time in this case. Including a
## smoother (`s()` term in the formula) as a random effect causes the
## associated variance term to "want" to be near zero, inducing substantial
## curvature in the posterior and divergent transitions. These could be
## problematic becuase they indicate an area of the posterior that cannot be
## sampled efficiently, potentially biasing the results. Divergences can be
## eliminated by treating the smoother as a fixed effect, fixing its marginal
## variance, or by eliminating it altogether. The former solution leaves the
## posterior sensitive to the chosen variance and the latter removes some
## flexibility.
pcap <- list(
  formula = ~ 1 + log_discharge_scl * temperature_scl,
  priors = fce_priors(
    fixed = list(
      Intercept = normal(-2, 1),
      log_discharge_scl = normal(0, 2.5),
      temperature_scl = normal(0, 2.5),
      `log_discharge_scl:temperature_scl` = normal(0, 2.5)
    ),
    random = list()
  )
)
rs <- list(
  formula = ~ 1 +
    factor(doy) +
    s(doy, k = 64, bs = "gp", m = c(3, 7)),
  priors = fce_priors(
    fixed = list(
      Intercept = normal(3, 1),
      doy_null = normal(0, 1e-2)
    ),
    random = list(
      `factor(doy)` = normal(0, 0.5),
      doy = normal(0, 4)
    )
  )
)

## Create data structure for Stan model. Safe to ignore warning "by variables
## lowersite, yellowgate, jackscreek are not factors"
data <- fce_data(
  avail$formula, avail$priors,
  pcap$formula, pcap$priors,
  rs$formula, rs$priors,
  sth_prep$release_df, sth_prep$recap_df,
  sth_prep$unmarked, sth_prep$pred_df
)

fit <- stan(
  here::here("fce-model.stan"),
  data = data,
  pars = c(
    "avail_coef_fixed", "avail_coef_random", "avail_rand_scale",
    "pcap_coef_fixed", "pcap_coef_random", "pcap_rand_scale",
    "rs_coef_fixed", "rs_coef_random", "rs_rand_scale",
    ## Derived quantities
    "avail", "pcap", "log_lambda",
    ## Posterior predictive quantities
    "recap_rep", "lost_rep", "unm_rep",
    ## Log-liklihood, used for `loo` model selection
    "log_lik"
  ),
  chains = 4, iter = 4000, warmup = 2000,
  init_r = 1,
  control = list(adapt_delta = 0.8)
)

write_rds(
  list(data = data, fit = fit),
  here::here(base_dir, "steelhead-fit.rds")
)
