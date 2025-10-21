library(tidyverse)
library(rstan)
library(bayesplot)

species <- "steelhead"
## species <- "cutthroat"

base_dir <- here::here("results", species)

fit_obj <- read_rds(here::here(base_dir, paste0(species, "-fit.rds")))

## Extract NUTS parameters
np <- nuts_params(fit_obj$fit)

## Make divergent transitions wider, less transparent, and red
npsty <- parcoord_style_np("red", 0.4, 1)

## The `mcmc_parcoord` function needs exact parameter names in the `pars`
## argument, including any indices. Use `regex_pars` to get all parameters with
## a given name. The function will error if only a single parameter is
## included. In that case consider combining it with another parameter (e.g.
## *_coef_random and *_rand_scale parameters). Another option is to include
## "__lp" as a parameter. Using `transform = scale` puts all the parameters on
## the same scale.
mcmc_parcoord(
  post,
  regex_pars = "avail_coef_fixed",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "aval_coef_random",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "avail_rand_scale",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "pcap_coef_fixed",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "pcap_coef_random",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "pcap_rand_scale",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "rs_coef_fixed",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "rs_coef_random",
  transform = scale,
  np = np, np_style = npsty
)
mcmc_parcoord(
  post,
  regex_pars = "rs_rand_scale", transform = scale, np = np, np_style = npsty
)
