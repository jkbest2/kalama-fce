library(tidyverse)
library(rstan)
library(posterior)
library(ggdist) ## For plotting posteriors and distributions
library(patchwork)

source(here::here("helper-functions.R"))
source(here::here("data-functions.R"))
source(here::here("post-functions.R"))

base_dir <- here::here("results", "cutthroat")

## Read prepared data and fit -------------------------------------------------
prep <- read_rds(here::here(base_dir, "cutthroat-prep.rds"))

fit_obj <- read_rds(here::here(base_dir, "cutthroat-fit.rds"))
data <- fit_obj$data

## Check posterior diagnostics ------------------------------------------------
diag_df <- convergence_diagnostics(fit_obj$fit)

## Convert to rvars draws makes it easier to manipulate the posterior
post <- as_draws_rvars(fit_obj$fit) |>
  thin_draws(3)

n_div <- get_num_divergent(fit_obj$fit)
if (n_div > 0) {
  warning("There were ", n_div, " divergent transitions in this fit.")
}
if (n_div > 0 && FALSE) { # Only run these manually
  np <- nuts_params(fit_obj$fit)
  library(bayesplot)
  npsty <- parcoord_style_np("red", 0.4, 1)
  mcmc_parcoord(post, regex_pars = "avail_coef_fixed", np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "avail_coef_random", transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "avail_rand_scale", transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "pcap_coef_fixed", transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "pcap_coef_random", transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, pars = c("pcap_rand_scale", "__lp"), transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "rs_coef_fixed", transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "rs_coef_random", transform = scale, np = np, np_style = npsty)
  mcmc_parcoord(post, regex_pars = "rs_rand_scale", transform = scale, np = np, np_style = npsty)
}

## Plots
avail_plot(post, prep, data, "cutthroat", base_dir)
# Warnings about missing data are from NA's substituted in when trap is not
# operating. Alternatively, can use `plot_noop = TRUE` to show pcap going to
# zero during trap closures.
pcap_plot(
  post, prep, data,
  "cutthroat", base_dir,
  .width = 0.95,
  plot_noop = FALSE, noop_bg = FALSE
) +
  scale_size_continuous(breaks = seq(1, 13, 4))
runsize_plot(
  post, prep, data,
  "cutthroat", base_dir,
  .width = 0.95,
  noop_bg = TRUE
)
fce_abund_plot(post, prep, data, "cutthroat", base_dir)
fce_abundance(post, prep, data) |>
  point_interval(
    .width = c(0.8, 0.95),
    .point = median,
    .interval = qi
  )
fce_abundance(post, prep, data) |>
  mutate(
    mean = mean(abundance),
    sd = sd(abundance),
    cv = sd / mean
  )

## Covariate effects
fce_coveff_plot(
  # use_pars = c("Intercept", "log_discharge_scl", "doy_cent"),
  use_pars = c("Intercept", "discharge_scl", "doy_cent"),
  mod = "pcap",
  post, prep, data,
  incl_prior = TRUE
)

eff_df <- fce_create_eff_df(
  c(
    "factor(doy)",
    "discharge_scl", "doy_cent"
  ),
  "pcap", post, prep, data
)

eff_df |>
  point_interval(`factor(doy)_eff`) |>
  ggplot(aes(x = date, y = `factor(doy)_eff`, ymin = .lower, ymax = .upper)) +
  geom_hline(yintercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  eff_df |>
  point_interval(doy_cent_eff) |>
  ggplot(aes(x = date, y = doy_cent_eff, ymin = .lower, ymax = .upper)) +
  geom_hline(yintercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  eff_df |>
  point_interval(discharge_scl_eff) |>
  ggplot(aes(x = date, y = discharge_scl_eff, ymin = .lower, ymax = .upper)) +
  geom_hline(yintercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  plot_layout(ncol = 1)

eff_df |>
  point_interval(`factor(doy)_eff`) |>
  ggplot(aes(x = date, y = `factor(doy)_eff`, ymin = .lower, ymax = .upper)) +
  geom_hline(yintercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  eff_df |>
  point_interval(discharge_scl_eff) |>
  ggplot(aes(x = discharge, y = discharge_scl_eff, ymin = .lower, ymax = .upper)) +
  geom_hline(yintercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  eff_df |>
  point_interval(doy_cent_eff) |>
  ggplot(aes(x = date, y = doy_cent_eff, ymin = .lower, ymax = .upper)) +
  geom_hline(yintercept = 0, alpha = 0.25, linetype = "dashed") +
  geom_ribbon(alpha = 0.4) +
  geom_line() +
  scale_x_date(date_breaks = "1 month", date_labels = "%B") +
  plot_layout(ncol = 1)

## How much variation in the random effects? If posterior is mostly close to
## zero might consider eliminating that random effect.
## Note that if availability smoothers are completely separate these parameters
## cannot be compared directly to each other.
fce_randscale_plot(post, prep, data, mod = "avail", incl_prior = TRUE)
fce_randscale_plot(post, prep, data, mod = "pcap", incl_prior = TRUE)
fce_randscale_plot(post, prep, data, mod = "rs", incl_prior = TRUE)
