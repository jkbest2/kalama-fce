library(tidyverse)
library(rstan) # `get_num_divergent`
library(posterior)
library(ggdist) # `point_interval`

source(here::here("helper-functions.R"))
source(here::here("data-functions.R"))
source(here::here("post-functions.R"))

## Steelhead results -----------------------------------------------------------
sth_dir <- here::here("results", "steelhead")
sth_prep <- read_rds(here::here(sth_dir, "steelhead-prep.rds"))
sth_fit <- read_rds(here::here(sth_dir, "steelhead-fit.rds"))
sth_diag <- convergence_diagnostics(sth_fit$fit)
sth_post <- as_draws_rvars(sth_fit$fit)
sth_res <- sth_prep$pred_df |>
  select(date, species, op) |>
  mutate(
    pcap = fce_post_pcap(sth_post, sth_fit$data),
    runsize = fce_rarrivals(
      fce_post_rs(sth_post, sth_fit$data),
      sth_prep$unmarked$count
    )
  )
sth_res_summ <- sth_res |>
  point_interval(pcap, runsize, .width = 0.95)
sth_abund <- sth_res |>
  summarize(abundance = rvar_sum(runsize))
sth_abund_df <- as_draws_df(sth_abund$abundance) |>
  rename(abundance = `x`)
sth_abund_summ <- sth_abund |>
  transmute(
    mean = mean(abundance),
    median = median(abundance),
    sd = sd(abundance),
    cv = sd / mean,
    p95.lower = quantile(abundance, 0.025),
    p95.upper = quantile(abundance, 0.975)
  )
sth_coveff <- fce_covariate_effect(
  use_pars = "fixed",
  sth_post, sth_fit$data,
  mod = "pcap"
)
sth_coveff_summ <- sth_coveff |>
  point_interval(post)
write_rds(
  list(
    sth_res_summ = sth_res_summ,
    sth_abund_df = sth_abund_df,
    sth_abund_summ = sth_abund_summ,
    sth_coveff_summ = sth_coveff_summ
  ),
  here::here(sth_dir, "steelhead-summary.rds")
)

## Cutthroat results -----------------------------------------------------------
cth_dir <- here::here("results", "cutthroat")
cth_prep <- read_rds(here::here(cth_dir, "cutthroat-prep.rds"))
cth_fit <- read_rds(here::here(cth_dir, "cutthroat-fit.rds"))
cth_diag <- convergence_diagnostics(cth_fit$fit)
cth_post <- as_draws_rvars(cth_fit$fit)
cth_res <- cth_prep$pred_df |>
  select(date, species, op) |>
  mutate(
    pcap = fce_post_pcap(cth_post, cth_fit$data),
    runsize = fce_rarrivals(
      fce_post_rs(cth_post, cth_fit$data),
      cth_prep$unmarked$count
    )
  )
cth_res_summ <- cth_res |>
  point_interval(pcap, runsize, .width = 0.95)
cth_abund <- cth_res |>
  summarize(abundance = rvar_sum(runsize))
cth_abund_df <- as_draws_df(cth_abund$abundance) |>
  rename(abundance = `x`)
cth_abund_summ <- cth_abund |>
  transmute(
    mean = mean(abundance),
    median = median(abundance),
    sd = sd(abundance),
    cv = sd / mean,
    p95.lower = quantile(abundance, 0.025),
    p95.upper = quantile(abundance, 0.975)
  )
cth_coveff <- fce_covariate_effect(
  use_pars = "fixed",
  cth_post, cth_fit$data,
  mod = "pcap"
)
cth_coveff_summ <- cth_coveff |>
  point_interval(post)
write_rds(
  list(
    cth_res_summ = cth_res_summ,
    cth_abund_df = cth_abund_df,
    cth_abund_summ = cth_abund_summ,
    cth_coveff_summ = cth_coveff_summ
  ),
  here::here(cth_dir, "cutthroat-summary.rds")
)
