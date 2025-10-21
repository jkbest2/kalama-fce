library(tidyverse)

source(here::here("helper-functions.R"))
source(here::here("data-functions.R"))

base_dir <- here::here("results", "steelhead")
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

## Read preprocessed data
kalr <- read_rds(here::here("data", "kalr.rds"))
trap <- read_rds(here::here("data", "trap.rds"))

## Data preparation -----------------------------------------------------------
trap_op <- trap |>
  select(date, op)
trap_cov <- trap |>
  select(date, discharge) |>
  mutate(discharge_scl = scale_to(discharge, 0.5))

kalr_sth <- kalr |>
  filter(species == "Steelhead")

pred_df <- fce_pred_df(kalr_sth, trap_op, species, lifestage) |>
  mutate(doy = yday(date)) |>
  left_join(trap_cov, by = join_by(date)) |>
  fce_check_pred_df()

## Tagged fish are released the same day
rel_ref <- tibble(
  collect_date = trap$date,
  release_date = trap$date
)

release_ind <- fce_ind_release(
  kalr_sth, rel_ref,
  mark_col = pit,
  capture_type == "Maiden",
  !is.na(release_site),
  !is.na(pit)
)
recap_ind <- fce_ind_recap(
  kalr_sth,
  mark_col = pit,
  capture_type == "Recapture",
  !is.na(pit)
)

rw_df2 <- right_join(release_ind, recap_ind, by = join_by(pit)) |>
  mutate(travel_time = as.numeric(recapture_date - release_date))

ggplot(rw_df2, aes(x = factor(travel_time), fill = release_site)) +
  geom_bar(aes(y = after_stat(count))) +
  facet_wrap(~release_site, ncol = 1)
