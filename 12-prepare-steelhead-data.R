library(tidyverse)

source(here::here("helper-functions.R"))
source(here::here("data-functions.R"))

base_dir <- here::here("results", "steelhead")
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

## Read preprocessed data
kalr <- read_rds(here::here("data", "kalr.rds"))
trap <- read_rds(here::here("data", "trap.rds"))

## Exclude outlier that took 39 days to be recaptured?
exclude_outlier <- FALSE

## Data preparation -----------------------------------------------------------
trap_op <- trap |>
  select(date, op)
trap_cov <- trap |>
  select(date, discharge, temperature) |>
  mutate(
    log_discharge_scl = scale_to(log(discharge), 0.5),
    ## Use linear interpolation to add values to temperature time series. All
    ## missing values are during trap outages. Because pcap during trap outages
    ## is fixed at 0, these interpolated values don't enter the likelihood in
    ## any way. The interpolation is done because Stan doesn't accept missing
    ## values.
    temp_interp = is.na(temperature),
    temperature = approx(
      x = seq_along(temperature),
      y = temperature,
      xout = seq_along(temperature)
    )$y,
    temperature_scl = scale_to(temp_interp, 0.5)
  )

kalr_sth <- kalr |>
  filter(species == "Steelhead")

pred_df <- fce_pred_df(kalr_sth, trap_op, species, lifestage) |>
  mutate(
    doy = yday(date),
    doy_cent = (yday(date) - mean(yday(date))) / 7 # Centered and week-scaled
  ) |>
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

if (exclude_outlier) {
  outl_ind <- release_ind |>
    right_join(recap_ind, by = join_by(pit)) |>
    mutate(travel_time = recapture_date - release_date) |>
    select(release_date, recapture_date, release_site, pit, travel_time) |>
    filter(travel_time > 20)
  release_ind <- anti_join(release_ind, outl_ind, by = join_by(pit))
  recap_ind <- anti_join(recap_ind, outl_ind, by = join_by(pit))
}

release_df <- fce_release_df(release_ind)
recap_df <- fce_recap_df(
  release_ind, recap_ind, pred_df,
  recap_window_fun = fce_recap_window_site_obs,
  drop_notrap = TRUE
) |>
  mutate(
    ## Add 0/1 columns for each release site to facilitate individual smoothers
    ## for availability from each release site
    lowersite = as.numeric(release_site == "Lower Site"),
    yellowgate = as.numeric(release_site == "Yellow Gate"),
    jackscreek = as.numeric(release_site == "Jacks Creek"),
    ## Add duplicate travel_time columns to avoid mgcv complaining about using
    ## the same variable twice
    tt_lower = travel_time,
    tt_yellow = travel_time,
    tt_jacks = travel_time
  )
unmarked <- fce_unmarked(
  kalr_sth, pred_df,
  capture_type == "Maiden",
  drop_notrap = TRUE
)

write_rds(
  list(
    release_df = release_df,
    recap_df = recap_df,
    unmarked = unmarked,
    pred_df = pred_df
  ),
  here::here(base_dir, "steelhead-prep.rds")
)
