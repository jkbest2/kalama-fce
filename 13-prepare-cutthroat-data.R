library(tidyverse)

source(here::here("helper-functions.R"))
source(here::here("data-functions.R"))

base_dir <- here::here("results", "cutthroat")
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

## Read preprocessed data
kalr <- read_rds(here::here("data", "kalr.rds"))
trap <- read_rds(here::here("data", "trap.rds"))

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
    temperature_scl = scale_to(temperature, 0.5)
  )

kalr_cth <- kalr |>
  filter(species == "Cutthroat")

pred_df <- fce_pred_df(kalr_cth, trap_op, species, lifestage) |>
  mutate(
    doy = yday(date),
    doy_cent = (doy - mean(doy)) / 7
  ) |>
  left_join(trap_cov, by = join_by(date)) |>
  fce_check_pred_df()

## Tagged fish are released the same day
rel_ref <- tibble(
  collect_date = trap$date,
  release_date = trap$date
)

release_ind <- fce_ind_release(
  kalr_cth, rel_ref,
  mark_col = pit,
  capture_type == "Maiden",
  !is.na(release_site),
  !is.na(pit),
  ## One fish was tagged and released the last day of trap operation, so there
  ## was no way for it to be recaptured. This drops the release but the
  ## individual is still counted as an unmarked capture.
  date != ymd("2025-06-22")
)
recap_ind <- fce_ind_recap(
  kalr_cth,
  mark_col = pit,
  capture_type == "Recapture",
  !is.na(pit)
)

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
  kalr_cth, pred_df,
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
  here::here(base_dir, "cutthroat-prep.rds")
)
