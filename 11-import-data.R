library(tidyverse)

kalr <- read_rds(here::here("rawdata", "kalr_fish.rds")) |>
  as_tibble() |>
  rename(
    date = Date,
    species = Species,
    pit = PIT_Code,
    capture_type = CaptureType,
    count = Count,
    length = Length,
    release_site = ReleaseUpstream
  ) |>
  mutate(
    species = factor(species),
    lifestage = factor("Smolt"), # Placeholder lifestage
    capture_type = factor(capture_type),
    release_site = factor(
      release_site,
      levels = c("Lower Site", "Yellow Gate", "Jacks Creek")
    )
  )
trap <- read_rds(here::here("rawdata", "kalr_trap.rds")) |>
  as_tibble() |>
  rename(
    date = Date,
    op = TrapStatus,
    discharge = Disc
  ) |>
  mutate(op = as.logical(op))

write_rds(kalr, here::here("data", "kalr.rds"))
write_rds(trap, here::here("data", "trap.rds"))
