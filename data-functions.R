library(tidyverse)
library(mgcv) # Smoothers
library(extraDistr)
library(distributional) # dist_normal and friends

source(here::here("helper-functions.R"))

## Finds the trap season for each year and fills in all years. Covariates
## should be joined by date to ensure they line up correctly.
fce_pred_df <- function(data, trap_op, ...) {
  expansion <- data |>
    select(species, lifestage, ...) |>
    distinct()

  ## Using trap_season allows us to look at a different season in each year,
  ## depending on when the trap was run. In some cases this can and should be
  ## extended.
  trap_op |>
    unnest(date) |>
    expand_grid(expansion) |>
    mutate(
      op = replace_na(op, FALSE),
      idx = seq_along(date)
    ) |>
    select(idx, date, op, ...)
}

fce_check_pred_df <- function(pred_df) {
  chk <- pred_df |>
    summarize(across(everything(), \(x) any(is.na(x)))) |>
    pivot_longer(everything()) |>
    filter(value)
  if (nrow(chk) == 0) {
    return(invisible(pred_df))
  }
  stop(
    "The following columns include missing entries: ",
    paste(chk$name, collapse = ", ")
  )
}

fce_ind_release <- function(data, release_sched, mark_col = pit, ...) {
  data |>
    filter(
      !is.na({{ mark_col }}),
      ...
    ) |>
    left_join(release_sched, by = join_by(date == collect_date)) |>
    select(
      release_date, species, lifestage,
      release_site, {{ mark_col }}, count
    )
}


fce_ind_release_summ <- function(release, trap_op) {
  release |>
    group_by(release_date, species, lifestage) |>
    summarize(
      count = sum(count),
      .groups = "drop"
    ) |>
    left_join(trap_op, by = join_by(release_date == date)) |>
    arrange(release_date, species, lifestage) |>
    select(release_date, species, lifestage, count, doy) |>
    mutate(release_idx = seq_along(release_date))
}

fce_ind_recap <- function(data, mark_col = pit, ...) {
  data |>
    filter(
      !is.na({{ mark_col }}),
      ...
    ) |>
    select(recapture_date = date, {{ mark_col }})
}

## This only makes sense to use on data frames where the release and recapture
## groups have already been separated.
fce_get_dup_mark <- function(df_ind, mark_col) {
  df_ind |>
    summarize(n = n(), .by = {{ mark_col }}) |>
    filter(n > 1)
}
fce_check_dup_mark <- function(df_ind, mark_col) {
  df_dup <- fce_get_dup_mark(df_ind, {{ mark_col }})
  if (nrow(df_dup) > 0) {
    stop("Duplicated marks found")
  }
  df_ind
}

## FIXME
fce_unmarked <- function(data, pred_df, ..., drop_notrap = FALSE) {
  cdat <- data |>
    filter(...) |>
    select(date, species, lifestage, count) |>
    summarize(
      count = sum(count),
      .by = c(date, species, lifestage)
    )
  unm <- pred_df |>
    left_join(cdat, by = join_by(date, species, lifestage)) |>
    mutate(count = replace_na(count, 0))

  if (drop_notrap) {
    unm <- unm |>
      mutate(count = ifelse(!op, 0, count))
  }

  return(unm)
}

fce_release_df <- function(release_ind) {
  release_ind |>
    summarize(
      count = sum(count),
      .by = c(release_date, species, lifestage, release_site)
    ) |>
    mutate(
      release_idx = seq_along(release_date)
    )
}

## Find the first and last day that the trap operates each year
fce_trap_season <- function(data) {
  data |>
    mutate(year = year(date)) |>
    group_by(year) |>
    summarize(
      start = min(date),
      end = max(date)
    )
}

## Extract the operational history of the trap.
fce_trap_op <- function(data) {
  ## Extract the trap season for each year. May be different lengths per year,
  ## will get the row indices later.
  trap_season <- fce_trap_season(data)

  ## Assume that collection is occurring any day when a row is present in the
  ## data set
  trap_coll <- data |>
    summarize(op = TRUE, .by = date)

  trap_season |>
    mutate(date = map2(start, end, seq, by = "1 day")) |>
    select(date) |>
    unnest(date) |>
    left_join(trap_coll, by = join_by(date)) |>
    mutate(
      op = replace_na(op, FALSE),
    )
}

## Generate the data frame with individual mark-recapture histories
fce_ind <- function(release_ind, recap_ind, mark_col = pit) {
  left_join(release_ind, recap_ind, join_by({{ mark_col }})) |>
    select(release_date, recapture_date, species, lifestage, pit, count) |>
    mutate(travel_time = recapture_date - release_date)
}

## Define the species-specific default recapture windows so that it doesn't have
## to be passed around a bunch.
fce_recap_window_default <- function(release_ind, recap_ind) {
  fce_release_df(release_ind) |>
    select(release_date, species) |>
    left_join(
      tribble(
        ~species, ~recap_window,
        "Steelhead", 35,
        "Coho", 35,
        "Chinook", 105
      ),
      by = join_by(species)
    )
}

## Define the recapture window based on the maximum observed recapture delay
fce_recap_window_obs <- function(release_ind, recap_ind) {
  rw_df <- left_join(release_ind, recap_ind, by = join_by(pit)) |>
    filter(!is.na(recapture_date)) |>
    mutate(recap_window = recapture_date - release_date) |>
    summarize(
      recap_window = max(recap_window),
      .by = c(species, lifestage)
    ) |>
    mutate(recap_window = as.numeric(recap_window) + 3)

  fce_release_df(release_ind) |>
    select(release_date, species, lifestage) |>
    left_join(rw_df, by = join_by(species, lifestage))
}

fce_recap_window_site_obs <- function(release_ind, recap_ind) {
  rw_df <- left_join(release_ind, recap_ind, by = join_by(pit)) |>
    filter(!is.na(recapture_date)) |>
    mutate(recap_window = recapture_date - release_date) |>
    summarize(
      recap_window = max(recap_window),
      .by = c(species, lifestage, release_site)
    ) |>
    mutate(recap_window = as.numeric(recap_window) + 3)

  fce_release_df(release_ind) |>
    select(release_date, species, lifestage, release_site) |>
    left_join(rw_df, by = join_by(species, lifestage, release_site))
}

fce_ragged_index <- function(lengths, additional = 0) {
  cumsum(c(1, lengths + additional))
}

fce_recap_df <- function(release_ind, recap_ind,
                         pred_df,
                         recap_window_fun = fce_recap_window_site_obs,
                         drop_notrap = FALSE) {
  ## FIXME add checks for required column names!

  ## Summarize release information. Probably already did this outside the
  ## function, but need the individual histories so it's not a big deal to do it
  ## again. Also attach the recapture window for each release.
  release_df <- fce_release_df(release_ind) |>
    left_join(recap_window_fun(release_ind, recap_ind),
      by = join_by(release_date, species, lifestage, release_site)
    )

  ## Create the skeleton for recapture information; fill in the potential dates
  ## of recapture based on the release date and the release window
  recap_df <- release_df |>
    mutate(
      recapture_date = map2(
        release_date, recap_window,
        \(rd, rw) c(rd + 1:rw)
      )
    ) |>
    select(release_date, release_site, species, lifestage, recapture_date) |>
    unnest(recapture_date) |>
    mutate(
      travel_time = as.numeric(recapture_date - release_date)
    ) |>
    left_join(
      select(pred_df, date, idx, op, species, lifestage),
      by = join_by(recapture_date == date, species, lifestage),
    ) |>
    mutate(op = replace_na(op, FALSE))

  ## Summarize the recapture information as number recaptured each recapture day
  recap_summ <- left_join(release_ind, recap_ind, by = join_by(pit)) |>
    summarize(
      count = sum(count),
      .by = c(species, lifestage, release_date, recapture_date)
    )

  ## Associate the recapture counts by day with the specific days, filling in
  ## zeros where necessary
  rdf <- recap_df |>
    left_join(
      recap_summ,
      by = join_by(
        species,
        lifestage,
        release_date,
        recapture_date
      )
    ) |>
    mutate(
      count = replace_na(count, 0)
    )

  if (drop_notrap) {
    rdf <- rdf |>
      mutate(count = ifelse(!op, 0, count))
  }

  ## Throw a warning if any recaptures occur on days when the trap is marked
  ## as not operational. Don't throw an error so that we can return the data
  ## frame and its rows for debugging. Should probably be an error though.
  chk <- !rdf$op & rdf$count > 0
  if (any(chk)) {
    warning(
      "Recaptures detected when trap not operating in rows ",
      paste(which(chk), collapse = ", ")
    )
  }

  return(rdf)
}

fce_num_lost <- function(release_df, recap_df) {
  n_recap <- recap_df |>
    summarize(
      recap_count = sum(count),
      .by = c(release_date, species, lifestage)
    )
  left_join(release_df, n_recap,
    by = join_by(release_date, species, lifestage)
  ) |>
    mutate(
      n_lost = count - recap_count
    ) |>
    pluck("n_lost")
}

fce_knot_vec <- function(inner, outer, expand = 1.1) {
  inr <- range(inner)
  inr <- inr + c(-1, 1) * diff(inr) * (expand - 1)
  outr <- range(outer)
  outr <- outr + c(-1, 1) * diff(outr) * (expand - 1)
  c(outr[1], inr, outr[2])
}

fce_parse_formula <- function(formula, data, knots = NULL) {
  if (attr(terms(formula), "response")) {
    warning("Must use a one-sided formula, response removed")
    formula <- formula[-2]
  }

  gp <- interpret.gam(formula)

  ## Check that `by` variables are factors. Occasionally may be numeric for
  ## varying-coefficient models. Character variables will not be coerced to
  ## factor and will result in "by variable not found" errors.
  nonfct_by <- map_chr(gp$smooth.spec, pluck, "by") |>
    discard(
      \(x) x == "NA"
    ) |>
    discard(
      \(x) isa(data[[x]], "factor") || isa(data[[x]], "numeric")
    )
  if (length(nonfct_by) > 0) {
    warning(
      "by variables ",
      paste(nonfct_by, collapse = " "),
      " are not factors or numeric"
    )
  }
  ## Check that all predictors are present in the provided data
  if (!all(gp$pred.names %in% names(data))) {
    np <- gp$pred.names[!gp$pred.names %in% names(data)]
    stop(
      ## Provide the formula so that it can be matched with the process
      "Formula provided is ", formula, " but ",
      paste(np, collapse = ", "), " is/are not present in data provided."
    )
  }

  sm <- lapply(gp$smooth.spec,
    smoothCon,
    data = data,
    knots = knots,
    absorb.cons = TRUE,
    scale.penalty = FALSE,
    diagonal.penalty = TRUE
  )
  list(
    linear = gp$pf,
    smooth = sm,
    data = data,
    pred_names = gp$pred.names
  )
}

fce_construct_mm <- function(pform, data = NULL,
                             check_fullrank = TRUE) {
  if (is.null(data)) {
    data <- pform$data
  } else {
    ## Make sure that all required columns are present
    if (!all(pform$pred_names %in% names(data))) {
      mc <- setdiff(pform$pred_names, names(data))
      stop("data is missing column(s): ", mc)
    }
    ## Makes sure that any columns that are factors in the original data set are
    ## also factors with the same set of levels in the new data.
    fct_cols <- names(Filter(is.factor, pform$data))
    data <- data |>
      mutate(
        across(
          all_of(fct_cols),
          \(col) factor(col, levels = levels(pform$data[[cur_column()]]))
        )
      )
  }

  lin_mm <- model.matrix(pform$linear, data = data, na.action = na.fail)
  ## The intercept in attr(terms(pform$linear)) has a numeric value of zero, so
  ## concatenate it at the start of the vector and then add one to the indices
  ## used below
  lin_labs <- c("Intercept", attr(terms(pform$linear), "term.labels")) |>
    set_names()
  lin_scl <- lin_labs[attr(lin_mm, "assign") + 1]

  sm <- flatten(pform$smooth)
  sm_lab <- map_chr(sm, pluck, "label")
  sm_term <- map(sm, pluck, "term")
  sm_term <- map_chr(sm_term, paste0, collapse = ":")
  sm_pnull <- map(sm, \(s) c(rep(FALSE, s$rank), rep(TRUE, s$null.space.dim)))

  sm_scl <- map2(
    sm_term, sm_pnull,
    \(t, p) paste0(t, ifelse(p, "_null", ""))
  )

  sm_mm <- map(sm, PredictMat, data = data)
  sm_mm <- map2(
    sm_mm, sm_lab,
    \(s, l) {
      colnames(s) <- rep(l, ncol(s))
      s
    }
  )

  ## Need `.init = vector()` for the case where there are no smoothers to
  ## `cbind`
  mm <- cbind(lin_mm, reduce(sm_mm, cbind, .init = vector()))

  if (check_fullrank && Matrix::rankMatrix(mm) < ncol(mm)) {
    warning("Model matrix is not full rank")
  }

  scl_idx <- flatten_chr(c(lin_scl, sm_scl))
  attr(mm, "scale_index") <- scl_idx

  mm
}

## Create the mgcv smoother objects that we can then use to get design matrices
fce_spl <- function(spec, data, knots = NULL) {
  smoothCon(
    spec,
    data = data,
    knots = knots,
    diagonal.penalty = TRUE, # Don't need to pass penalty matrix to Stan
    scale.penalty = FALSE
  )
}

normal <- function(mu = 0, sigma = 1) {
  stopifnot(is.finite(mu), is.finite(sigma), sigma > 0)
  structure(
    list(mu = mu, sigma = sigma),
    class = "normal_prior"
  )
}

fce_priors <- function(fixed, random) {
  fnames <- names(fixed)
  rnames <- names(random)
  if (length(intersect(fnames, rnames)) != 0) {
    stop("Cannot specify priors in both fixed and random")
  }
  if (any(c(fnames, rnames) == "")) {
    stop("All elements must be named")
  }
  walk(
    c(fixed, random),
    ~ all(is(., "normal_prior")) ||
      stop("All priors should be specified using the \`normal\` function")
  )
  walk(random, ~ .$mu == 0 || stop("random hyperpriors must have mean 0"))

  structure(
    list(
      prior_mean = map(fixed, pluck, "mu"),
      prior_sd = map(fixed, pluck, "sigma"),
      hprior_sd = map(random, pluck, "sigma")
    ),
    class = "fce_priors",
    terms = c(fnames, rnames)
  )
}

fce_basis <- function(pform, priors, data) {
  isa(priors, "fce_priors") ||
    stop("priors must be constructed using the \`fce_priors\` function")

  mm <- fce_construct_mm(pform, data, FALSE)

  ## Construct the coefficient scaling parameter vector
  if (!all(unique(attr(mm, "scale_index") %in% attr(priors, "terms")))) {
    miss <- setdiff(attr(mm, "scale_index"), attr(priors, "terms"))
    stop("Prior missing for ", paste(miss, collapse = ", "))
  }

  ## Separate out the fixed effects design matrix and priors
  idx_fixed <- attr(mm, "scale_index") %in% names(priors$prior_mean)
  ## Extract the fixed effects columns
  basis_fixed <- mm[, idx_fixed, drop = FALSE]
  ## Be sure to retain the scale associations
  attr(basis_fixed, "scale_index") <- attr(mm, "scale_index")[idx_fixed]
  prior_mean <- priors$prior_mean[
    names(priors$prior_mean) %in% attr(basis_fixed, "scale_index")
  ] |> unlist()
  prior_sd <- priors$prior_sd[
    names(priors$prior_sd) %in% attr(basis_fixed, "scale_index")
  ] |> unlist()
  prior_idx <- match(
    attr(basis_fixed, "scale_index"),
    names(prior_mean)
  )

  all(prior_idx == match(attr(basis_fixed, "scale_index"), names(prior_sd))) ||
    stop("Prior means and sds not indexed the same")

  ## Separate out the random effects design matrix and scales
  idx_rand <- attr(mm, "scale_index") %in% names(priors$hprior_sd)
  ## Extract the random effects columns
  basis_rand <- mm[, idx_rand, drop = FALSE]
  ## Be sure to retain the scale associations
  attr(basis_rand, "scale_index") <- attr(mm, "scale_index")[idx_rand]
  hprior_rand <- priors$hprior_sd[
    names(priors$hprior_sd) %in% attr(basis_rand, "scale_index")
  ] |> unlist()
  if (is.null(hprior_rand)) {
    hprior_rand <- numeric(0)
  }
  hprior_idx <- match(
    attr(basis_rand, "scale_index"),
    names(hprior_rand)
  )

  stopifnot(
    ## Make sure you use all the columns
    ncol(mm) == ncol(basis_fixed) + ncol(basis_rand),
    ## And rows
    nrow(mm) == nrow(basis_fixed),
    ## And both bases have the same number of rows
    nrow(mm) == nrow(basis_rand),
    ## And we have priors and hyperpriors assigned for each column/coefficient
    ncol(basis_fixed) == length(prior_idx),
    ncol(basis_rand) == length(hprior_idx)
  )

  ## Use `array` here to ensure that single-element entries will not be treated
  ## as scalars by Stan
  structure(
    list(
      basis_fixed = basis_fixed,
      prior_mean = array(prior_mean, dim = length(prior_mean)),
      prior_sd = array(prior_sd, dim = length(prior_sd)),
      prior_idx = array(prior_idx, dim = length(prior_idx)),
      basis_rand = basis_rand,
      hprior_sd = array(hprior_rand, dim = length(hprior_rand)),
      hprior_idx = array(hprior_idx, dim = length(hprior_idx))
    ),
    terms = attr(scale, "terms"),
    class = "fce_basis"
  )
}

## Full model functions --------------------------------------------------------
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Construct the data argument for the full model
##' @param avail_formula One-sided formula specifying the covariates and/or
##'   smoothers for the availability model.
##' @param avail_priors Priors for parameters in availability model constructed
##'   using the \code{fce_priors} function.
##' @param pcap_formula One-sided formula specifying the covariates and/or
##'   smoothers for the probability of capture model.
##' @param pcap_priors Priors for parameters in probability of capture model
##'   constructed using the \code{fce_priors} function.
##' @param rs_formula One-sided formula specifying the covariates and/or
##'   smoothers for the run size model.
##' @param rs_priors Priors for parameters in run size model constructed
##'   using the \code{fce_priors} function.
##' @param release_df Data frame with one row per release
##' @param recap_df Data frame with one row per recapture opportunity, including
##'   all covariates used in \code{avail_formula}
##' @param unmarked Data frame with one row per unmarked capture day, including
##'   all covariates used in \code{rs_formula}
##' @param pred_df Data frame with all covariates used in \code{pcap_formula}
##    and \code{rs_formula} as well as a logical column indicating when the
##    trap was operating
##' @param knots Knot locations for the \code{pcap} model. Not required.
##' @return A named list suitable for use in the FCE model.
##' @author John Best
##' @export
fce_data <- function(
    avail_formula, avail_priors,
    pcap_formula, pcap_priors,
    rs_formula, rs_priors,
    release_df, recap_df, unmarked, pred_df,
    knots = NULL) {
  ## Add an index for the availability model
  recap_df <- recap_df |>
    mutate(avail_idx = seq_along(release_date))

  ## Extract the start index of each release
  # Filter out no-trap days
  recap_trap_df <- recap_df |>
    filter(op)

  unm_trap_df <- unmarked |>
    filter(op)

  avail_idx <- recap_df |>
    summarize(
      n_avail = n(),
      .by = c(species, lifestage, release_date)
    ) |>
    pluck("n_avail") |>
    fce_ragged_index()

  ## This index drops no-trap days. We'll need the index column to index into
  ## pred_df-based design matrices
  recap_idx <- recap_trap_df |>
    summarize(
      n_recap = n(),
      .by = c(species, lifestage, release_date, release_site)
    ) |>
    pluck("n_recap") |>
    fce_ragged_index()

  ## But we need the no-trap days when looking at availability
  avail_pform <- fce_parse_formula(avail_formula, recap_df)
  avail_basis <- fce_basis(avail_pform, avail_priors, recap_df)

  ## Construct smoother and design matrix for probability of capture
  pcap_pform <- fce_parse_formula(pcap_formula, pred_df, knots)
  pcap_basis <- fce_basis(pcap_pform, pcap_priors, pred_df)

  ## Construct smoother and design matrix for runsize estimation
  rs_pform <- fce_parse_formula(rs_formula, pred_df)
  rs_basis <- fce_basis(rs_pform, rs_priors, pred_df)

  structure(
    list(
      ## Sizes
      N_release = nrow(release_df),
      N_pred = nrow(pred_df),
      N_avail_obs = nrow(avail_basis$basis_fixed),
      N_rec_obs = nrow(recap_trap_df),
      N_unm_obs = nrow(filter(unmarked, op)),

      ## Is the trap operating?
      op = pred_df$op,

      ## Release information
      num_released = release_df$count,
      num_lost = fce_num_lost(release_df, recap_df),

      ## Availability information
      avail_idx = avail_idx,
      avail_rec_idx = recap_trap_df$avail_idx,

      ## Recapture information
      rec_idx = recap_idx,
      pred_rec_idx = recap_trap_df$idx,
      num_recaptured = recap_trap_df$count,

      ## Unmarked capture observations
      unm_idx = unm_trap_df$idx,
      unm_count = unm_trap_df$count,

      ## Availabilitiy model dimensions
      N_avail_fixed = ncol(avail_basis$basis_fixed),
      N_avail_prior = length(avail_basis$prior_mean),
      N_avail_random = ncol(avail_basis$basis_rand),
      N_avail_hprior = length(avail_basis$hprior_sd),
      ## Availability model
      avail_fixed = avail_basis$basis_fixed,
      avail_prior_sd = avail_basis$prior_sd,
      avail_prior_mean = avail_basis$prior_mean,
      avail_prior_idx = avail_basis$prior_idx,
      avail_random = avail_basis$basis_rand,
      avail_hprior = avail_basis$hprior_sd,
      avail_hprior_idx = avail_basis$hprior_idx,

      ## Probability of capture model dimensions
      N_pcap_fixed = ncol(pcap_basis$basis_fixed),
      N_pcap_prior = length(pcap_basis$prior_mean),
      N_pcap_random = ncol(pcap_basis$basis_rand),
      N_pcap_hprior = length(pcap_basis$hprior_sd),
      ## Probability of capture model data
      pcap_fixed = pcap_basis$basis_fixed,
      pcap_prior_sd = pcap_basis$prior_sd,
      pcap_prior_mean = pcap_basis$prior_mean,
      pcap_prior_idx = pcap_basis$prior_idx,
      pcap_random = pcap_basis$basis_rand,
      pcap_hprior = pcap_basis$hprior_sd,
      pcap_hprior_idx = pcap_basis$hprior_idx,

      ## Runsize model dimensions
      N_rs_fixed = ncol(rs_basis$basis_fixed),
      N_rs_prior = length(rs_basis$prior_mean),
      N_rs_random = ncol(rs_basis$basis_rand),
      N_rs_hprior = length(rs_basis$hprior_sd),
      ## Probability of capture model data
      rs_fixed = rs_basis$basis_fixed,
      rs_prior_sd = rs_basis$prior_sd,
      rs_prior_mean = rs_basis$prior_mean,
      rs_prior_idx = rs_basis$prior_idx,
      rs_random = rs_basis$basis_rand,
      rs_hprior = rs_basis$hprior_sd,
      rs_hprior_idx = rs_basis$hprior_idx,

      ## Mix IS flag
      mixis = 0
    ),
    class = "fce_data",
    avail_pform = avail_pform,
    avail_basis = avail_basis,
    pcap_pform = pcap_pform,
    pcap_basis = pcap_basis,
    rs_pform = rs_pform,
    rs_basis = rs_basis
  )
}

## Check the data object for consistent dimensions etc.
fce_validate_data <- function(data) {
  ## TODO Add tests for prior and hyperprior indices to ensure all are used
  stopifnot(
    ## Sizes
    data$N_release > 0,
    data$N_pred > 0,
    data$N_avail_obs > 0,
    data$N_rec_obs > 0,
    data$N_unm_obs > 0,

    ## Trap operation
    length(data$op) == data$N_pred,
    all(data$op %in% c(FALSE, TRUE)),

    ## Release information
    length(data$num_released) == data$N_release,
    all(data$num_released >= 1),
    length(data$num_lost) == data$N_release,
    all(data$num_lost >= 0),
    all(data$num_lost <= data$num_released),

    ## Availability information
    length(data$avail_idx) == data$N_release + 1,
    all(data$avail_idx > 0),
    !is.unsorted(data$avail_idx),
    length(data$avail_rec_idx) == data$N_rec_obs,
    all(data$avail_rec_idx >= 1 & data$avail_rec_idx <= data$N_avail_obs),

    ## Recapture information
    length(data$rec_idx) == data$N_release + 1,
    all(data$rec_idx > 0),
    !is.unsorted(data$rec_idx),
    length(data$pred_rec_idx) == data$N_rec_obs,
    all(data$pred_rec_idx >= 1 & data$pred_rec_idx <= data$N_pred),
    length(data$num_recaptured) == data$N_rec_obs,
    all(data$num_recaptured >= 0),

    ## Unmarked capture observations
    length(data$unm_idx) == data$N_unm_obs,
    all(data$unm_idx >= 1 & data$unm_idx <= data$N_pred),
    length(data$unm_count) == data$N_unmarked_obs,
    all(data$unm_count >= 0),

    ## Availability model dimensions
    data$N_avail_fixed >= 0,
    data$N_avail_prior >= 0,
    data$N_avail_random >= 0,
    data$N_avail_hprior >= 0,
    ## Availability model
    nrow(data$avail_fixed) == data$N_avail_obs,
    ncol(data$avail_fixed) == data$N_avail_fixed,
    length(data$avail_prior_sd) == data$N_avail_prior,
    all(data$avail_prior_sd > 0),
    length(data$avail_prior_mean) == data$N_avail_prior,
    length(data$avail_prior_idx) == data$N_avail_fixed,
    all(data$avail_prior_idx > 0),
    all(data$avail_prior_idx <= data$N_avail_fixed),
    nrow(data$avail_random) == data$N_avail_obs,
    ncol(data$avail_random) == data$N_avail_random,
    length(data$avail_hprior) == data$N_avail_hprior,
    all(data$avail_hprior > 0),
    length(data$avail_hprior_idx) == data$N_avail_random,
    all(data$avail_hprior_idx > 0),
    all(data$avail_hprior_idx <= data$N_avail_random),

    ## Probability of capture model dimensions
    data$N_pcap_fixed >= 0,
    data$N_pcap_prior >= 0,
    data$N_pcap_random >= 0,
    data$N_pcap_hprior >= 0,
    ## Probability of capture model
    nrow(data$pcap_fixed) == data$N_pcap_obs,
    ncol(data$pcap_fixed) == data$N_pcap_fixed,
    length(data$pcap_prior_sd) == data$N_pcap_prior,
    all(data$pcap_prior_sd > 0),
    length(data$pcap_prior_mean) == data$N_pcap_prior,
    length(data$pcap_prior_idx) == data$N_pcap_fixed,
    all(data$pcap_prior_idx > 0),
    all(data$pcap_prior_idx <= data$N_pcap_fixed),
    nrow(data$pcap_random) == data$N_pcap_obs,
    ncol(data$pcap_random) == data$N_pcap_random,
    length(data$pcap_hprior) == data$N_pcap_hprior,
    all(data$pcap_hprior > 0),
    length(data$pcap_hprior_idx) == data$N_pcap_random,
    all(data$pcap_hprior_idx > 0),
    all(data$pcap_hprior_idx <= data$N_pcap_random),

    ## Runsize model dimensions
    data$N_rs_fixed >= 0,
    data$N_rs_prior >= 0,
    data$N_rs_random >= 0,
    data$N_rs_hprior >= 0,
    ## Runsize model
    nrow(data$rs_fixed) == data$N_rs_obs,
    ncol(data$rs_fixed) == data$N_rs_fixed,
    length(data$rs_prior_sd) == data$N_rs_prior,
    all(data$rs_prior_sd > 0),
    length(data$rs_prior_mean) == data$N_rs_prior,
    length(data$rs_prior_idx) == data$N_rs_fixed,
    all(data$rs_prior_idx > 0),
    all(data$rs_prior_idx <= data$N_rs_fixed),
    nrow(data$rs_random) == data$N_rs_obs,
    ncol(data$rs_random) == data$N_rs_random,
    length(data$rs_hprior) == data$N_rs_hprior,
    all(data$rs_hprior > 0),
    length(data$rs_hprior_idx) == data$N_rs_random,
    all(data$rs_hprior_idx > 0),
    all(data$rs_hprior_idx <= data$N_rs_random),

    ## Mix IS flag
    !data$mixis
  )
  invisible(data)
}

fce_covariate_effect <- function(
    use_pars, post, data, mod = c("pcap", "rs", "avail")) {
  mod <- match.arg(mod)
  mlab <- function(s) paste0(mod, s)

  fix_covs <- attr(data[[mlab("_fixed")]], "scale_index")
  rand_covs <- attr(data[[mlab("_random")]], "scale_index")

  ## Convenience to get all fixed, random, or all parameter values
  if (length(use_pars) == 1 && use_pars == "fixed") {
    use_pars <- fix_covs
  } else if (length(use_pars) == 1 && use_pars == "random") {
    use_pars <- rand_covs
  } else if (length(use_pars) == 1 && use_pars == "all") {
    use_pars <- c(fix_covs, rand_covs)
  }

  ## Throw an error if named parameters are missing
  chk <- !use_pars %in% c(fix_covs, rand_covs)
  if (any(chk)) {
    stop(
      paste(use_pars[chk], collapse = ", "),
      " is/are not used in this model. Available covariates are: ",
      paste(c(fix_covs, rand_covs), collapse = ", ")
    )
  }

  fix_idx <- fix_covs %in% use_pars
  rand_idx <- rand_covs %in% use_pars

  if (any(fix_idx)) {
    fix_df <- tibble(
      par = fix_covs[fix_idx],
      post = post[[mlab("_coef_fixed")]][fix_idx]
    )
  } else {
    fix_df <- tibble(
      par = NULL,
      post = NULL
    )
  }
  if (any(rand_idx)) {
    rand_df <- tibble(
      par = rand_covs[rand_idx],
      post = post[[mlab("_coef_random")]][rand_idx]
    )
  } else {
    rand_df <- tibble(
      par = NULL,
      post = NULL
    )
  }
  bind_rows(fix_df, rand_df)
}

fce_get_priors <- function(data, mod = c("pcap", "rs", "avail")) {
  mod <- match.arg(mod)
  mlab <- function(s) paste0(mod, s)

  tibble(
    par = attr(data[[mlab("_fixed")]], "scale_index"),
    idx = data[[mlab("_prior_idx")]]
  ) |>
    distinct() |>
    mutate(
      prior = dist_normal(
        data[[mlab("_prior_mean")]],
        data[[mlab("_prior_sd")]]
      )
    ) |>
    select(-idx)
}

fce_get_hpriors <- function(data, mod = c("pcap", "rs", "avail")) {
  mod <- match.arg(mod)
  mlab <- function(s) paste0(mod, s)

  tibble(
    par = attr(data[[mlab("_random")]], "scale_index"),
    idx = data[[mlab("_hprior_idx")]]
  ) |>
    distinct() |>
    mutate(
      prior = dist_truncated(dist_normal(0, data[[mlab("_hprior")]]), 0)
    ) |>
    select(-idx)
}

fce_random_scale <- function(post, data, mod = c("pcap", "rs", "avail")) {
  mod <- match.arg(mod)
  mlab <- function(s) paste0(mod, s)
  tibble(
    par = attr(data[[mlab("_random")]], "scale_index"),
    idx = data[[mlab("_hprior_idx")]]
  ) |>
    distinct() |>
    mutate(post = post[[mlab("_rand_scale")]][idx]) |>
    select(-idx)
}

fce_partial_effect <- function(
    use_pars, post, data, mod = c("pcap", "rs", "avail")) {
  mod <- match.arg(mod)
  ## Create a function for indexing for each model
  mlab <- function(s) paste0(mod, s)
  fixed <- data[[mlab("_fixed")]]
  random <- data[[mlab("_random")]]

  ## Throw an error if appropriate columns are missing
  chk <- !use_pars %in% c(
    attr(fixed, "scale_index"),
    attr(random, "scale_index")
  )
  if (any(chk)) {
    stop(
      paste(use_pars[chk], collapse = ", "),
      " is/are not used in this model. Available covariates are: ",
      paste(
        c(
          unique(attr(fixed, "scale_index")),
          unique(attr(random, "scale_index"))
        ),
        collapse = ", "
      )
    )
  }

  eff <- matrix(rep(0, nrow(fixed)), ncol = 1)

  keep_fixed <- attr(fixed, "scale_index") %in% use_pars
  if (sum(keep_fixed) > 0) {
    fixed <- fixed[, keep_fixed, drop = FALSE]
    coef_fixed <- post[[mlab("_coef_fixed")]][keep_fixed]
    eff <- eff + fixed %*% coef_fixed
  }

  keep_rand <- attr(random, "scale_index") %in% use_pars
  if (sum(keep_rand) > 0) {
    random <- random[, keep_rand, drop = FALSE]
    coef_rand <- post[[mlab("_coef_random")]][keep_rand]
    eff <- eff + random %*% coef_rand
  }

  drop(eff)
}

fce_linpred <- function(
    post, data,
    pars = c("all", "fixed", "random"),
    mod = c("pcap", "rs", "avail")) {
  pars <- match.arg(pars)
  mod <- match.arg(mod)
  fixpars <- data[[paste0(mod, "_fixed")]] |>
    attr("scale_index") |>
    unique()
  randpars <- data[[paste0(mod, "_random")]] |>
    attr("scale_index") |>
    unique()
  use_pars <- switch(pars,
    "all" = c(fixpars, randpars),
    "fixed" = fixpars,
    "random" = randpars
  )
  fce_partial_effect(use_pars, post, data, mod)
}

fce_lpred_compare_plot <- function(
    pars1, mod1,
    pars2, mod2,
    pred_df, post, data) {
  if (pars1 %in% c("all", "fixed", "random")) {
    lp1 <- fce_linpred(post, data, pars1, mod1)
  } else {
    lp1 <- fce_partial_effect(pars1, post, data, mod1)
  }
  if (pars2 %in% c("all", "fixed", "random")) {
    lp2 <- fce_linpred(post, data, pars2, mod2)
  } else {
    lp2 <- fce_partial_effect(pars2, post, data, mod2)
  }

  pred_df |>
    mutate(
      lp1 = lp1,
      lp2 = lp2
    ) |>
    ggplot(aes(color = month(date, label = TRUE))) +
    geom_path(
      aes(x = median(lp1), y = median(lp2)),
      group = "a"
    ) +
    stat_pointinterval(
      aes(xdist = lp1, y = median(lp2)),
      interval_alpha = 0.25, .width = 0.8
    ) +
    stat_pointinterval(
      aes(x = median(lp1), ydist = lp2),
      interval_alpha = 0.25, .width = 0.8
    ) +
    scale_color_discrete() +
    labs(x = mod1, y = mod2, color = "Month")
}

fce_post_avail <- function(post, data) {
  ## Only calculate what isn't saved as model output
  if (!is.null(post$avail)) {
    avail <- post$avail
  } else {
    if (!is.null(post$avail_uc)) {
      avail_uc <- post$avail_uc
    } else {
      avail_uc <- data$avail_fixed %*% post$avail_coef_fixed +
        data$avail_random %*% post$avail_coef_rand
    }

    avail <- map2(
      head(data$avail_idx, -1),
      tail(data$avail_idx, -1) - 1,
      function(s, e) {
        rfun(softmax)(avail_uc[s:e])
      }
    ) |>
      reduce(c)
  }
  avail
}

fce_post_pcap <- function(post, data, use_op = TRUE) {
  if (!is.null(post$pcap)) {
    pcap <- post$pcap
  } else {
    if (!is.null(post$pcap)) {
      pcap_uc <- post$pcap_uc
    } else {
      pcap_uc <- drop(
        data$pcap_fixed %*% post$pcap_coef_fixed +
          data$pcap_random %*% pcap_coef_rand
      )
    }
    ## Need to use `rfun(plogis)` or just `plogis` based on whether `pcap_uc`
    ## is an `rvar` or not
    if (inherits(pcap_uc, "rvar")) {
      pcap <- rfun(plogis)(pcap_uc)
    } else {
      pcap <- plogis(pcap_uc)
    }
  }
  if (use_op) {
    pcap <- pcap * data$op
  }
  pcap
}

fce_post_rs <- function(post, data) {
  if (!is.null(post$log_lambda)) {
    rs_uc <- post$log_lambda
  } else {
    ## Use drop here to make the result a vector instead of a 1-column matrix
    rs_uc <- drop(
      data$rs_fixed %*% post$rs_coef_fixed +
        data$rs_random %*% post$rs_coef_rand
    )
  }
  ## Don't need to use `rfun` here, dispatches correctly by default
  exp(rs_uc)
}

fce_rarrivals <- function(rs, unm_count) {
  rvar_rng(rtpois, length(rs), rs, a = unm_count - 0.5, b = Inf)
}

fce_collection_efficiency <- function(arr, unm_count) {
  sum(unmarked) / rvar_sum(arr)
}

fce_sample_mixis <- function(data, ...) {
  data$mixis <- 1
  stan(
    here::here("fce-model.stan"),
    data = data,
    pars = "log_lik",
    ...
  )
}
