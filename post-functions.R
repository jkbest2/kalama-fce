### MCMC Diagnostics ----------------------------------------------------------
check_rhat <- function(post) {
  post |>
    keep_at(c(
      "avail_coef_fixed", "avail_coef_random", "avail_rand_scale",
      "pcap_coef_fixed", "pcap_coef_random", "pcap_rand_scale",
      "rs_coef_fixed", "rs_coef_random", "rs_rand_scale",
      "lp__"
    )) |>
    map(posterior::rhat) |>
    map(\(.x) c(rhat_min = min(.x), rhat_max = max(.x)))
}

check_ess_bulk <- function(post) {
  post |>
    keep_at(c(
      "avail_coef_fixed", "avail_coef_random", "avail_rand_scale",
      "pcap_coef_fixed", "pcap_coef_random", "pcap_rand_scale",
      "rs_coef_fixed", "rs_coef_random", "rs_rand_scale",
      "lp__"
    )) |>
    map(posterior::ess_bulk) |>
    map_dbl(min)
}

check_ess_tail <- function(post) {
  post |>
    keep_at(c(
      "avail_coef_fixed", "avail_coef_random", "avail_rand_scale",
      "pcap_coef_fixed", "pcap_coef_random", "pcap_rand_scale",
      "rs_coef_fixed", "rs_coef_random", "rs_rand_scale",
      "lp__"
    )) |>
    map(posterior::ess_tail) |>
    map_dbl(min)
}

convergence_diagnostics <- function(fit) {
  post <- as_draws_rvars(fit)
  div <- get_num_divergent(fit)
  rhat_range_df <- check_rhat(post) |>
    enframe(name = "par", value = "rhat") |>
    unnest_wider(rhat)
  ess_bulk_df <- check_ess_bulk(post) |>
    enframe(name = "par", value = "ess_bulk")
  ess_tail_df <- check_ess_tail(post) |>
    enframe(name = "par", value = "ess_tail")
  diag_df <- tibble(
    par = c(
      "avail_coef_fixed", "avail_coef_random", "avail_rand_scale",
      "pcap_coef_fixed", "pcap_coef_random", "pcap_rand_scale",
      "rs_coef_fixed", "rs_coef_random", "rs_rand_scale",
      "lp__"
    )
  ) |>
    left_join(rhat_range_df, by = join_by(par)) |>
    left_join(ess_bulk_df, by = join_by(par)) |>
    left_join(ess_tail_df, by = join_by(par))
  # structure(
  new_tibble(
    diag_df,
    class = "mcmc_diagnostics",
    n_divergent = div,
    name = fit@model_name
  )
}

n_divergent <- function(diag_df) {
  attr(diag_df, "n_divergent")
}

## Posterior plots ------------------------------------------------------------
## Plot availabilitiy over time
avail_plot <- function(post, prep, data, name, base_dir) {
  release_df <- prep$release_df
  recap_df <- prep$recap_df
  unmarked <- prep$unmarked
  pred_df <- prep$pred_df

  avail_vec <- fce_post_avail(post, data)

  plt <- attr(data, "avail_pform")$data |>
    mutate(
      avail_prob = avail_vec,
      month = factor(month(release_date, label = TRUE, abbr = FALSE))
    ) |>
    curve_interval(avail_prob, .width = 0.8) |>
    ggplot(aes(
      x = travel_time,
      y = avail_prob,
      group = release_site,
      color = release_site,
      fill = release_site
    )) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.25, color = NA) +
    geom_line() +
    scale_y_continuous(
      breaks = waiver(),
      minor_breaks = seq(0, 1, 0.05),
      labels = scales::percent
    ) +
    scale_color_discrete() +
    labs(
      x = "Travel Time (days)",
      y = "Probability of arrival",
      color = "Release Site",
      fill = "Release Site"
    ) +
    theme_minimal()

  ggsave(
    here::here(base_dir, paste0(name, "-avail.png")),
    plt,
    width = 12, height = 8
  )

  plt
}

## Get observed recaptures for each release
recap_obs <- function(prep, recap_end = NULL) {
  release_df <- prep$release_df
  recap_df <- prep$recap_df
  unmarked <- prep$unmarked
  pred_df <- prep$pred_df

  trap_end <- recap_end %||% max(pred_df$date)

  recap_df |>
    summarize(
      recap_end = max(recapture_date, na.rm = TRUE),
      recap_end = pmin(recap_end, trap_end),
      recap_count = sum(count),
      .by = c(release_date, species)
    ) |>
    left_join(release_df, by = join_by(release_date, species)) |>
    mutate(recap_rate = recap_count / count)
}

## Plot probability of recapture
pcap_plot <- function(
    post, prep, data,
    name, base_dir,
    .width = 0.95,
    plot_noop = FALSE, noop_bg = FALSE) {
  release_df <- prep$release_df
  recap_df <- prep$recap_df
  unmarked <- prep$unmarked
  pred_df <- prep$pred_df

  pred_plot_df <- pred_df |>
    mutate(
      pcap = fce_post_pcap(post, data),
      unmarked = unmarked$count,
      year = factor(year(date))
    )

  if (!plot_noop) {
    pred_plot_df <- pred_plot_df |>
      mutate(pcap = if_else(op, pcap, NA))
  }

  ## Summarize the recapture information so the observed data can be compared to
  ## the fitted capture probabilities.
  recap_obs <- recap_obs(prep)

  plt <- pred_plot_df |>
    point_interval(pcap, .width = .width) |>
    ggplot(aes(x = date, y = pcap))
  if (noop_bg) {
    op_df <- pred_df |>
      select(date, op) |>
      filter(!op)
    plt <- plt +
      geom_rect(
        aes(x = date),
        ymin = 0, ymax = Inf, width = 1,
        color = NA, fill = "red", alpha = 0.1,
        data = op_df, inherit.aes = FALSE
      )
  }
  plt <- plt +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), color = NA, alpha = 0.4) +
    geom_segment(
      data = recap_obs,
      aes(
        x = release_date, xend = recap_end,
        y = recap_rate, yend = recap_rate,
        color = release_site
      ), alpha = 0.8
    ) +
    geom_point(
      data = recap_obs,
      aes(
        x = release_date, y = recap_rate,
        size = count, color = release_site
        # shape = release_site
      ), alpha = 0.8
    ) +
    geom_line() +
    facet_wrap(~year, ncol = 1) +
    scale_x_date(date_breaks = "1 month", date_labels = "%B") +
    scale_y_continuous(
      labels = scales::percent,
      expand = expansion(c(0, 0.02))
    ) +
    labs(
      x = "Date",
      y = "Probability of capture",
      color = "Release_site"
    ) +
    theme_minimal()
  ggsave(
    here::here(base_dir, paste0(name, "-pcap.png")),
    plt,
    width = 12, height = 8
  )
  plt
}

## Plot predicted run size
runsize_plot <- function(
    post, prep, data,
    name, base_dir,
    .width = 0.95, noop_bg = FALSE) {
  release_df <- prep$release_df
  recap_df <- prep$recap_df
  unmarked <- prep$unmarked
  pred_df <- prep$pred_df

  pred_plot_df <- pred_df |>
    mutate(
      rs = fce_post_rs(post, data),
      unmarked = unmarked$count,
      arr = fce_rarrivals(rs, unmarked)
    )

  plt <- pred_plot_df |>
    point_interval(arr, .width = 0.95) |>
    ggplot(aes(x = date))
  if (noop_bg) {
    op_df <- pred_df |>
      select(date, op) |>
      filter(!op)
    plt <- plt + geom_rect(
      aes(x = date),
      ymin = 0, ymax = Inf, width = 1,
      color = NA, fill = "red", alpha = 0.1,
      data = op_df, inherit.aes = FALSE
    )
  }
  plt <- plt +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
    geom_line(aes(y = arr)) +
    scale_x_date(date_breaks = "1 month", date_labels = "%B") +
    scale_y_continuous(
      breaks = seq(0, 10e3, 100),
      minor_breaks = seq(50, 10e3, 100),
      expand = expansion(c(0, 0.02))
    ) +
    labs(
      x = "Date",
      y = "Number of Smolts"
    ) +
    theme_minimal()
  ggsave(
    here::here(base_dir, paste0(name, "-runsize.png")),
    plt,
    width = 12, height = 8
  )
  plt
}

fce_abundance <- function(post, prep, data, ...) {
  release_df <- prep$release_df
  recap_df <- prep$recap_df
  unmarked <- prep$unmarked
  pred_df <- prep$pred_df

  pred_df |>
    mutate(
      rs = fce_post_rs(post, data),
      unmarked = unmarked$count,
      arr = fce_rarrivals(rs, unmarked)
    ) |>
    group_by(...) |>
    summarize(
      abundance = rvar_sum(arr)
    )
}

## Plot total abundance
fce_abund_plot <- function(
    post, prep, data,
    name, base_dir,
    .width = c(0.66, 0.95)) {
  pred_plot_df <- fce_abundance(post, prep, data)

  abund_plt <- ggplot(pred_plot_df, aes(xdist = abundance)) +
    stat_halfeye() +
    scale_x_continuous(
      labels = scales::comma,
      breaks = seq(0, 100e3, 1000),
      minor_breaks = seq(0, 100e3, 250)
    ) +
    labs(
      x = "Seasonal abundance",
      y = ""
    ) +
    theme_minimal()

  ggsave(
    here::here(base_dir, paste0(name, "-seasonal.png")),
    abund_plt,
    width = 12, height = 8
  )

  abund_plt
}

## Plot covariate effects
fce_coveff_plot <- function(
    use_pars = "fixed", mod = "pcap",
    post, prep, data, incl_prior = TRUE) {
  par_prior <- fce_get_priors(data, mod)
  par_post <- fce_covariate_effect(use_pars, post, data, mod)

  xlims <- range(quantile(par_post$post, c(0.001, 0.999)))

  if (incl_prior) {
    df <- left_join(par_post, par_prior, by = join_by(par))
    plt <- ggplot(df) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      stat_slab(
        aes(xdist = prior, y = par),
        alpha = 0.25
      ) +
      stat_halfeye(
        aes(xdist = post, y = par, color = par, fill = par),
        alpha = 0.75
      ) +
      scale_thickness_shared() +
      labs(x = "Value", y = "Parameter") +
      guides(color = "none", fill = "none") +
      coord_cartesian(xlim = xlims)
  } else {
    plt <- ggplot(par_post) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      stat_halfeye(
        aes(xdist = post, y = par, color = par, fill = par),
        alpha = 0.75
      ) +
      scale_thickness_shared() +
      labs(x = "Value", y = "Parameter") +
      guides(color = "none", fill = "none") +
      coord_cartesian(xlim = xlims)
  }
  plt
}

fce_create_eff_df <- function(
    use_pars, mod = c("pcap", "rs", "avail"),
    post, prep, data) {
  mod <- match.arg(mod)
  eff_df <- prep$pred_df
  for (par in use_pars) {
    eff_df <- eff_df |>
      mutate("{par}_eff" := fce_partial_effect(
        use_pars = par, post, data, mod
      ))
  }
  eff_df
}

## Plot fitted random effect scale parameters
fce_randscale_plot <- function(
    post, prep, data,
    mod = c("pcap", "rs", "avail"), incl_prior = TRUE) {
  mod <- match.arg(mod)
  df <- full_join(
    fce_get_hpriors(data, mod = mod),
    fce_random_scale(post, data, mod),
    by = join_by(par)
  )

  ## If looking at the availabiltiy model, recode the parameters to increasing
  ## upstream release site distance
  if (mod == "avail" &&
    all(c("tt_lower", "tt_yellow", "tt_jacks") %in% df$par)) {
    df <- df |>
      mutate(par = factor(par, levels = c("tt_lower", "tt_yellow", "tt_jacks")))
  }

  xlims <- c(0, max(quantile(df$post, c(0.001, 0.999))))

  if (incl_prior) {
    plt <- ggplot(df) +
      stat_slab(
        aes(xdist = prior, y = par),
        alpha = 0.25
      ) +
      stat_halfeye(
        aes(xdist = post, y = par, color = par, fill = par),
        alpha = 0.75
      ) +
      scale_thickness_shared() +
      labs(x = "Value", y = "Parameter") +
      guides(color = "none", fill = "none") +
      coord_cartesian(xlim = xlims)
  }
  plt
}
