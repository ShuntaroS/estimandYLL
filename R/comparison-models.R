# Prepare age-split Lexis data for the Poisson rate model.
# This mirrors the epi-demographic Poisson approach in the reference
# code, but keeps user column names configurable for package use.
#' @noRd
yll_split_age_intervals <- function(data,
                                       time_var,
                                       event_var,
                                       exposure_var,
                                       age_at_entry_var,
                                       reference_level,
                                       exposed_level,
                                       age_end,
                                       split_interval) {
  d <- data |>
    mutate(
      .age_in = .data[[age_at_entry_var]],
      .age_out = .data[[age_at_entry_var]] + .data[[time_var]],
      .event = as.integer(.data[[event_var]]),
      .expo = factor(
        as.character(.data[[exposure_var]]),
        levels = as.character(c(reference_level, exposed_level))
      )
    )

  breaks <- seq(
    from = floor(min(d$.age_in, na.rm = TRUE)),
    to = ceiling(max(age_end, max(d$.age_out, na.rm = TRUE))),
    by = split_interval
  )
  if (tail(breaks, 1) < max(d$.age_out, na.rm = TRUE)) {
    breaks <- c(breaks, ceiling(max(d$.age_out, na.rm = TRUE)))
  }

  lx <- Epi::Lexis(
    entry = list(.lex_age = .age_in),
    exit = list(.lex_age = .age_out),
    entry.status = factor(rep("Alive", nrow(d)), levels = c("Alive", "Dead")),
    exit.status = factor(.event, levels = c(0, 1), labels = c("Alive", "Dead")),
    data = d
  )
  split <- Epi::splitLexis(lx, breaks = breaks)
  split |>
    mutate(
      .death = as.integer(.data[["lex.Xst"]] == "Dead"),
      .age_mid = .data[[".lex_age"]] + .data[["lex.dur"]] / 2
    )
}

#' @noRd
yll_estimate_poisson_once <- function(data,
                                           id_var,
                                           time_var,
                                           event_var,
                                           exposure_var,
                                           reference_level,
                                           exposed_level,
                                           age_at_entry_var,
                                           confounders_baseline,
                                           age_start,
                                           age_end,
                                           age_interval,
                                           prediction_interval,
                                           integration) {
  split <- yll_split_age_intervals(
    data = data,
    time_var = time_var,
    event_var = event_var,
    exposure_var = exposure_var,
    age_at_entry_var = age_at_entry_var,
    reference_level = reference_level,
    exposed_level = exposed_level,
    age_end = age_end,
    split_interval = prediction_interval
  )

  knot_p <- yll_knot_quantiles(split$.age_mid)
  rhs <- yll_rhs_confounders(confounders_baseline)
  form <- stats::as.formula(paste0(
    ".death ~ .expo + Epi::Ns(.age_mid, knots = knot_p) + ",
    ".expo * Epi::Ns(.age_mid, knots = knot_p) + ",
    rhs,
    " + offset(log(lex.dur))"
  ))
  fit <- stats::glm(formula = form, family = stats::poisson(), data = split)

  age_grid <- yll_prediction_ages(age_start, age_end, age_interval, prediction_interval)
  pred_base <- yll_prediction_data(
    data = data,
    age_grid = age_grid,
    age_var = ".age_mid",
    exposure_var = ".expo",
    reference_level = reference_level,
    exposed_level = exposed_level,
    confounders_baseline = confounders_baseline
  )

  # Each prediction row represents the interval from its grid age to the next
  # grid point. `age_grid` can be unevenly spaced (it merges the prediction
  # grid with the reporting grid), so the Poisson offset must use the actual
  # width of each interval, not a constant `prediction_interval`. The width of
  # the final interval never enters the survival curves (the cumulative
  # product is lagged by one), so repeating the last width is harmless.
  interval_widths <- diff(age_grid)
  if (length(interval_widths) == 0L) {
    interval_widths <- prediction_interval
  } else {
    interval_widths <- c(interval_widths, tail(interval_widths, 1))
  }
  pred_base$lex.dur <- rep(interval_widths, times = nrow(data))

  # The model relates the death rate to `.age_mid`, the midpoint of each
  # Lexis split interval. To predict the rate over [a, a + w) we therefore
  # evaluate at a + w/2; `pred_base$.age_mid` itself stays at the interval
  # start because it serves as the age axis downstream.
  level_set <- as.character(c(reference_level, exposed_level))
  pred_ref <- pred_base
  pred_exp <- pred_base
  pred_ref$.expo <- factor(as.character(reference_level), levels = level_set)
  pred_exp$.expo <- factor(as.character(exposed_level), levels = level_set)
  pred_ref$.age_mid <- pred_ref$.age_mid + pred_ref$lex.dur / 2
  pred_exp$.age_mid <- pred_exp$.age_mid + pred_exp$lex.dur / 2

  # predict(..., type="response") returns expected deaths over
  # lex.dur. Convert this to an interval hazard probability.
  mu_ref <- stats::predict(fit, newdata = pred_ref, type = "response")
  mu_exp <- stats::predict(fit, newdata = pred_exp, type = "response")
  pred_base$hazard0 <- 1 - exp(-pmax(as.numeric(mu_ref), 0))
  pred_base$hazard1 <- 1 - exp(-pmax(as.numeric(mu_exp), 0))

  curves <- yll_curves_from_hazards(
    pred_base,
    age_var = ".age_mid",
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    integration = integration
  )
  yll_format_estimate(curves$results, curves$survival)
}

#' @noRd
yll_estimate_royston_parmar_once <- function(data,
                                      id_var,
                                      time_var,
                                      event_var,
                                      exposure_var,
                                      reference_level,
                                      exposed_level,
                                      age_at_entry_var,
                                      confounders_baseline,
                                      age_start,
                                      age_end,
                                      age_interval,
                                      prediction_interval,
                                      integration,
                                      rp_df,
                                      rp_model) {
  if (!requireNamespace("rstpm2", quietly = TRUE)) {
    stop(
      "Package 'rstpm2' is required for estimate_yll_royston_parmar(). ",
      "Install it or use estimate_yll_poisson().",
      call. = FALSE
    )
  }

  # `rstpm2::stpm2()` rewrites the call to `gsm()` and evaluates it in the
  # caller's frame. The model frame can also need `nsx()` for the baseline
  # spline. Keeping these functions in this local frame lets users call
  # `estimate_yll_royston_parmar()` without first attaching rstpm2 with library(rstpm2).
  stpm2 <- get("stpm2", envir = asNamespace("rstpm2"))
  gsm <- get("gsm", envir = asNamespace("rstpm2"))
  nsx <- get("nsx", envir = asNamespace("rstpm2"))
  predict_stpm2 <- methods::selectMethod("predict", signature = "stpm2")

  d <- data |>
    mutate(
      .age_in = .data[[age_at_entry_var]],
      .age_out = .data[[age_at_entry_var]] + .data[[time_var]],
      .event = as.integer(.data[[event_var]]),
      .expo = factor(
        as.character(.data[[exposure_var]]),
        levels = as.character(c(reference_level, exposed_level))
      )
    )

  rhs_cov <- yll_rhs_confounders(confounders_baseline)
  base_rhs <- if (identical(rhs_cov, "1")) ".expo" else paste(".expo", rhs_cov, sep = " + ")
  combined_formula <- stats::as.formula(paste0(
    "survival::Surv(.age_in, .age_out, .event) ~ ",
    base_rhs
  ))
  strat_formula <- stats::as.formula(paste0(
    "survival::Surv(.age_in, .age_out, .event) ~ ",
    rhs_cov
  ))

  if (identical(rp_model, "combined")) {
    fit <- stpm2(combined_formula, data = d, df = rp_df)
  } else {
    fit_ref <- stpm2(strat_formula, data = d[d$.expo == as.character(reference_level), , drop = FALSE], df = rp_df)
    fit_exp <- stpm2(strat_formula, data = d[d$.expo == as.character(exposed_level), , drop = FALSE], df = rp_df)
  }

  age_grid <- yll_prediction_ages(age_start, age_end, age_interval, prediction_interval)
  pred_base <- yll_prediction_data(
    data = data,
    age_grid = age_grid,
    age_var = ".age_out",
    exposure_var = ".expo",
    reference_level = reference_level,
    exposed_level = exposed_level,
    confounders_baseline = confounders_baseline
  )
  pred_base$.age_in <- 0

  pred_ref <- pred_base
  pred_exp <- pred_base
  pred_ref$.expo <- factor(as.character(reference_level), levels = as.character(c(reference_level, exposed_level)))
  pred_exp$.expo <- factor(as.character(exposed_level), levels = as.character(c(reference_level, exposed_level)))

  if (identical(rp_model, "combined")) {
    pred_base$surv0_ind <- as.numeric(predict_stpm2(fit, newdata = pred_ref, type = "surv"))
    pred_base$surv1_ind <- as.numeric(predict_stpm2(fit, newdata = pred_exp, type = "surv"))
  } else {
    pred_base$surv0_ind <- as.numeric(predict_stpm2(fit_ref, newdata = pred_ref, type = "surv"))
    pred_base$surv1_ind <- as.numeric(predict_stpm2(fit_exp, newdata = pred_exp, type = "surv"))
  }

  # Per-subject conditioning, as in the g-formula engine: for each starting
  # age, divide each person's S_i(t) by their own S_i(a_start) before
  # averaging, so the estimand is identical across starting ages.
  age_list <- seq(from = age_start, to = age_end, by = age_interval)

  survival <- map_dfr(age_list, function(a_start) {
    sub <- pred_base[pred_base[[".age_out"]] >= a_start, , drop = FALSE]

    # Normalise at the first grid age >= a_start (avoids exact floating-point
    # matching against a_start itself).
    cond_age <- min(sub[[".age_out"]])
    at_start <- sub[sub[[".age_out"]] == cond_age,
                    c(".prediction_id", "surv0_ind", "surv1_ind"), drop = FALSE]
    names(at_start) <- c(".prediction_id", ".s0_at_start", ".s1_at_start")

    sub |>
      left_join(at_start, by = ".prediction_id") |>
      mutate(
        surv0_cond = .data[["surv0_ind"]] / pmax(.data[[".s0_at_start"]], 1e-15),
        surv1_cond = .data[["surv1_ind"]] / pmax(.data[[".s1_at_start"]], 1e-15)
      ) |>
      group_by(age_temp = .data[[".age_out"]]) |>
      summarise(
        surv0 = mean(.data[["surv0_cond"]]),
        surv1 = mean(.data[["surv1_cond"]]),
        .groups = "drop"
      ) |>
      mutate(age_start = a_start, .before = 1)
  })

  results <- survival |>
    group_by(age_start) |>
    summarise(
      le_m0 = yll_integrate_le(age_temp, surv0, integration = integration),
      le_m1 = yll_integrate_le(age_temp, surv1, integration = integration),
      .groups = "drop"
    ) |>
    mutate(yll = le_m0 - le_m1) |>
    select(age_start, yll, le_m0, le_m1)

  yll_format_estimate(results, survival)
}

#' @noRd
yll_prediction_ages <- function(age_start, age_end, age_interval, prediction_interval) {
  # Prediction may use a fine grid, while the reported YLL can use coarser
  # starting ages. Include both sets so every requested starting age has a
  # survival curve available.
  prediction_grid <- seq(from = age_start, to = age_end, by = prediction_interval)
  reporting_grid <- seq(from = age_start, to = age_end, by = age_interval)
  sort(unique(c(prediction_grid, reporting_grid, age_end)))
}

# Build a marginal-prediction grid.
#
# To condition on arbitrary baseline covariates without choosing one
# artificial "typical" patient, we copy the empirical covariate distribution
# for every prediction age and then average predictions over people.
#' @noRd
yll_prediction_data <- function(data,
                                            age_grid,
                                            age_var,
                                            exposure_var,
                                            reference_level,
                                            exposed_level,
                                            confounders_baseline) {
  n_ages <- length(age_grid)
  base <- data |>
    dplyr::mutate(.prediction_id = dplyr::row_number()) |>
    tidyr::uncount(weights = n_ages, .remove = FALSE) |>
    dplyr::mutate(.prediction_row = dplyr::row_number())

  base[[age_var]] <- rep(age_grid, times = nrow(data))
  base[[exposure_var]] <- factor(as.character(reference_level), levels = as.character(c(reference_level, exposed_level)))

  keep <- unique(c(".prediction_id", ".prediction_row", age_var, exposure_var, confounders_baseline))
  base[, keep, drop = FALSE]
}

# Turn the per-row interval hazards on the prediction grid into per-starting-
# age conditional survival curves and the YLL table.
#
# For each requested starting age the per-subject cumulative product restarts
# at that age (conditioning each subject on being alive there), the resulting
# conditional curves are averaged over subjects, and residual life expectancy
# is the area under each averaged curve. This mirrors the g-formula engine,
# so results do not depend on the grid lower bound.
#' @noRd
yll_curves_from_hazards <- function(df,
                                                        age_var,
                                                        age_start,
                                                        age_end,
                                                        age_interval,
                                                        integration) {
  age_list <- seq(from = age_start, to = age_end, by = age_interval)

  survival <- map_dfr(age_list, function(a_start) {
    df[df[[age_var]] >= a_start, , drop = FALSE] |>
      arrange(.data[[".prediction_id"]], .data[[age_var]]) |>
      group_by(.data[[".prediction_id"]]) |>
      mutate(
        .surv0_after_interval = cumprod(1 - .data[["hazard0"]]),
        .surv1_after_interval = cumprod(1 - .data[["hazard1"]]),
        surv0_ind = lag(.surv0_after_interval, default = 1),
        surv1_ind = lag(.surv1_after_interval, default = 1)
      ) |>
      ungroup() |>
      group_by(age_temp = .data[[age_var]]) |>
      summarise(
        surv0 = mean(.data[["surv0_ind"]]),
        surv1 = mean(.data[["surv1_ind"]]),
        .groups = "drop"
      ) |>
      mutate(age_start = a_start, .before = 1)
  })

  results <- survival |>
    group_by(age_start) |>
    summarise(
      le_m0 = yll_integrate_le(age_temp, surv0, integration = integration),
      le_m1 = yll_integrate_le(age_temp, surv1, integration = integration),
      .groups = "drop"
    ) |>
    mutate(yll = le_m0 - le_m1) |>
    select(age_start, yll, le_m0, le_m1)

  list(results = results, survival = survival)
}

