#' Estimate conditional years of life lost by regression methods
#'
#' `conditional_yll()` implements model-based, observed-group YLL estimators.
#' Unlike [estimand_yll()], this function does not define an intervention or a
#' target population. It estimates remaining life expectancy under observed
#' exposure groups after conditioning on optional baseline covariates.
#'
#' @param data A `data.frame` with one row per individual.
#' @param method Regression method. `"poisson"` uses an age-split Poisson rate
#'   model. `"flexible_parametric"` uses a Royston-Parmar model via
#'   `rstpm2::stpm2()` when `rstpm2` is installed.
#' @param id_var,time_var,event_var,exposure_var Column names in `data`.
#' @param reference_level,exposed_level Values labelling the binary exposure
#'   levels. If `NULL` they are inferred from the data.
#' @param age_at_entry_var Column name for age at study entry.
#' @param confounders_baseline Optional baseline covariates to condition on.
#'   Predictions are averaged over the empirical covariate distribution.
#' @param age_start,age_end,age_interval Starting ages for reported YLL.
#' @param prediction_interval Age-grid spacing used when integrating survival
#'   curves. Smaller values are smoother but slower.
#' @param B Number of bootstrap iterations. Set to `0` to skip intervals.
#' @param seed Integer seed for reproducibility.
#' @param conf_level Confidence level for intervals.
#' @param integration Numerical integration rule.
#' @param rp_df Degrees of freedom for `rstpm2::stpm2()` when
#'   `method = "flexible_parametric"`.
#' @param rp_model One of `"combined"` or `"stratified"`. `"combined"` fits one
#'   model including exposure and covariates; `"stratified"` fits separate
#'   models by exposure group.
#' @param show_progress Show a simple bootstrap progress bar.
#'
#' @return A list with `summary`, `detailed_results`, `meta`, and
#'   `conditional_survival_point`.
#' @export
conditional_yll <- function(
    data,
    method = c("poisson", "flexible_parametric"),
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    confounders_baseline = NULL,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    prediction_interval = 1,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    integration = c("left_rectangle", "trapezoidal"),
    rp_df = 4,
    rp_model = c("combined", "stratified"),
    show_progress = TRUE
) {
  method <- match.arg(method)
  integration <- match.arg(integration)
  rp_model <- match.arg(rp_model)

  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  yll_warn_entry_age_in_confounders(age_at_entry_var, confounders_baseline)

  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level

  if (prediction_interval <= 0) {
    stop("`prediction_interval` must be positive.", call. = FALSE)
  }
  if (age_end <= age_start) {
    stop("`age_end` must be greater than `age_start`.", call. = FALSE)
  }

  set.seed(seed)

  point <- yll_conditional_engine_single(
    data = data,
    method = method,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    exposure_var = exposure_var,
    reference_level = reference_level,
    exposed_level = exposed_level,
    age_at_entry_var = age_at_entry_var,
    confounders_baseline = confounders_baseline,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    prediction_interval = prediction_interval,
    integration = integration,
    rp_df = rp_df,
    rp_model = rp_model
  )
  point_surv <- attr(point, "conditional_survival")

  boot_df <- tibble()
  boot_surv <- NULL
  if (B > 0L) {
    ids <- unique(data[[id_var]])

    one_boot <- function(b) {
      sampled_ids <- yll_bootstrap_ids(ids)
      boot_data <- yll_make_boot_data(data, id_var, sampled_ids)

      out <- tryCatch(
        yll_conditional_engine_single(
          data = boot_data,
          method = method,
          id_var = id_var,
          time_var = time_var,
          event_var = event_var,
          exposure_var = exposure_var,
          reference_level = reference_level,
          exposed_level = exposed_level,
          age_at_entry_var = age_at_entry_var,
          confounders_baseline = confounders_baseline,
          age_start = age_start,
          age_end = age_end,
          age_interval = age_interval,
          prediction_interval = prediction_interval,
          integration = integration,
          rp_df = rp_df,
          rp_model = rp_model
        ),
        error = function(e) NULL
      )

      if (is.null(out)) {
        return(NULL)
      }

      surv <- attr(out, "conditional_survival")
      if (!is.null(surv)) {
        surv[["b"]] <- b
      }
      list(yll = mutate(out, b = b), curves = surv)
    }

    boot_results <- yll_run_conditional_bootstrap(B, one_boot, show_progress)
    boot_results <- boot_results[!vapply(boot_results, is.null, logical(1))]
    if (length(boot_results) > 0) {
      boot_df <- bind_rows(lapply(boot_results, `[[`, "yll"))
      boot_surv <- yll_collect_bootstrap_curves(boot_results)
    }
  }

  yll_build_conditional_result_object(
    point_est = point,
    boot_df = boot_df,
    method = method,
    conf_level = conf_level,
    meta = list(
      method = method,
      B = B,
      seed = seed,
      conf_level = conf_level,
      age_start = age_start,
      age_end = age_end,
      age_interval = age_interval,
      prediction_interval = prediction_interval,
      integration = integration,
      reference_level = reference_level,
      exposed_level = exposed_level,
      confounders_baseline = confounders_baseline,
      rp_df = rp_df,
      rp_model = rp_model
    ),
    conditional_survival_point = point_surv,
    conditional_survival_boot = boot_surv
  )
}

#' @noRd
yll_conditional_engine_single <- function(data,
                                          method,
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
  if (identical(method, "poisson")) {
    return(yll_conditional_poisson_single(
      data = data,
      id_var = id_var,
      time_var = time_var,
      event_var = event_var,
      exposure_var = exposure_var,
      reference_level = reference_level,
      exposed_level = exposed_level,
      age_at_entry_var = age_at_entry_var,
      confounders_baseline = confounders_baseline,
      age_start = age_start,
      age_end = age_end,
      age_interval = age_interval,
      prediction_interval = prediction_interval,
      integration = integration
    ))
  }

  yll_conditional_rp_single(
    data = data,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    exposure_var = exposure_var,
    reference_level = reference_level,
    exposed_level = exposed_level,
    age_at_entry_var = age_at_entry_var,
    confounders_baseline = confounders_baseline,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    prediction_interval = prediction_interval,
    integration = integration,
    rp_df = rp_df,
    rp_model = rp_model
  )
}

# Prepare age-split Lexis data for the Poisson rate model.
# English: This mirrors the epi-demographic Poisson approach in the reference
# code, but keeps user column names configurable for package use.
# 日本語: 参考コードのLexis分割を、package用に変数名可変で実装する。
#' @noRd
yll_conditional_make_lexis <- function(data,
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
yll_conditional_poisson_single <- function(data,
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
  split <- yll_conditional_make_lexis(
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

  age_grid <- yll_conditional_age_grid(age_start, age_end, age_interval, prediction_interval)
  pred_base <- yll_conditional_prediction_base(
    data = data,
    age_grid = age_grid,
    age_var = ".age_mid",
    exposure_var = ".expo",
    reference_level = reference_level,
    exposed_level = exposed_level,
    confounders_baseline = confounders_baseline
  )
  pred_base$lex.dur <- prediction_interval

  pred_ref <- pred_base
  pred_exp <- pred_base
  pred_ref$.expo <- factor(as.character(reference_level), levels = as.character(c(reference_level, exposed_level)))
  pred_exp$.expo <- factor(as.character(exposed_level), levels = as.character(c(reference_level, exposed_level)))

  # English: predict(..., type="response") returns expected deaths over
  # lex.dur. Convert this to an interval hazard probability.
  # 日本語: Poisson予測値は区間内死亡数の期待値なので、1-exp(-rate*duration)で
  # 区間死亡確率へ変換する。
  mu_ref <- stats::predict(fit, newdata = pred_ref, type = "response")
  mu_exp <- stats::predict(fit, newdata = pred_exp, type = "response")
  pred_base$hazard0 <- 1 - exp(-pmax(as.numeric(mu_ref), 0))
  pred_base$hazard1 <- 1 - exp(-pmax(as.numeric(mu_exp), 0))

  curves <- yll_conditional_curves_from_prediction_grid(
    pred_base,
    id_var = id_var,
    age_var = ".age_mid",
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    integration = integration
  )
  attr(curves$results, "conditional_survival") <- curves$survival
  curves$results
}

#' @noRd
yll_conditional_rp_single <- function(data,
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
      "Package 'rstpm2' is required for `method = \"flexible_parametric\"`. ",
      "Install it or use `method = \"poisson\"`.",
      call. = FALSE
    )
  }

  # `rstpm2::stpm2()` rewrites the call to `gsm()` and evaluates it in the
  # caller's frame. The model frame can also need `nsx()` for the baseline
  # spline. Keeping these functions in this local frame lets users call
  # `conditional_yll()` without first attaching rstpm2 with library(rstpm2).
  # `rstpm2::stpm2()` は内部で `gsm()` 呼び出しに置き換えて評価します。
  # ベースラインスプラインのために `nsx()` も必要になるため、ここで両方を
  # ローカル環境に置き、利用者が library(rstpm2) を先に実行していなくても
  # flexible parametric model を推定できるようにします。
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

  age_grid <- yll_conditional_age_grid(age_start, age_end, age_interval, prediction_interval)
  pred_base <- yll_conditional_prediction_base(
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

  surv <- pred_base |>
    group_by(.data[[".age_out"]]) |>
    summarise(
      age_temp = .data[[".age_out"]][1],
      surv0 = mean(.data[["surv0_ind"]], na.rm = TRUE),
      surv1 = mean(.data[["surv1_ind"]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    select(age_temp, surv0, surv1)

  results <- yll_conditional_results_from_survival(
    surv,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    integration = integration
  )
  attr(results, "conditional_survival") <- yll_conditional_survival_stack(surv, age_start, age_end, age_interval)
  results
}

#' @noRd
yll_conditional_age_grid <- function(age_start, age_end, age_interval, prediction_interval) {
  # Prediction may use a fine grid, while the reported YLL can use coarser
  # starting ages. Include both sets so every requested starting age has a
  # survival curve available.
  # 予測用の細かいグリッドと、結果表示用の開始年齢グリッドは一致しないことがある。
  # 両方を含めることで、指定された開始年齢の YLL が NA にならないようにする。
  prediction_grid <- seq(from = age_start, to = age_end, by = prediction_interval)
  reporting_grid <- seq(from = age_start, to = age_end, by = age_interval)
  sort(unique(c(prediction_grid, reporting_grid, age_end)))
}

# Build a marginal-prediction grid.
#
# English: To condition on arbitrary baseline covariates without choosing one
# artificial "typical" patient, we copy the empirical covariate distribution
# for every prediction age and then average predictions over people.
# 日本語: 任意の共変量で条件づけるため、代表値1人ではなく全対象者の共変量分布に
# 予測して平均する。
#' @noRd
yll_conditional_prediction_base <- function(data,
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

#' @noRd
yll_conditional_curves_from_prediction_grid <- function(df,
                                                       id_var,
                                                       age_var,
                                                       age_start,
                                                       age_end,
                                                       age_interval,
                                                       integration) {
  df <- df |>
    arrange(.data[[".prediction_id"]], .data[[age_var]]) |>
    group_by(.data[[".prediction_id"]]) |>
    mutate(
      .surv0_after_interval = cumprod(1 - .data[["hazard0"]]),
      .surv1_after_interval = cumprod(1 - .data[["hazard1"]]),
      surv0_ind = lag(.surv0_after_interval, default = 1),
      surv1_ind = lag(.surv1_after_interval, default = 1)
    ) |>
    ungroup()

  surv <- df |>
    group_by(.data[[age_var]]) |>
    summarise(
      age_temp = .data[[age_var]][1],
      surv0 = mean(.data[["surv0_ind"]], na.rm = TRUE),
      surv1 = mean(.data[["surv1_ind"]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    select(age_temp, surv0, surv1)

  results <- yll_conditional_results_from_survival(
    surv,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    integration = integration
  )

  list(
    results = results,
    survival = yll_conditional_survival_stack(surv, age_start, age_end, age_interval)
  )
}

#' @noRd
yll_conditional_results_from_survival <- function(surv,
                                                  age_start,
                                                  age_end,
                                                  age_interval,
                                                  integration) {
  age_list <- seq(from = age_start, to = age_end, by = age_interval)

  map_dfr(age_list, function(a_start) {
    df_cond <- yll_conditional_survival_from_age(surv, a_start)
    if (is.null(df_cond)) {
      return(tibble(age_start = a_start, yll = NA_real_, le_m0 = NA_real_, le_m1 = NA_real_))
    }

    le_m0 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g0, integration = integration)
    le_m1 <- yll_integrate_le(df_cond$age_temp, df_cond$surv_cond_g1, integration = integration)

    tibble(
      age_start = a_start,
      yll = le_m0 - le_m1,
      le_m0 = le_m0,
      le_m1 = le_m1
    )
  })
}

#' @noRd
yll_conditional_survival_stack <- function(surv, age_start, age_end, age_interval) {
  age_list <- seq(from = age_start, to = age_end, by = age_interval)
  out <- lapply(age_list, function(a_start) {
    df_cond <- yll_conditional_survival_from_age(surv, a_start)
    if (is.null(df_cond)) {
      return(NULL)
    }
    df_cond$starting_age <- a_start
    df_cond
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    return(NULL)
  }
  bind_rows(out)
}

#' @noRd
yll_run_conditional_bootstrap <- function(B, one_boot, show_progress) {
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    on.exit(close(pb), add = TRUE)
    return(lapply(seq_len(B), function(b) {
      setTxtProgressBar(pb, b)
      one_boot(b)
    }))
  }

  lapply(seq_len(B), one_boot)
}

#' @noRd
yll_build_conditional_result_object <- function(point_est,
                                                boot_df,
                                                method,
                                                conf_level,
                                                meta,
                                                conditional_survival_point,
                                                conditional_survival_boot) {
  summary_df <- point_est

  if (nrow(boot_df) > 0) {
    norm <- yll_ci_normal(point_est, boot_df, conf_level) |>
      rename_with(~ paste0(.x, "_norm"), -age_start)
    summary_df <- left_join(summary_df, norm, by = "age_start")
  }

  readable <- yll_make_readable_results(
    point_est = point_est,
    boot_df = boot_df,
    summary_df = arrange(summary_df, age_start),
    method = "normal"
  )

  detailed <- readable$detailed_results |>
    mutate(
      estimate = .data[["yll"]],
      method = method,
      measure = "yll",
      comparison = "conditional"
    )

  summary <- readable$summary |>
    mutate(
      estimate = .data[["yll"]],
      method = method,
      measure = "yll",
      comparison = "conditional"
    )

  out <- list(
    detailed_results = detailed,
    summary = summary,
    meta = meta,
    conditional_survival_point = conditional_survival_point,
    conditional_survival_boot = conditional_survival_boot
  )
  class(out) <- c("conditional_yll", "yll_estimate", class(out))
  out
}
