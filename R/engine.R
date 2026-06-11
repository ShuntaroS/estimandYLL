# Core estimation engine.
#
# Performs one YLL estimation under arbitrary `intervention_reference` /
# `intervention_exposed` specifications and an optional `target_population`.
# This is the function that all user-facing wrappers (`estimate_yll_gformula`,
# `*_intervention`, `*_binary_stochastic`, the ATE/ATT/ATC wrappers, ...)
# eventually call. Doing one bootstrap iteration also calls this function.
#
# The pipeline is:
#
#   1. Validate inputs and discretise follow-up time.
#   2. Fit the discrete-time pooled logistic hazard model on the observed
#      data (full sample — see the note in `intervention.R` about why we do
#      not subset by estimand here).
#   3. Build a counterfactual prediction grid (one row per subject per
#      integer age in [age_start, age_end]).
#   4. Resolve the two intervention specifications into per-row P(A* = exp).
#   5. For each requested starting age a_start, build per-subject survival
#      curves starting from that age (cumprod restarts from 1), average across
#      subjects within the chosen target population to get the conditional
#      marginal curves, and integrate to get residual life expectancies.
#   6. YLL is the difference (LE_ref - LE_exp).
#
# Returns a tibble with one row per starting age, plus the conditional curves
# and the marginal curves (= conditional from the minimum age) as attributes.
#' @noRd
estimate_yll_gformula_engine_single <- function(
    data,
    id_var,
    time_var,
    event_var,
    exposure_var,
    reference_level,
    exposed_level,
    age_at_entry_var,
    age_start,
    age_end,
    age_interval,
    confounders_baseline = NULL,
    intervention_reference,
    intervention_exposed,
    target_population = NULL,
    integration = c("left_rectangle", "trapezoidal")
) {
  integration <- match.arg(integration)

  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )

  cut_points <- yll_make_cut_points(data, time_var = time_var, by = 1L)

  # Person-period long form, then derive `age_temp` (= attained age within
  # the interval) and `expo` (a stable name for the exposure) so the rest of
  # the pipeline doesn't have to know the user's column names.
  df_long <- yll_survsplit(
    data       = data,
    time_var   = time_var,
    event_var  = event_var,
    cut_points = cut_points
  ) |>
    mutate(
      age_temp = .data[[age_at_entry_var]] + tstart,
      expo     = .data[[exposure_var]]
    )

  knot_p <- yll_knot_quantiles(df_long$age_temp)

  fit <- yll_fit_discrete_hazard_glm(
    df_long              = df_long,
    event_var            = event_var,
    confounders_baseline = confounders_baseline,
    knot_p               = knot_p
  )

  # Counterfactual data: build rows for ages in [age_start, age_end].
  cf_base <- yll_expand_counterfactual_data(
    data           = data,
    id_var         = id_var,
    exposure_var   = exposure_var,
    age_temp_start = as.integer(age_start),
    age_temp_end   = as.integer(age_end)
  )

  # Per-subject P(A* = exposed) under each intervention arm.
  p_exposed_reference <- yll_resolve_exposed_probability(
    cf_base,
    intervention_reference,
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  p_exposed_exposed <- yll_resolve_exposed_probability(
    cf_base,
    intervention_exposed,
    reference_level = reference_level,
    exposed_level = exposed_level
  )

  # Predict hazards under "everyone reference" and "everyone exposed" copies.
  cf0 <- cf_base
  cf1 <- cf_base
  cf0$expo <- reference_level
  cf1$expo <- exposed_level
  cf_base$hazard0 <- yll_predict_hazard(fit, cf0)
  cf_base$hazard1 <- yll_predict_hazard(fit, cf1)

  # Apply the target-population mask *before* the per-a_start loop.
  pop_idx <- yll_resolve_population_index(cf_base, target_population)
  cf_pop <- cf_base[pop_idx, , drop = FALSE]
  p_exposed_reference <- p_exposed_reference[pop_idx]
  p_exposed_exposed <- p_exposed_exposed[pop_idx]

  age_list <- seq(from = age_start, to = age_end, by = age_interval)

  # For each requested starting age: subset to ages >= a_start and restart the
  # per-subject cumulative product there, so each subject's curve is
  # conditional on being alive at a_start. Averaging these per-subject
  # conditional curves keeps the estimand identical across starting ages and
  # makes the result independent of the grid lower bound. `cf_pop` is sorted
  # by (id, age), so logical subsetting keeps the intervention-probability
  # vectors aligned with the rows.
  conditional_curves <- map_dfr(age_list, function(a_start) {
    keep <- cf_pop$age_temp >= a_start
    sub_surv <- yll_compute_individual_survival_curves(
      cf_pop[keep, , drop = FALSE],
      id_var = id_var
    )
    s0 <- yll_mean_survival_under_intervention(sub_surv, p_exposed_reference[keep], surv_name = "surv0")
    s1 <- yll_mean_survival_under_intervention(sub_surv, p_exposed_exposed[keep], surv_name = "surv1")
    inner_join(s0, s1, by = "age_temp") |>
      mutate(age_start = a_start, .before = 1)
  })

  # Integrate each conditional curve to get residual life expectancy under
  # each arm; YLL is the difference.
  yll_table <- conditional_curves |>
    group_by(age_start) |>
    summarise(
      le_m0 = yll_integrate_le(age_temp, surv0, integration = integration),
      le_m1 = yll_integrate_le(age_temp, surv1, integration = integration),
      .groups = "drop"
    ) |>
    mutate(yll = le_m0 - le_m1) |>
    select(age_start, yll, le_m0, le_m1)

  # Stash the curves as attributes for the bootstrap collector and the
  # plotting helpers. The "marginal" curves are the a_start = age_start case:
  # the full-grid curves conditional on being alive at age_start.
  marginal_curves <- conditional_curves |>
    filter(age_start == age_list[[1]]) |>
    select(-"age_start")
  attr(yll_table, "marginal_curves") <- as_tibble(marginal_curves)
  attr(yll_table, "conditional_curves") <- as_tibble(conditional_curves)
  yll_table
}

# Backward-compatible wrapper around the engine for the named-estimand
# (ATE / ATT / ATC) interface.
#
# It resolves binary exposure levels, then turns the estimand into a
# `target_population` rule and forwards everything to the engine. The
# `intervention_reference` / `intervention_exposed` arguments are filled in
# with the deterministic "everyone reference" vs "everyone exposed" pair —
# which is exactly the contrast the named estimands require.
#' @noRd
estimate_yll_gformula_single <- function(
    data,
    id_var,
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start,
    age_end,
    age_interval,
    confounders_baseline = NULL,
    estimand = c("ATT", "ATC", "ATE"),
    integration = c("left_rectangle", "trapezoidal")
) {
  estimand <- match.arg(estimand)
  integration <- match.arg(integration)

  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level

  estimate_yll_gformula_engine_single(
    data = data,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    exposure_var = exposure_var,
    reference_level = reference_level,
    exposed_level = exposed_level,
    age_at_entry_var = age_at_entry_var,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    confounders_baseline = confounders_baseline,
    intervention_reference = reference_level,
    intervention_exposed = exposed_level,
    target_population = yll_make_population_rule_from_estimand(
      estimand = estimand,
      reference_level = reference_level,
      exposed_level = exposed_level
    ),
    integration = integration
  )
}
