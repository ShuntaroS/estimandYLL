# Fit the hazard model on everyone, then standardize over the requested group.
# Selecting the group before fitting would change the fitted hazard model.
yll_estimate_logistic_once <- function(data, id_var, time_var, event_var,
                                      exposure_var, reference_level, exposed_level,
                                      age_at_entry_var, confounders_baseline,
                                      age_start, age_end, age_interval,
                                      target_population, integration) {
  cut_points <- yll_make_cut_points(data, time_var)
  person_periods <- yll_survsplit(data, time_var, event_var, cut_points)
  person_periods$age_temp <- person_periods[[age_at_entry_var]] + person_periods$tstart
  person_periods$expo <- person_periods[[exposure_var]]
  knots <- yll_knot_quantiles(person_periods$age_temp)
  fit <- yll_fit_discrete_hazard_glm(person_periods, event_var,
                                   confounders_baseline, knots)

  population <- data
  if (target_population == "exposed") {
    population <- data[data[[exposure_var]] == exposed_level, , drop = FALSE]
  } else if (target_population == "unexposed") {
    population <- data[data[[exposure_var]] == reference_level, , drop = FALSE]
  }

  # Only prediction variables are expanded; unrelated input columns are not copied.
  columns <- unique(c(id_var, exposure_var, confounders_baseline))
  predictions <- yll_expand_counterfactual_data(
    population[, columns, drop = FALSE], id_var, exposure_var, age_start, age_end
  )
  predictions$expo <- reference_level
  predictions$hazard0 <- yll_predict_hazard(fit, predictions)
  predictions$expo <- exposed_level
  predictions$hazard1 <- yll_predict_hazard(fit, predictions)

  starting_ages <- seq(age_start, age_end, by = age_interval)
  curves <- lapply(starting_ages, function(starting_age) {
    prediction_window <- predictions[predictions$age_temp >= starting_age, , drop = FALSE]
    # Restart each person's survival at one before taking the population mean.
    individual_curves <- yll_compute_individual_survival_curves(prediction_window, id_var)
    population_curves <- individual_curves |>
      group_by(age_temp) |>
      summarise(surv0 = mean(surv0), surv1 = mean(surv1), .groups = "drop")
    population_curves$age_start <- starting_age
    population_curves
  })
  curves <- bind_rows(curves)
  estimates <- curves |>
    group_by(age_start) |>
    summarise(
      le_m0 = yll_integrate_le(age_temp, surv0, integration),
      le_m1 = yll_integrate_le(age_temp, surv1, integration),
      .groups = "drop"
    )
  estimates$yll <- estimates$le_m0 - estimates$le_m1
  yll_format_estimate(estimates, curves)
}

# Give all model engines the same descriptive output names.
yll_format_estimate <- function(estimates, curves) {
  estimates <- tibble::tibble(
    starting_age = estimates$age_start,
    erl_reference = estimates$le_m0,
    erl_exposed = estimates$le_m1,
    yll = estimates$yll
  )
  curves <- tibble::tibble(
    starting_age = curves$age_start,
    age = curves$age_temp,
    survival_reference = curves$surv0,
    survival_exposed = curves$surv1
  )
  if (any(!is.finite(as.matrix(estimates))) || any(!is.finite(as.matrix(curves)))) {
    stop("The fitted model produced non-finite estimates or survival probabilities.", call. = FALSE)
  }
  if (any(curves$survival_reference < -1e-8 | curves$survival_reference > 1 + 1e-8) ||
      any(curves$survival_exposed < -1e-8 | curves$survival_exposed > 1 + 1e-8)) {
    warning("The fitted model produced survival probabilities outside zero and one. ",
            "These estimates are unreliable; inspect the model and data.", call. = FALSE)
  }
  list(estimates = estimates, curves = curves)
}
