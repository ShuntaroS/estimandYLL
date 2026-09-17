# Shared orchestration keeps resampling and result definitions identical across models.
yll_estimate <- function(data, model, id_var, time_var, event_var, exposure_var,
                         reference_level, exposed_level, age_at_entry_var,
                         target_population, age_start, age_end, age_interval,
                         confounders_baseline, B, seed, conf_level, ci_method,
                         integration, show_progress, use_future,
                         prediction_interval = 1, rp_df = 4, rp_model = "combined") {
  data <- yll_check_inputs(
    data, id_var, time_var, event_var, exposure_var, age_at_entry_var,
    confounders_baseline, age_start, age_end, age_interval, B, seed, conf_level,
    show_progress, use_future, model, prediction_interval, rp_df
  )
  levels <- yll_exposure_levels(data[[exposure_var]], reference_level, exposed_level)
  reference_level <- levels$reference
  exposed_level <- levels$exposed
  # Restoring the caller's RNG state also covers errors during fitting.
  withr::local_seed(seed)

  estimate_once <- function(sample_data) {
    arguments <- list(
      data = sample_data, id_var = id_var, time_var = time_var, event_var = event_var,
      exposure_var = exposure_var, reference_level = reference_level,
      exposed_level = exposed_level, age_at_entry_var = age_at_entry_var,
      confounders_baseline = confounders_baseline,
      age_start = age_start, age_end = age_end, age_interval = age_interval,
      integration = integration
    )
    if (model == "pooled_logistic") {
      arguments$target_population <- target_population
      return(do.call(yll_estimate_logistic_once, arguments))
    }
    arguments$prediction_interval <- prediction_interval
    if (model == "poisson") return(do.call(yll_estimate_poisson_once, arguments))
    arguments$rp_df <- rp_df
    arguments$rp_model <- rp_model
    do.call(yll_estimate_royston_parmar_once, arguments)
  }

  point <- estimate_once(data)
  replicates <- yll_run_bootstrap(B, estimate_once, data, id_var, use_future, show_progress)
  target_size <- switch(target_population,
    all = nrow(data),
    exposed = sum(data[[exposure_var]] == exposed_level),
    unexposed = sum(data[[exposure_var]] == reference_level)
  )
  meta <- list(
    model = model, target_population = target_population,
    reference_level = reference_level, exposed_level = exposed_level,
    id_var = id_var, time_var = time_var, event_var = event_var,
    exposure_var = exposure_var, age_at_entry_var = age_at_entry_var,
    confounders_baseline = confounders_baseline,
    n = nrow(data), n_target = target_size,
    age_start = age_start, age_end = age_end, age_interval = age_interval,
    prediction_interval = prediction_interval, integration = integration,
    B = B, seed = seed, conf_level = conf_level, ci_method = ci_method,
    use_future = use_future
  )
  if (model == "royston_parmar") {
    meta$rp_df <- rp_df
    meta$rp_model <- rp_model
  }
  yll_build_result(point, replicates, meta)
}
