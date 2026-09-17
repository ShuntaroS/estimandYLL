# Run after both run-comparison.R sessions finish. No model is refitted here.
library(dplyr)
settings <- c("all", "exposed", "unexposed", "poisson", "royston_parmar")
comparison_tables <- list()
curve_tables <- list()
failure_tables <- list()
status_tables <- list()
for (setting in settings) {
  old_run <- readRDS(paste0("local/validation/legacy-", setting, ".rds"))
  new_run <- readRDS(paste0("local/validation/current-", setting, ".rds"))
  old <- old_run$result
  new <- new_run$result
  if (!is.null(old$error) || !is.null(new$error)) {
    stop("A point estimation failed for ", setting, ": ", old$error, " ", new$error)
  }
  old_bootstrap <- old$validation_bootstrap
  old_se <- tapply(old_bootstrap$yll, old_bootstrap$age_start, sd, na.rm = TRUE)
  old_summary <- data.frame(
    starting_age = old$summary$starting_age,
    erl_reference = old$summary$le_reference,
    erl_exposed = old$summary$le_exposed,
    yll = old$summary$yll,
    yll_se = unname(old_se[as.character(old$summary$starting_age)]),
    ci_low = old$summary$ci_low, ci_high = old$summary$ci_high
  )
  methods <- if (setting %in% c("all", "exposed", "unexposed")) c("normal", "percentile") else "normal"
  for (method in methods) {
    old_values <- old_summary
    new_values <- new$summary
    if (method == "percentile") {
      # Reuse the 100 participant resamples; changing interval construction
      # does not require a second bootstrap run.
      for (row in seq_len(nrow(old_values))) {
        age <- old_values$starting_age[row]
        old_values[row, c("ci_low", "ci_high")] <- quantile(
          old_bootstrap$yll[old_bootstrap$age_start == age], c(.025, .975), na.rm = TRUE)
        new_values[row, c("ci_low", "ci_high")] <- as.list(quantile(
          new$bootstrap_estimates$yll[new$bootstrap_estimates$starting_age == age], c(.025, .975)))
      }
    }
    for (metric in setdiff(names(old_values), "starting_age")) {
      comparison_tables[[length(comparison_tables) + 1L]] <- data.frame(
        setting = setting, ci_method = method, starting_age = old_values$starting_age,
        metric = metric, legacy = old_values[[metric]], current = new_values[[metric]],
        difference = new_values[[metric]] - old_values[[metric]]
      )
    }
  }
  old_curves <- old$conditional_survival_point
  curves <- merge(
    data.frame(starting_age = old_curves$age_start, age = old_curves$age_temp,
               reference_legacy = old_curves$surv0, exposed_legacy = old_curves$surv1),
    new$survival_curves, by = c("starting_age", "age")
  )
  curves$setting <- setting
  curves$reference_difference <- curves$survival_reference - curves$reference_legacy
  curves$exposed_difference <- curves$survival_exposed - curves$exposed_legacy
  curve_tables[[setting]] <- curves
  failed_old <- setdiff(seq_len(100), unique(old_bootstrap$b))
  if (length(failed_old)) {
    failure_tables[[paste0(setting, "-old")]] <- data.frame(
      version = "legacy", setting = setting, iteration = failed_old,
      reason = "Not returned by legacy estimator; its internal error handler discarded the reason."
    )
  }
  failed_new <- new$meta$bootstrap_failures
  if (nrow(failed_new)) {
    failure_tables[[paste0(setting, "-new")]] <- data.frame(
      version = "current", setting = setting, failed_new
    )
  }
  status_tables[[setting]] <- data.frame(
    setting = setting, legacy_successful = length(unique(old_bootstrap$b)),
    current_successful = new$meta$bootstrap_successful,
    legacy_warnings = length(old_run$warnings), current_warnings = length(new_run$warnings),
    legacy_seconds = old_run$elapsed_seconds, current_seconds = new_run$elapsed_seconds
  )
}
comparison <- bind_rows(comparison_tables)
curves <- bind_rows(curve_tables)
status <- bind_rows(status_tables)
failures <- bind_rows(failure_tables)
write.csv(comparison, "notes/validation/estimate-comparison.csv", row.names = FALSE)
write.csv(curves, "notes/validation/survival-comparison.csv", row.names = FALSE)
write.csv(status, "notes/validation/run-status.csv", row.names = FALSE)
write.csv(failures, "notes/validation/bootstrap-failures.csv", row.names = FALSE)
summary <- comparison |>
  group_by(setting, ci_method, metric) |>
  summarise(max_absolute_difference = max(abs(difference)), .groups = "drop")
write.csv(summary, "notes/validation/difference-summary.csv", row.names = FALSE)
print(status)
print(summary, n = Inf)
cat("Maximum survival-curve difference:",
    max(abs(c(curves$reference_difference, curves$exposed_difference))), "\n")
