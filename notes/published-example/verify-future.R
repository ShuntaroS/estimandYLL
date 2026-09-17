# Use an installed package so each worker loads the same namespace.
# Rscript notes/published-example/verify-future.R
library(estimandYLL)

verify_future <- function() {
  previous_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(previous_plan), add = TRUE)
  data("yll_toy", package = "estimandYLL")

  # A small run checks the backend without repeating the 1,000-replicate fit.
  estimate_small_example <- function() {
    estimate_yll(
      yll_toy[1:1000, ], time_var = "period", event_var = "event",
      exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
      age_at_entry_var = "age", target_population = "unexposed",
      age_start = 50, age_end = 70, age_interval = 5,
      confounders_baseline = "sex", B = 4, seed = 20260917,
      use_future = TRUE, show_progress = FALSE
    )
  }
  first <- estimate_small_example()
  repeated <- estimate_small_example()
  future::plan(future::sequential)
  sequential <- estimate_small_example()

  # Check both numerical outputs and the complete per-replicate survival curves.
  for (result in list(first, repeated, sequential)) {
    stopifnot(result$meta$bootstrap_successful == 4,
              nrow(result$meta$bootstrap_failures) == 0)
  }
  for (element in c("summary", "bootstrap_estimates", "bootstrap_survival_curves")) {
    stopifnot(identical(first[[element]], repeated[[element]]),
              identical(first[[element]], sequential[[element]]))
  }
  message("PASS: two workers, repeated run, and sequential future backend; four successful replicates each.")
}

verify_future()
