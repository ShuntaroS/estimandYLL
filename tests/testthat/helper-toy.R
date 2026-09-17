toy_arguments <- function() {
  list(data = yll_toy[1:1000, ], time_var = "period", event_var = "event",
       exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
       age_at_entry_var = "age", target_population = "all", age_start = 50,
       age_end = 70, age_interval = 10, confounders_baseline = "sex",
       B = 0, seed = 47, use_future = FALSE, show_progress = FALSE)
}

estimate_toy <- function(...) {
  arguments <- utils::modifyList(toy_arguments(), list(...), keep.null = TRUE)
  do.call(estimate_yll, arguments)
}
