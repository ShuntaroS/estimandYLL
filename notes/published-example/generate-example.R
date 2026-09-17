# Run from the repository root with the current package installed.
# Rscript notes/published-example/generate-example.R
# The saved result is reused across the website, README, and vignettes.
library(estimandYLL)

generate_example <- function() {
  dir.create("local/published-example", recursive = TRUE, showWarnings = FALSE)
  destination <- "local/published-example/result.rds"
  if (file.exists(destination)) {
    message("The completed example already exists: ", destination)
    return(invisible(readRDS(destination)))
  }

  # Use separate R sessions and restore the caller's plan when finished.
  previous_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(previous_plan), add = TRUE)
  data("yll_toy", package = "estimandYLL")
  started <- Sys.time()
  message("Starting 1,000 bootstrap replicates with two workers: ", started)
  result <- estimate_yll(
    data = yll_toy,
    time_var = "period", event_var = "event",
    exposure_var = "hypertension",
    reference_level = "No", exposed_level = "Yes",
    age_at_entry_var = "age", target_population = "unexposed",
    age_start = 40, age_end = 90, age_interval = 5,
    confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
    B = 1000, seed = 20260917, ci_method = "normal",
    use_future = TRUE, show_progress = TRUE
  )
  saveRDS(result, destination, compress = "xz")
  write.csv(result$summary, "notes/published-example/summary.csv", row.names = FALSE)
  capture.output(sessionInfo(), file = "local/published-example/session-info.txt")
  print(result$summary)
  message("Successful replicates: ", result$meta$bootstrap_successful,
          "; failed: ", nrow(result$meta$bootstrap_failures))
  message("Elapsed minutes: ", round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2))
  invisible(result)
}

generate_example()
