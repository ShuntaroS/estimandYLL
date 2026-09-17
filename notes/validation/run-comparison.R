# Run from the package root in separate sessions:
# Rscript notes/validation/run-comparison.R legacy
# Rscript notes/validation/run-comparison.R current
# The baseline source is git archive aa9973e, extracted to local/legacy.
version <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(version %in% c("legacy", "current"))
dir.create("local/validation", recursive = TRUE, showWarnings = FALSE)
if (version == "legacy") {
  pkgload::load_all("local/legacy", quiet = TRUE)
  namespace <- asNamespace("estimandYLL")
  # Capture replicate estimates that the old public result did not expose.
  for (name in c("yll_build_result_object", "yll_build_conditional_result_object")) {
    original <- get(name, namespace)
    replacement <- local({
      original_function <- original
      function(...) {
        arguments <- list(...)
        result <- do.call(original_function, arguments)
        result$validation_bootstrap <- arguments$boot_df
        result
      }
    })
    unlockBinding(name, namespace)
    assign(name, replacement, namespace)
    lockBinding(name, namespace)
  }
} else {
  pkgload::load_all(".", quiet = TRUE)
}
data("yll_toy", package = "estimandYLL")
arguments <- list(
  data = yll_toy, id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age", age_start = 40, age_end = 90, age_interval = 10,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  B = 100, seed = 20260917, conf_level = 0.95,
  integration = "left_rectangle", show_progress = FALSE
)
saveRDS(arguments[names(arguments) != "data"], "local/validation/settings.rds")
for (setting in c("all", "exposed", "unexposed", "poisson", "royston_parmar")) {
  destination <- file.path("local/validation", paste0(version, "-", setting, ".rds"))
  if (file.exists(destination)) next
  cat(format(Sys.time()), version, setting, "started\n")
  warnings <- character()
  started <- proc.time()[[3]]
  result <- tryCatch(withCallingHandlers({
    call_arguments <- arguments
    if (setting %in% c("all", "exposed", "unexposed")) {
      call_arguments$target_population <- setting
      call_arguments$use_future <- FALSE
      if (version == "legacy") {
        call_arguments$intervention <- "full_contrast"
        call_arguments$measure <- "yll"
        do.call(estimand_yll, call_arguments)
      } else do.call(estimate_yll, call_arguments)
    } else if (version == "legacy") {
      call_arguments$method <- if (setting == "poisson") "poisson" else "flexible_parametric"
      do.call(conditional_yll, call_arguments)
    } else {
      estimator <- if (setting == "poisson") estimate_yll_poisson else estimate_yll_royston_parmar
      do.call(estimator, call_arguments)
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(error = conditionMessage(e)))
  saveRDS(list(result = result, warnings = warnings,
               elapsed_seconds = proc.time()[[3]] - started,
               session = sessionInfo()), destination)
  cat(format(Sys.time()), version, setting, "saved\n")
}
