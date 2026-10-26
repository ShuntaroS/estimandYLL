# Run from the package root:
# Rscript notes/sandwich-validation/run-validation.R
# Compares variance = "sandwich" with the participant bootstrap on yll_toy and
# checks the sandwich standard error against repeated synthetic datasets.
pkgload::load_all(".", quiet = TRUE)
output <- "notes/sandwich-validation"
workers <- max(1L, parallel::detectCores() - 2L)
data("yll_toy", package = "estimandYLL")

common <- list(time_var = "period", event_var = "event", age_at_entry_var = "age",
               age_start = 40, age_end = 90, age_interval = 10,
               confounders_baseline = c("sex", "education_years", "bmi"))
settings <- list(
  smoking_exposed = list(exposure_var = "smoke_binary", reference_level = "Never",
                         exposed_level = "Current/Ever", target_population = "exposed"),
  smoking_all = list(exposure_var = "smoke_binary", reference_level = "Never",
                     exposed_level = "Current/Ever", target_population = "all"),
  hypertension_unexposed = list(exposure_var = "hypertension", reference_level = "No",
                                exposed_level = "Yes", target_population = "unexposed")
)

# 1. Sandwich versus 1,000 bootstrap replicates --------------------------------
# Forked workers inherit the functions loaded by pkgload::load_all().
previous_plan <- future::plan(future::multicore, workers = workers)
comparison <- dplyr::bind_rows(lapply(names(settings), function(name) {
  arguments <- c(list(data = yll_toy), common, settings[[name]])
  sandwich_seconds <- system.time(
    sandwich <- do.call(estimate_yll, c(arguments, list(variance = "sandwich")))
  )[["elapsed"]]
  bootstrap_seconds <- system.time(
    bootstrap <- do.call(estimate_yll, c(arguments, list(
      B = 1000, seed = 20261026, use_future = TRUE, show_progress = FALSE
    )))
  )[["elapsed"]]
  tibble::tibble(
    setting = name, starting_age = sandwich$summary$starting_age,
    yll = sandwich$summary$yll,
    se_sandwich = sandwich$summary$yll_se,
    se_bootstrap = bootstrap$summary$yll_se,
    ratio = se_sandwich / se_bootstrap,
    seconds_sandwich = sandwich_seconds,
    seconds_bootstrap = bootstrap_seconds,
    bootstrap_workers = workers
  )
}))
future::plan(previous_plan)
write.csv(comparison, file.path(output, "bootstrap-comparison.csv"), row.names = FALSE)

# 2. The published example, against its saved 1,000-replicate bootstrap -------
saved <- readRDS(system.file("extdata", "yll-example.rds", package = "estimandYLL"))
example_seconds <- system.time(example <- estimate_yll(
  data = yll_toy, time_var = "period", event_var = "event",
  exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age", target_population = "unexposed",
  age_start = 40, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  variance = "sandwich"
))[["elapsed"]]
published <- tibble::tibble(
  starting_age = example$summary$starting_age,
  yll = example$summary$yll,
  yll_saved = saved$summary$yll,
  se_sandwich = example$summary$yll_se,
  se_bootstrap = saved$summary$yll_se,
  ratio = se_sandwich / se_bootstrap,
  seconds_sandwich = example_seconds
)
write.csv(published, file.path(output, "published-example-comparison.csv"), row.names = FALSE)

# 3. Mean sandwich SE versus the SD of estimates over repeated datasets -------
# Same data-generating process as data-raw/yll_toy.R.
simulate_toy <- function(n = 5000L, max_follow = 15L) {
  age <- pmin(pmax(round(rnorm(n, 54.5, 7)), 35), 69)
  sex <- factor(rbinom(n, 1, 0.52), levels = c(0, 1), labels = c("Female", "Male"))
  education_years <- pmin(pmax(round(rnorm(
    n, 14.2 - 0.3 * (sex == "Male") - 0.05 * (age - 55), 2.2)), 9), 20)
  bmi <- pmin(pmax(round(rnorm(
    n, 23.5 + 0.06 * (age - 55) + 0.9 * (sex == "Male") - 0.12 * (education_years - 14),
    3.2), 1), 16), 40)
  smoke_binary <- factor(rbinom(n, 1, plogis(
    -0.35 + 0.65 * (sex == "Male") - 0.14 * (education_years - 14) +
      0.03 * (age - 55) + 0.02 * (bmi - 25))),
    levels = c(0, 1), labels = c("Never", "Current/Ever"))
  hypertension <- factor(rbinom(n, 1, plogis(
    -0.8 + 0.09 * (age - 55) + 0.55 * (sex == "Male") + 0.09 * (bmi - 25) -
      0.08 * (education_years - 14) + 0.30 * (smoke_binary == "Current/Ever"))),
    levels = c(0, 1), labels = c("No", "Yes"))
  age_mat <- outer(age, 0:(max_follow - 1L), "+")
  hazard <- plogis(-5.8 + 0.095 * (age_mat - 55) +
    0.42 * (smoke_binary == "Current/Ever") + 0.60 * (hypertension == "Yes") +
    0.16 * (sex == "Male") + 0.025 * (bmi - 25) - 0.05 * (education_years - 14) +
    0.12 * (smoke_binary == "Current/Ever") * (hypertension == "Yes"))
  event_mat <- matrix(runif(n * max_follow), n, max_follow) < hazard
  first <- apply(event_mat, 1, function(x) which(x)[1])
  event <- as.integer(!is.na(first))
  period <- rep(max_follow, n)
  period[event == 1L] <- first[event == 1L] - 1 + runif(sum(event))
  data.frame(id = seq_len(n), period = period, event = event,
             smoke_binary = smoke_binary, hypertension = hypertension, sex = sex,
             education_years = education_years, bmi = bmi, age = age)
}

replicates <- 500L
simulation <- dplyr::bind_rows(parallel::mclapply(seq_len(replicates), function(replicate) {
  set.seed(20261026 + replicate)
  simulated <- simulate_toy()
  dplyr::bind_rows(lapply(names(settings), function(name) {
    result <- do.call(estimate_yll, c(list(data = simulated), common, settings[[name]],
                                      list(variance = "sandwich")))
    tibble::tibble(setting = name, replicate = replicate,
                   starting_age = result$summary$starting_age,
                   yll = result$summary$yll, se = result$summary$yll_se)
  }))
}, mc.cores = workers))
simulation_summary <- simulation |>
  dplyr::filter(starting_age < common$age_end) |>
  dplyr::group_by(setting, starting_age) |>
  dplyr::summarise(replicates = dplyr::n(), mean_yll = mean(yll),
                   empirical_se = sd(yll), mean_sandwich_se = mean(se),
                   se_ratio = mean_sandwich_se / empirical_se, .groups = "drop")
write.csv(simulation_summary, file.path(output, "simulation-summary.csv"), row.names = FALSE)

print(comparison, n = Inf)
print(published, n = Inf)
print(simulation_summary, n = Inf)
