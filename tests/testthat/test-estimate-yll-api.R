test_that("estimand_yll reproduces legacy full-contrast estimands", {
  skip_if_not_installed("Epi")
  data(yll_toy, envir = environment())

  common_args <- list(
    data = yll_toy,
    B = 0,
    seed = 1,
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "event",
    exposure_var = "hypertension",
    reference_level = "No",
    exposed_level = "Yes",
    age_at_entry_var = "age",
    age_start = 60,
    age_end = 65,
    age_interval = 5,
    confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
    use_future = FALSE
  )

  legacy_att <- suppressWarnings(do.call(estimate_yll_gformula_att, common_args))
  new_att <- suppressWarnings(do.call(estimand_yll, c(
    common_args,
    list(target_population = "exposed", intervention = "full_contrast", measure = "yll")
  )))

  legacy_atc <- suppressWarnings(do.call(estimate_yll_gformula_atc, common_args))
  new_atc <- suppressWarnings(do.call(estimand_yll, c(
    common_args,
    list(target_population = "unexposed", intervention = "full_contrast", measure = "yll")
  )))

  expect_equal(new_att$summary$yll, legacy_att$summary$yll, tolerance = 1e-10)
  expect_equal(new_atc$summary$yll, legacy_atc$summary$yll, tolerance = 1e-10)
})

test_that("partial_change reproduces the legacy stochastic wrapper with explicit target population", {
  skip_if_not_installed("Epi")
  data(yll_toy, envir = environment())

  common_args <- list(
    data = yll_toy,
    B = 0,
    seed = 1,
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "event",
    exposure_var = "hypertension",
    reference_level = "No",
    exposed_level = "Yes",
    age_at_entry_var = "age",
    age_start = 60,
    age_end = 65,
    age_interval = 5,
    confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
    use_future = FALSE
  )

  target_unexposed <- function(data) data$expo_original == "No"

  legacy <- suppressWarnings(do.call(estimate_yll_gformula_binary_stochastic_vs_natural, c(
    common_args,
    list(
      prob_exposed_if_unexposed = 0.2,
      prob_exposed_if_exposed = 1,
      target_population = target_unexposed
    )
  )))

  new <- suppressWarnings(do.call(estimand_yll, c(
    common_args,
    list(
      target_population = "unexposed",
      intervention = "partial_change",
      change_from = "reference",
      change_to = "exposed",
      change_probability = 0.2,
      measure = "yll"
    )
  )))

  expect_equal(new$summary$yll, legacy$summary$yll, tolerance = 1e-10)
  expect_equal(new$summary$estimate, legacy$summary$yll, tolerance = 1e-10)
})

test_that("life_year_change makes ATC-style harmful changes negative", {
  fake <- list(
    summary = tibble::tibble(starting_age = 40, yll = 2, ci_low = 1, ci_high = 3),
    detailed_results = tibble::tibble(
      starting_age = 40,
      yll = 2,
      le_reference = 40,
      le_exposed = 38
    ),
    meta = list(conf_level = 0.95)
  )

  out <- yll_apply_measure(
    res = fake,
    measure = "life_year_change",
    target_population_label = "unexposed",
    intervention_label = "full_contrast",
    intervention_direction = "full_contrast"
  )

  expect_equal(out$summary$estimate, -2)
  expect_equal(out$summary$ci_low, -3)
  expect_equal(out$summary$ci_high, -1)
  expect_equal(out$detailed_results$le_before, 40)
  expect_equal(out$detailed_results$le_after, 38)
})

test_that("estimate_yll remains a compatibility alias", {
  skip_if_not_installed("Epi")
  data(yll_toy, envir = environment())

  args <- list(
    data = yll_toy,
    B = 0,
    seed = 1,
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "event",
    exposure_var = "hypertension",
    reference_level = "No",
    exposed_level = "Yes",
    age_at_entry_var = "age",
    age_start = 60,
    age_end = 65,
    age_interval = 5,
    confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
    target_population = "all",
    intervention = "full_contrast",
    measure = "yll",
    use_future = FALSE
  )

  expect_equal(
    suppressWarnings(do.call(estimate_yll, args))$summary$estimate,
    suppressWarnings(do.call(estimand_yll, args))$summary$estimate
  )
})
