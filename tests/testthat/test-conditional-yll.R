test_that("conditional_yll poisson returns a readable result on toy data", {
  skip_if_not_installed("Epi")
  data(yll_toy, envir = environment())

  res <- suppressWarnings(conditional_yll(
    data = yll_toy,
    method = "poisson",
    B = 0,
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
    confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
  ))

  expect_named(res, c(
    "detailed_results", "summary", "meta",
    "conditional_survival_point", "conditional_survival_boot"
  ))
  expect_true(all(c("starting_age", "yll", "estimate", "le_reference", "le_exposed") %in%
                    names(res$summary)))
  expect_true(all(is.finite(res$summary$estimate)))
  expect_equal(res$meta$method, "poisson")
})

test_that("conditional_yll flexible parametric gives an informative dependency error", {
  skip_if(requireNamespace("rstpm2", quietly = TRUE), "rstpm2 is installed")
  data(yll_toy, envir = environment())

  expect_error(
    conditional_yll(
      data = yll_toy,
      method = "flexible_parametric",
      B = 0,
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
      age_interval = 5
    ),
    "rstpm2"
  )
})

test_that("conditional_yll flexible parametric returns estimates when rstpm2 is available", {
  skip_if_not_installed("rstpm2")
  data(yll_toy, envir = environment())

  res <- suppressWarnings(conditional_yll(
    data = yll_toy,
    method = "flexible_parametric",
    B = 0,
    show_progress = FALSE,
    id_var = "id",
    time_var = "period",
    event_var = "event",
    exposure_var = "hypertension",
    reference_level = "No",
    exposed_level = "Yes",
    age_at_entry_var = "age",
    age_start = 60,
    age_end = 70,
    age_interval = 5,
    prediction_interval = 2,
    rp_df = 3,
    confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
  ))

  expect_equal(res$meta$method, "flexible_parametric")
  expect_true(all(is.finite(res$summary$estimate)))
  expect_true(all(c(60, 65, 70) %in% res$summary$starting_age))
})
