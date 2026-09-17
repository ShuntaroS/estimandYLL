test_that("all target populations use the same fitted model and sign", {
  all <- estimate_toy()
  exposed <- estimate_toy(target_population = "exposed")
  unexposed <- estimate_toy(target_population = "unexposed")
  prevalence <- mean(toy_arguments()$data$hypertension == "Yes")
  expect_equal(all$summary$yll,
               prevalence * exposed$summary$yll + (1 - prevalence) * unexposed$summary$yll,
               tolerance = 1e-10)
  for (result in list(all, exposed, unexposed)) {
    expect_named(result, c("summary", "bootstrap_estimates", "survival_curves",
                           "bootstrap_survival_curves", "meta"))
    expect_named(result$summary, c("starting_age", "erl_reference", "erl_exposed",
                                   "yll", "yll_se", "ci_low", "ci_high"))
    expect_equal(result$summary$yll, result$summary$erl_reference - result$summary$erl_exposed)
    expect_true(all(is.na(result$summary$ci_low)))
    expect_equal(nrow(result$bootstrap_estimates), 0)
    curves <- result$survival_curves
    at_start <- curves$age == curves$starting_age
    expect_true(all(curves$survival_reference[at_start] == 1))
    expect_true(all(curves$survival_exposed[at_start] == 1))
    end <- result$summary[result$summary$starting_age == 70, ]
    expect_equal(end$yll, 0)
    expect_equal(end$erl_reference, 0)
    expect_equal(end$erl_exposed, 0)
  }
  expect_equal(estimate_toy(age_start = 70)$summary$yll, 0)
  expect_equal(estimate_toy(age_start = 60)$summary$yll, all$summary$yll[2:3])
})

test_that("reversing the exposure labels reverses the all-population contrast", {
  forward <- estimate_toy()
  reverse <- estimate_toy(reference_level = "Yes", exposed_level = "No")
  expect_equal(forward$summary$yll, -reverse$summary$yll, tolerance = 1e-10)
})

test_that("bootstrap intervals use the selected method and restore random state", {
  set.seed(83)
  state <- .Random.seed
  normal <- estimate_toy(B = 4)
  expect_identical(.Random.seed, state)
  percentile <- estimate_toy(B = 4, ci_method = "percentile")
  expect_identical(normal$bootstrap_estimates, percentile$bootstrap_estimates)
  expect_identical(normal$bootstrap_estimates, estimate_toy(B = 4)$bootstrap_estimates)
  first_age <- subset(normal$bootstrap_estimates, starting_age == 50)$yll
  expect_equal(normal$summary$ci_low[1], normal$summary$yll[1] - qnorm(.975) * sd(first_age))
  expect_equal(percentile$summary$ci_low[1], unname(quantile(first_age, .025)))
  expect_equal(percentile$summary$yll_se[1], sd(first_age))
  expect_equal(normal$meta$bootstrap_successful, 4)
  expect_equal(nrow(normal$meta$bootstrap_failures), 0)
})

test_that("resampled copies have unique identifiers", {
  data <- toy_arguments()$data
  set.seed(4)
  sampled <- yll_resample_people(data, "id")
  expect_equal(nrow(sampled), nrow(data))
  expect_equal(anyDuplicated(sampled$id), 0L)
  expect_true(any(grepl("_rep2", sampled$id)))
  data$id <- factor(data$id)
  sampled_factor_ids <- yll_resample_people(data, "id")
  expect_equal(anyDuplicated(sampled_factor_ids$id), 0L)
  expect_false(anyNA(sampled_factor_ids$id))
})

test_that("failed bootstrap replicates are reported and retained in metadata", {
  point <- list(estimates = tibble::tibble(starting_age = 50, erl_reference = 10,
                                         erl_exposed = 8, yll = 2),
                curves = tibble::tibble(starting_age = 50, age = 50,
                                       survival_reference = 1, survival_exposed = 1))
  good <- point
  good$estimates$iteration <- 1L
  good$curves$iteration <- 1L
  failure <- list(failure = tibble::tibble(iteration = 2L, reason = "Example fit failure"))
  expect_warning(result <- yll_build_result(point, list(good, failure, good),
                                            list(B = 3, age_end = 70, conf_level = .95, ci_method = "normal")),
                 "1 of 3 bootstrap replicates failed")
  expect_equal(result$meta$bootstrap_failures$reason, "Example fit failure")
  expect_equal(result$meta$bootstrap_successful, 2)
})

test_that("finite but impossible bootstrap lifetimes are not silently accepted", {
  point <- list(estimates = tibble::tibble(starting_age = 50, erl_reference = 10,
                                         erl_exposed = 8, yll = 2),
                curves = tibble::tibble(starting_age = 50, age = 50,
                                       survival_reference = 1, survival_exposed = 1))
  bad <- point
  bad$estimates$erl_reference <- 1000
  bad$estimates$yll <- 992
  bad$estimates$iteration <- 1L
  bad$curves$iteration <- 1L
  expect_warning(result <- yll_build_result(point, list(bad, bad),
    list(B = 2, age_end = 70, conf_level = .95, ci_method = "normal")), "numerically unstable")
  expect_equal(nrow(result$meta$bootstrap_unreliable), 2)
  expect_equal(result$bootstrap_estimates$yll, c(992, 992))
})

test_that("invalid data and removed arguments fail explicitly", {
  arguments <- toy_arguments()
  arguments$target_population <- NULL
  expect_error(do.call(estimate_yll, arguments), "target_population|arg")
  expect_error(estimate_toy(intervention = "partial_change"), "unused argument")
  expect_error(estimate_toy(measure = "life_year_change"), "unused argument")
  expect_error(estimate_toy(B = 1), "at least 2")
  expect_error(estimate_toy(age_start = 50.5), "integer")
  expect_error(estimate_toy(age_interval = 0), "age_interval")
  expect_error(estimate_toy(age_end = 40), "age_end")
  data <- toy_arguments()$data
  data$id[2] <- data$id[1]
  expect_error(estimate_toy(data = data), "unique person")
  data <- toy_arguments()$data
  data$sex[1] <- NA
  expect_error(estimate_toy(data = data), "missing values")
  data <- toy_arguments()$data
  data$hypertension <- "No"
  expect_error(estimate_toy(data = data), "two observed exposure levels")
})

test_that("future bootstrap is reproducible without changing the user's plan", {
  skip_if_not_installed("future")
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  result1 <- estimate_toy(B = 3, use_future = TRUE)
  result2 <- estimate_toy(B = 3, use_future = TRUE)
  expect_identical(result1$bootstrap_estimates, result2$bootstrap_estimates)
  expect_s3_class(future::plan(), "sequential")
})
