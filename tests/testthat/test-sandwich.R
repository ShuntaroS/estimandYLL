# Sandwich variance: compare the chunked covariance with a direct
# influence-function calculation that uses predict() and numerical derivatives.

direct_sandwich_se <- function(arguments) {
  data <- arguments$data
  cut_points <- yll_make_cut_points(data, arguments$time_var)
  person_periods <- yll_survsplit(data, arguments$time_var, arguments$event_var, cut_points)
  person_periods$age_temp <- person_periods[[arguments$age_at_entry_var]] + person_periods$tstart
  person_periods$expo <- person_periods[[arguments$exposure_var]]
  fit <- yll_fit_discrete_hazard_glm(person_periods, arguments$event_var,
                                     arguments$confounders_baseline,
                                     yll_knot_quantiles(person_periods$age_temp))
  beta <- coef(fit)
  target <- switch(arguments$target_population,
    all = seq_len(nrow(data)),
    exposed = which(data[[arguments$exposure_var]] == arguments$exposed_level),
    unexposed = which(data[[arguments$exposure_var]] == arguments$reference_level)
  )
  starting_ages <- seq(arguments$age_start, arguments$age_end, by = arguments$age_interval)
  ages <- arguments$age_start:arguments$age_end
  grid <- data[rep(target, each = length(ages)), , drop = FALSE]
  grid$age_temp <- rep(ages, times = length(target))

  # Person-level YLL at each starting age for a given coefficient vector.
  person_yll <- function(coefficients) {
    fit$coefficients <- coefficients
    erl <- lapply(list(arguments$reference_level, arguments$exposed_level), function(level) {
      grid$expo <- level
      hazard <- matrix(predict(fit, grid, type = "response"), ncol = length(ages), byrow = TRUE)
      sapply(starting_ages, function(starting_age) {
        window <- ages >= starting_age
        if (sum(window) < 2L) return(rep(0, nrow(hazard)))
        surv <- t(apply(1 - hazard[, window, drop = FALSE], 1, cumprod))
        rowSums(cbind(1, surv[, -ncol(surv), drop = FALSE])[, -sum(window), drop = FALSE])
      })
    })
    matrix(erl[[1]] - erl[[2]], nrow = length(target))
  }
  values <- person_yll(beta)
  gradient <- sapply(seq_along(beta), function(k) {
    step <- 1e-5 * max(1, abs(beta[k]))
    up <- beta
    down <- beta
    up[k] <- up[k] + step
    down[k] <- down[k] - step
    (colMeans(person_yll(up)) - colMeans(person_yll(down))) / (2 * step)
  })
  score <- rowsum(model.matrix(fit) * (fit$y - fitted(fit)), person_periods$id, reorder = FALSE)
  influence <- score[as.character(data$id), , drop = FALSE] %*% vcov(fit) %*% t(gradient)
  influence[target, ] <- influence[target, ] +
    sweep(values, 2, colMeans(values)) / length(target)
  sqrt(colSums(influence^2))
}

test_that("sandwich standard errors match a direct influence-function calculation", {
  for (population in c("unexposed", "all")) {
    arguments <- utils::modifyList(toy_arguments(), list(
      target_population = population, age_start = 50, age_end = 70, age_interval = 5
    ))
    result <- do.call(estimate_yll, c(arguments[setdiff(names(arguments), "B")],
                                      list(variance = "sandwich")))
    expect_equal(result$summary$yll_se, direct_sandwich_se(arguments), tolerance = 1e-5)
  }
})

test_that("sandwich results keep the common shape and point estimates", {
  for (population in c("all", "exposed", "unexposed")) {
    bootstrap <- estimate_toy(target_population = population)
    sandwich <- estimate_toy_sandwich(target_population = population)
    expect_named(sandwich, names(bootstrap))
    expect_named(sandwich$summary, names(bootstrap$summary))
    expect_equal(sandwich$summary$yll, bootstrap$summary$yll, tolerance = 1e-12)
    expect_equal(sandwich$summary$ci_low,
                 sandwich$summary$yll - qnorm(.975) * sandwich$summary$yll_se)
    expect_true(all(sandwich$summary$yll_se[1:2] > 0))
    expect_equal(tail(sandwich$summary$yll_se, 1), 0)
    expect_equal(nrow(sandwich$bootstrap_estimates), 0)
    expect_equal(nrow(sandwich$bootstrap_survival_curves), 0)
    curves <- sandwich$survival_curves
    expect_true(all(c("survival_reference_se", "survival_exposed_se") %in% names(curves)))
    at_start <- curves$age == curves$starting_age
    expect_true(all(curves$survival_reference_se[at_start] == 0))
    expect_true(all(curves$survival_exposed_se[!at_start] > 0))
    expect_equal(sandwich$meta$variance, "sandwich")
    expect_equal(sandwich$meta$B, 0)
  }
  expect_equal(bootstrap$meta$variance, "bootstrap")
  trapezoidal <- estimate_toy_sandwich(integration = "trapezoidal")
  expect_equal(trapezoidal$summary$yll, estimate_toy(integration = "trapezoidal")$summary$yll,
               tolerance = 1e-12)
  expect_equal(estimate_toy_sandwich(age_start = 70)$summary$yll_se, 0)
})

test_that("splitting the target population into chunks does not change the variance", {
  arguments <- toy_arguments()
  point <- do.call(yll_estimate_logistic_once, c(
    arguments[c("data", "time_var", "event_var", "exposure_var", "reference_level",
                "exposed_level", "age_at_entry_var", "confounders_baseline",
                "age_start", "age_end", "age_interval", "target_population")],
    list(id_var = "id", integration = "left_rectangle", keep_fit = TRUE)
  ))
  sandwich <- function(chunk_rows) {
    yll_sandwich_logistic(
      point$fit, point$person_ids, arguments$data, id_var = "id",
      exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
      confounders_baseline = "sex", target_population = "all", age_start = 50,
      age_end = 70, age_interval = 10, integration = "left_rectangle",
      chunk_rows = chunk_rows
    )
  }
  small <- sandwich(50)
  large <- sandwich(1e7)
  expect_equal(small$standard_errors, large$standard_errors, tolerance = 1e-8)
  expect_equal(small$curves, large$curves, tolerance = 1e-8)
  expect_equal(small$estimates$yll, point$estimates$yll, tolerance = 1e-9)
})

test_that("integration weights reproduce yll_integrate_le()", {
  points <- tibble::tibble(starting_age = c(50, 50, 50, 50, 52, 52, 53), age = c(50:53, 52:53, 53))
  surv <- c(1, 0.9, 0.7, 0.6, 1, 0.8, 1)
  for (integration in c("left_rectangle", "trapezoidal")) {
    weights <- yll_integration_weights(points, integration)
    for (starting_age in unique(points$starting_age)) {
      rows <- points$starting_age == starting_age
      expect_equal(sum(weights[rows] * surv[rows]),
                   yll_integrate_le(points$age[rows], surv[rows], integration))
    }
  }
})

test_that("sandwich works with numeric exposure and no covariates", {
  data <- toy_arguments()$data
  data$hypertension <- as.integer(data$hypertension == "Yes")
  bootstrap <- estimate_toy(data = data, reference_level = 0, exposed_level = 1,
                            confounders_baseline = NULL)
  sandwich <- estimate_toy_sandwich(data = data, reference_level = 0, exposed_level = 1,
                                    confounders_baseline = NULL)
  expect_equal(sandwich$summary$yll, bootstrap$summary$yll, tolerance = 1e-12)
  expect_true(all(is.finite(sandwich$summary$yll_se)))
})

test_that("sandwich settings that cannot apply fail or warn explicitly", {
  expect_error(estimate_toy_sandwich(ci_method = "percentile"),
               "normal approximation")
  expect_warning(estimate_toy(B = 100, variance = "sandwich"), "ignored")
  expect_error(estimate_toy(variance = "jackknife"), "arg")
})
