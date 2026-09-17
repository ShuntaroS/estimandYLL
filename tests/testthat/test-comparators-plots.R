test_that("comparator results have the same meaning and shape", {
  arguments <- toy_arguments()
  arguments$target_population <- NULL
  arguments$use_future <- NULL
  poisson <- do.call(estimate_yll_poisson, arguments)
  expect_equal(poisson$meta$target_population, "all")
  expect_equal(poisson$summary$yll, poisson$summary$erl_reference - poisson$summary$erl_exposed)
  arguments$age_start <- 70
  expect_equal(do.call(estimate_yll_poisson, arguments)$summary$yll, 0)
  skip_if_not_installed("rstpm2")
  arguments$age_start <- 50
  rp <- do.call(estimate_yll_royston_parmar, arguments)
  expect_named(rp, names(poisson))
  expect_equal(rp$summary$yll, rp$summary$erl_reference - rp$summary$erl_exposed)
  expect_equal(tail(rp$summary$yll, 1), 0)
})

test_that("plots are editable ggplot objects with optional confidence intervals", {
  skip_if_not_installed("ggplot2")
  for (B in c(0, 3)) {
    result <- estimate_toy(B = B)
    expect_s3_class(plot_yll(result), "ggplot")
    expect_s3_class(plot_survival(result, age_start = 50), "ggplot")
    expect_s3_class(plot_yll(result) + ggplot2::labs(title = "Edited"), "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot_survival(result, 50, conf_band = FALSE)))
    expect_no_error(ggplot2::ggplot_build(plot_yll(result)))
    expect_no_error(ggplot2::ggplot_build(plot_survival(result, 70)))
    expect_error(plot_survival(result, 55), "not a reported starting age")
  }
})
