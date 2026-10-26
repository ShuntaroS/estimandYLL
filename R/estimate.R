#' Estimate years of life lost in a specified target population
#'
#' Fit a pooled logistic hazard model on the full analytic sample, predict
#' survival under both exposure levels, and average individual conditional
#' survival curves over the specified baseline population. The contrast is
#' always restricted expected residual lifetime (ERL) under the reference
#' level minus ERL under the exposed level. Its sign never depends on the
#' target population. Confidence intervals come from a participant bootstrap
#' or from the empirical sandwich variance.
#'
#' @param data A data frame with one row per person and no missing values in
#'   the selected variables. Remove incomplete records explicitly before use.
#' @param id_var Name of the unique person identifier column.
#' @param time_var Name of the positive follow-up duration column, in years.
#' @param event_var Name of the event indicator column, coded 0/1.
#' @param exposure_var Name of the binary exposure column.
#' @param reference_level,exposed_level Observed exposure values corresponding
#'   to x = 0 and x = 1. Specify both to make the direction explicit. If omitted,
#'   factor order (or order of appearance for other types) is used and reported.
#' @param age_at_entry_var Name of the age-at-entry column, in years.
#' @param target_population Required. One of `"all"`, `"exposed"`, or
#'   `"unexposed"`. Selects the baseline covariate distribution used for
#'   standardization, not the data used to fit the hazard model.
#' @param age_start First starting age at which to report ERL and YLL.
#' @param age_end Upper age limiting every ERL integral. Also the upper bound
#'   of the reporting sequence. Must be at least `age_start`.
#' @param age_interval Spacing between reported starting ages. Results are
#'   reported at `seq(age_start, age_end, age_interval)`. The main estimator
#'   requires integer ages and uses a one-year integration grid regardless
#'   of this reporting interval.
#' @param confounders_baseline Character vector of baseline covariate columns.
#' @param B Number of participant bootstrap replicates, at least 2, or 0 to
#'   obtain point estimates only. Failed replicates are reported in a warning
#'   and recorded in the result metadata; intervals use successful replicates.
#' @param seed Integer random seed. The caller's random state is restored.
#' @param conf_level Confidence level, strictly between zero and one.
#' @param ci_method `"normal"` (default) uses the estimated standard error;
#'   `"percentile"` uses bootstrap quantiles and is available only with
#'   `variance = "bootstrap"`.
#' @param integration `"left_rectangle"` (default) or `"trapezoidal"`.
#' @param show_progress Whether to show bootstrap progress. Parallel progress
#'   uses the caller's progressr handlers.
#' @param use_future Use the caller's future plan for bootstrapping. The
#'   package never changes the plan or starts a worker pool itself.
#' @param variance `"bootstrap"` (default) resamples participants `B` times.
#'   `"sandwich"` fits the model once and uses the empirical sandwich variance
#'   of the stacked estimating equations; it ignores `B`, `seed`,
#'   `show_progress`, and `use_future`, and requires `ci_method = "normal"`.
#'
#' @return A list with five elements:
#' * `summary`: one row per `starting_age`, with `erl_reference`, `erl_exposed`,
#'   `yll`, `yll_se`, `ci_low`, and `ci_high`, all measured in years.
#' * `bootstrap_estimates`: the four point-estimate columns plus `iteration`.
#' * `survival_curves`: `starting_age`, `age`, `survival_reference`, and
#'   `survival_exposed`, averaged over the target population. With
#'   `variance = "sandwich"`, `survival_reference_se` and `survival_exposed_se`
#'   are added.
#' * `bootstrap_survival_curves`: the same curve columns plus `iteration`.
#' * `meta`: model, target population, exposure labels, age window, variance
#'   and interval methods, bootstrap settings, successful replicate count, and
#'   failure details.
#'   `bootstrap_unreliable` records finite bootstrap ERLs outside the possible
#'   age window. These trigger a warning and are retained, not silently removed.
#'
#' With `B = 0` and the bootstrap, intervals and standard errors are `NA` and
#' bootstrap tables are empty. With `variance = "sandwich"`, the bootstrap
#' tables are always empty. Curves are standardized after conditioning each person on
#' survival to the starting age. They are not curves obtained by first
#' averaging unconditional survival and then conditioning that average.
#' Intervals are pointwise, not simultaneous confidence bands.
#'
#' @details The fitted model includes exposure, a natural cubic spline in
#' attained age, their interaction, and the supplied baseline covariates.
#' Follow-up is split into one-year intervals after cohort entry. Predictions
#' represent sustained exposure scenarios, not exposure changes initiated at
#' the starting age. Causal interpretation requires consistency, exchangeability,
#' positivity, appropriate censoring and entry assumptions, and correct model
#' specification. Predictions beyond observed age support require extrapolation.
#' The same baseline covariate distribution is retained at every starting age.
#'
#' The sandwich variance stacks the pooled logistic score, summed within
#' person, with the estimating equations for the standardized survival
#' probabilities, following Zivich et al. (2025). It accounts for both the
#' fitted hazard model and the baseline covariate distribution of the target
#' population. Spline knots are treated as fixed, whereas each bootstrap
#' replicate selects them again. Predictions are processed in groups of people,
#' so their memory use does not grow with the size of the target population.
#'
#' @references
#' Zivich PN, Cole SR, Shook-Sa BE, DeMonte JB, Edwards JK (2025). Estimating
#' equations for causal survival analysis with pooled logistic regression.
#' arXiv:2504.13291. \doi{10.48550/arXiv.2504.13291}
#'
#' @examples
#' data(yll_toy)
#' result <- estimate_yll(
#'   yll_toy[1:500, ], time_var = "period", event_var = "event",
#'   exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
#'   age_at_entry_var = "age", target_population = "unexposed",
#'   age_start = 50, age_end = 70, age_interval = 10,
#'   confounders_baseline = "sex", B = 0, use_future = FALSE
#' )
#' result$summary
#'
#' # Sandwich confidence intervals from a single model fit
#' sandwich <- estimate_yll(
#'   yll_toy[1:500, ], time_var = "period", event_var = "event",
#'   exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
#'   age_at_entry_var = "age", target_population = "unexposed",
#'   age_start = 50, age_end = 70, age_interval = 10,
#'   confounders_baseline = "sex", variance = "sandwich"
#' )
#' sandwich$summary
#' @export
estimate_yll <- function(
    data,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    target_population,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    ci_method = c("normal", "percentile"),
    integration = c("left_rectangle", "trapezoidal"),
    show_progress = TRUE,
    use_future = TRUE,
    variance = c("bootstrap", "sandwich")
) {
  variance <- match.arg(variance)
  if (variance == "sandwich" && !missing(B)) {
    warning("`B` is ignored when variance = \"sandwich\".", call. = FALSE)
  }
  yll_estimate(
    data = data, id_var = id_var, time_var = time_var, event_var = event_var,
    exposure_var = exposure_var, reference_level = reference_level,
    exposed_level = exposed_level, age_at_entry_var = age_at_entry_var,
    age_start = age_start, age_end = age_end, age_interval = age_interval,
    confounders_baseline = confounders_baseline, B = B, seed = seed,
    conf_level = conf_level, integration = match.arg(integration),
    show_progress = show_progress,
    model = "pooled_logistic", target_population = match.arg(target_population,
      c("all", "exposed", "unexposed")),
    ci_method = match.arg(ci_method), use_future = use_future, variance = variance
  )
}

#' Estimate standardized YLL using a Poisson model
#'
#' Predict both exposure scenarios for every participant and average the
#' individual conditional survival curves over the full analytic sample.
#' This comparator changes the mortality model, not the direction of YLL.
#' Confidence intervals use the normal approximation to the participant
#' bootstrap distribution. They are pointwise intervals.
#' @inheritParams estimate_yll
#' @param prediction_interval Spacing of prediction ages, in years. For the
#'   Poisson estimator this also determines the age-splitting interval.
#' @return The same five-element list as [estimate_yll()], with
#'   the full sample as target population and normal-approximation intervals.
#' @seealso [estimate_yll()], [plot_yll()], [plot_survival()]
#' @export
estimate_yll_poisson <- function(
    data,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    prediction_interval = 1,
    integration = c("left_rectangle", "trapezoidal"),
    show_progress = TRUE
) {
  yll_estimate(
    data = data, id_var = id_var, time_var = time_var, event_var = event_var,
    exposure_var = exposure_var, reference_level = reference_level,
    exposed_level = exposed_level, age_at_entry_var = age_at_entry_var,
    age_start = age_start, age_end = age_end, age_interval = age_interval,
    confounders_baseline = confounders_baseline, B = B, seed = seed,
    conf_level = conf_level, integration = match.arg(integration),
    show_progress = show_progress,
    model = "poisson", target_population = "all",
    ci_method = "normal", use_future = FALSE, prediction_interval = prediction_interval
  )
}

#' Estimate standardized YLL using a Royston-Parmar model
#'
#' Predict both exposure scenarios for every participant and average the
#' individual conditional survival curves over the full analytic sample.
#' This comparator changes the mortality model, not the direction of YLL.
#' Confidence intervals use the normal approximation to the participant
#' bootstrap distribution. They are pointwise intervals.
#' @inheritParams estimate_yll
#' @param prediction_interval Spacing of prediction ages, in years. For the
#'   Poisson estimator this also determines the age-splitting interval.
#' @param rp_df Degrees of freedom for the baseline log cumulative hazard
#'   spline, passed to rstpm2::stpm2().
#' @param rp_model `"combined"` fits one model including exposure; `"stratified"`
#'   fits a separate model in each exposure group. The optional rstpm2 package
#'   is required. Inspect warnings and bootstrap failures for unstable fits.
#' @return The same five-element list as [estimate_yll()], with
#'   the full sample as target population and normal-approximation intervals.
#' @seealso [estimate_yll()], [plot_yll()], [plot_survival()]
#' @export
estimate_yll_royston_parmar <- function(
    data,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    prediction_interval = 1,
    rp_df = 4,
    rp_model = c("combined", "stratified"),
    integration = c("left_rectangle", "trapezoidal"),
    show_progress = TRUE
) {
  yll_estimate(
    data = data, id_var = id_var, time_var = time_var, event_var = event_var,
    exposure_var = exposure_var, reference_level = reference_level,
    exposed_level = exposed_level, age_at_entry_var = age_at_entry_var,
    age_start = age_start, age_end = age_end, age_interval = age_interval,
    confounders_baseline = confounders_baseline, B = B, seed = seed,
    conf_level = conf_level, integration = match.arg(integration),
    show_progress = show_progress,
    model = "royston_parmar", target_population = "all",
    ci_method = "normal", use_future = FALSE, prediction_interval = prediction_interval,
    rp_df = rp_df, rp_model = match.arg(rp_model)
  )
}
