# Estimate standardized YLL using a Royston-Parmar model

Predict both exposure scenarios for every participant and average the
individual conditional survival curves over the full analytic sample.
This comparator changes the mortality model, not the direction of YLL.
Confidence intervals use the normal approximation to the participant
bootstrap distribution. They are pointwise intervals.

## Usage

``` r
estimate_yll_royston_parmar(
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
)
```

## Arguments

- data:

  A data frame with one row per person and no missing values in the
  selected variables. Remove incomplete records explicitly before use.

- id_var:

  Name of the unique person identifier column.

- time_var:

  Name of the positive follow-up duration column, in years.

- event_var:

  Name of the event indicator column, coded 0/1.

- exposure_var:

  Name of the binary exposure column.

- reference_level, exposed_level:

  Observed exposure values corresponding to x = 0 and x = 1. Specify
  both to make the direction explicit. If omitted, factor order (or
  order of appearance for other types) is used and reported.

- age_at_entry_var:

  Name of the age-at-entry column, in years.

- age_start:

  First starting age at which to report ERL and YLL.

- age_end:

  Upper age limiting every ERL integral. Also the upper bound of the
  reporting sequence. Must be at least `age_start`.

- age_interval:

  Spacing between reported starting ages. Results are reported at
  `seq(age_start, age_end, age_interval)`. The main estimator requires
  integer ages and uses a one-year integration grid regardless of this
  reporting interval.

- confounders_baseline:

  Character vector of baseline covariate columns.

- B:

  Number of participant bootstrap replicates, at least 2, or 0 to obtain
  point estimates only. Failed replicates are reported in a warning and
  recorded in the result metadata; intervals use successful replicates.

- seed:

  Integer random seed. The caller's random state is restored.

- conf_level:

  Confidence level, strictly between zero and one.

- prediction_interval:

  Spacing of prediction ages, in years. For the Poisson estimator this
  also determines the age-splitting interval.

- rp_df:

  Degrees of freedom for the baseline log cumulative hazard spline,
  passed to rstpm2::stpm2().

- rp_model:

  `"combined"` fits one model including exposure; `"stratified"` fits a
  separate model in each exposure group. The optional rstpm2 package is
  required. Inspect warnings and bootstrap failures for unstable fits.

- integration:

  `"left_rectangle"` (default) or `"trapezoidal"`.

- show_progress:

  Whether to show bootstrap progress. Parallel progress uses the
  caller's progressr handlers.

## Value

The same five-element list as
[`estimate_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll.md),
with the full sample as target population and normal-approximation
intervals.

## See also

[`estimate_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll.md),
[`plot_yll()`](https://shuntaros.github.io/estimandYLL/reference/plot_yll.md),
[`plot_survival()`](https://shuntaros.github.io/estimandYLL/reference/plot_survival.md)
