# Estimate years of life lost from a clearly stated estimand

Main user-facing interface for the package. The estimand is specified by
three explicit components: target population, intervention, and measure.

## Usage

``` r
estimand_yll(
  data,
  id_var = "id",
  time_var,
  event_var,
  exposure_var,
  reference_level = NULL,
  exposed_level = NULL,
  age_at_entry_var,
  target_population = c("all", "exposed", "unexposed"),
  intervention = c("full_contrast", "partial_change"),
  change_from = NULL,
  change_to = NULL,
  change_probability = NULL,
  measure = c("yll", "life_year_change"),
  age_start = 50,
  age_end = 100,
  age_interval = 5,
  confounders_baseline = NULL,
  B = 1000,
  seed = 1,
  conf_level = 0.95,
  method = c("normal", "percentile"),
  integration = c("left_rectangle", "trapezoidal"),
  show_progress = TRUE,
  use_future = TRUE
)
estimate_yll(...)
```

## Arguments

- data:

  A `data.frame` with one row per individual.

- id_var, time_var, event_var, exposure_var:

  Column names in `data`.

- reference_level, exposed_level:

  Values labelling the binary exposure levels.

- age_at_entry_var:

  Column name for age at study entry.

- target_population:

  Who to average over: `"all"`, `"exposed"`, `"unexposed"`, or a custom
  function.

- intervention:

  `"full_contrast"` for all-reference versus all-exposed, or
  `"partial_change"` for a policy-like partial change from the natural
  exposure distribution.

- change_from, change_to, change_probability:

  Partial-change intervention settings. Use `"reference"` or `"exposed"`
  for `change_from` and `change_to`.

- measure:

  `"yll"` for `LE_reference - LE_exposed`; `"life_year_change"` for the
  clinically oriented life expectancy change.

- age_start, age_end, age_interval:

  Starting ages for which estimates are reported.

- confounders_baseline:

  Character vector of baseline confounder column names.

- B, seed, conf_level, method, integration, show_progress, use_future:

  Bootstrap and numerical-integration settings.

- ...:

  Arguments passed from the compatibility alias `estimate_yll()` to
  `estimand_yll()`.

## Value

A list containing `summary`, `detailed_results`, `meta`, and marginal
survival curves.
