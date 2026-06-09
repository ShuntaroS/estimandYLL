# Estimate conditional years of life lost by regression methods

Implements model-based observed-group YLL estimators. This is separate
from
[`estimand_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimand_yll.md),
which targets intervention-based estimands. `conditional_yll()` compares
the predicted remaining life expectancy for the reference and exposed
groups, optionally conditioning on baseline covariates.

## Usage

``` r
conditional_yll(
  data,
  method = c("poisson", "flexible_parametric"),
  id_var = "id",
  time_var,
  event_var,
  exposure_var,
  reference_level = NULL,
  exposed_level = NULL,
  age_at_entry_var,
  confounders_baseline = NULL,
  age_start = 50,
  age_end = 100,
  age_interval = 5,
  prediction_interval = 1,
  B = 1000,
  seed = 1,
  conf_level = 0.95,
  integration = c("left_rectangle", "trapezoidal"),
  rp_df = 4,
  rp_model = c("combined", "stratified"),
  show_progress = TRUE
)
```

## Arguments

- data:

  A `data.frame` with one row per individual.

- method:

  `"poisson"` for an age-split Poisson rate model, or
  `"flexible_parametric"` for a Royston-Parmar model via `rstpm2`.

- id_var, time_var, event_var, exposure_var:

  Column names in `data`.

- reference_level, exposed_level:

  Values labelling the binary exposure levels.

- age_at_entry_var:

  Column name for age at study entry.

- confounders_baseline:

  Optional baseline covariates. Predictions are averaged over the
  empirical covariate distribution.

- age_start, age_end, age_interval:

  Starting ages for reported YLL.

- prediction_interval:

  Age-grid spacing used when integrating survival curves.

- B, seed, conf_level, show_progress:

  Bootstrap settings.

- integration:

  Numerical integration rule.

- rp_df, rp_model:

  Flexible parametric model settings used only for
  `method = "flexible_parametric"`.

## Value

A list with `summary`, `detailed_results`, `meta`,
`conditional_survival_point`, and `conditional_survival_boot`.
