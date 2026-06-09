# estimandYLL

<!-- badges: start -->
[![R-CMD-check](https://github.com/ShuntaroS/estimandYLL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ShuntaroS/estimandYLL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/ShuntaroS/estimandYLL/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ShuntaroS/estimandYLL/actions/workflows/pkgdown.yaml)
<!-- badges: end -->

Years of life lost (YLL) and counterfactual life expectancy under
user-specified interventions on a binary exposure, estimated with the
parametric **g-formula** on the attained-age time scale.

The package supports:

- Explicit target populations: all individuals, observed-exposed
  individuals, or observed-unexposed individuals.
- Full exposure contrasts and clinically phrased partial-change policies
  through one main function, `estimand_yll()`.
- Conditional, observed-group YLL by Poisson or Royston-Parmar regression
  through `conditional_yll()`.
- Bootstrap confidence intervals (normal-approximation or percentile),
  optionally parallelised with `future`.
- Plotting helpers for marginal and conditional survival curves and for
  the selected YLL measure across starting ages.

## Documentation

Full reference and vignettes are available on the package website:

**<https://shuntaros.github.io/estimandYLL/>**

- [Function reference](https://shuntaros.github.io/estimandYLL/reference/index.html)
- [Getting started](https://shuntaros.github.io/estimandYLL/articles/getting-started.html)
- [Generalized interventional effects](https://shuntaros.github.io/estimandYLL/articles/interventional-effects.html)

## Installation

The package is currently distributed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ShuntaroS/estimandYLL")
```

## Quick start

```r
library(estimandYLL)
data(yll_toy)

res <- estimand_yll(
  data = yll_toy,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  target_population = "all",
  intervention = "full_contrast",
  measure = "yll",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  B = 200,
  method = "normal",
  use_future = FALSE
)

res$summary
plot_yll_estimate(res)
```

See `vignette("getting-started", package = "estimandYLL")` for a tour of
the user-facing function and the estimand definition.

## Result object

`estimand_yll()` returns a list with:

| Element                    | Description                                                                  |
| -------------------------- | ---------------------------------------------------------------------------- |
| `summary`                  | One row per starting age: selected estimate, CI, measure, and CI method.     |
| `detailed_results`         | Same plus per-arm life expectancy and both percentile and normal CIs.        |
| `meta`                     | Call metadata (`B`, `seed`, `conf_level`, `target_population`, `intervention`, `measure`, …). |
| `marginal_survival_point`  | Population-level marginal survival curves under each intervention arm.       |
| `marginal_survival_boot`   | Same, one row per bootstrap iteration (used by the plotting helpers).        |

`method` selects between normal-approximation (`"normal"`, default) and
percentile (`"percentile"`) bootstrap CIs; `res$summary$ci_method`
records which one was used.

## Estimand at a glance

The main design is:

```text
YLL estimand = target_population + intervention + measure
```

`target_population` answers "who are we averaging over?":

- `"all"`: all individuals (ATE-like).
- `"exposed"`: individuals observed exposed (ATT-like).
- `"unexposed"`: individuals observed unexposed (ATC-like).

`intervention` answers "what exposure worlds are compared?":

- `"full_contrast"`: all reference versus all exposed.
- `"partial_change"`: the natural observed distribution versus a policy-like
  partial change, such as 30% of observed-exposed individuals moving to the
  reference level.

`measure` answers "how is the result displayed?":

- `"yll"`: `LE_reference - LE_exposed`.
- `"life_year_change"`: the clinically oriented change in life expectancy
  for the stated intervention direction. This makes harmful ATC-like changes
  negative.

Example: 30% of smokers quit, averaged over observed smokers:

```r
estimand_yll(
  data = yll_toy,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "smoke_binary",
  reference_level = "Never",
  exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  target_population = "exposed",
  intervention = "partial_change",
  change_from = "exposed",
  change_to = "reference",
  change_probability = 0.30,
  measure = "life_year_change",
  B = 0,
  use_future = FALSE
)
```

In statistical terminology, `partial_change` is a stochastic intervention,
but the user-facing API uses policy language.

Plot helpers:

- `plot_yll_estimate()` / `plot_yll()` — selected estimate across starting ages.
- `plot_marginal_survival()` and `plot_conditional_survival()` — survival curves.

## Conditional YLL

`conditional_yll()` implements regression-based observed-group comparisons
similar to the Poisson and flexible parametric Royston-Parmar approaches
used in the YLL methods literature. It is separate from `estimand_yll()`:
there is no intervention or target population argument. Instead, the function
compares predicted remaining life expectancy for the reference and exposed
groups, optionally conditioning on baseline covariates.

```r
res_cond <- conditional_yll(
  data = yll_toy,
  method = "poisson",
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  B = 0
)
```

`method = "flexible_parametric"` requires the optional `rstpm2` package.

## Caveats

- The current implementation is for **binary** exposures. Multi-level
  and continuous exposures are not yet supported.
- Variable entry age introduces a left-truncation / delayed-entry
  problem that the current implementation does not address rigorously;
  treat results from datasets with variable `age_at_entry_var` as
  preliminary.

## License

MIT © 2026 Shuntaro Sato.
