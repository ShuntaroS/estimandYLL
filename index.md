# estimandYLL

[![R-CMD-check](https://github.com/ShuntaroS/estimandYLL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ShuntaroS/estimandYLL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/ShuntaroS/estimandYLL/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ShuntaroS/estimandYLL/actions/workflows/pkgdown.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23653505.svg)](https://doi.org/10.5281/zenodo.23653505)

Estimate **years of life lost (YLL)** between two binary exposure
scenarios, for a stated baseline population and age window.

**YLL = ERL under the reference level − ERL under the exposed level.**
ERL is expected residual lifetime restricted at the upper age. The
contrast has the same direction for every target population; the
interpretation changes with the population being studied.

## Install

Version 0.2.0 is distributed from GitHub. It has not been submitted to
CRAN.

``` r

install.packages("remotes")
remotes::install_github("ShuntaroS/estimandYLL")
```

## First estimate

``` r

library(estimandYLL)
data(yll_toy)
previous_plan <- future::plan(future::multisession, workers = 2)
result <- estimate_yll(
  data = yll_toy,
  time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age", target_population = "unexposed",
  age_start = 40, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  B = 1000, seed = 20260917, ci_method = "normal", use_future = TRUE
)
future::plan(previous_plan)
result$summary
```

This example uses all 5,000 simulated people, 1,000 participant
bootstrap replicates, and starting ages from 40 to 90 in five-year
steps. The figures below include 95% pointwise confidence intervals. The
data are synthetic; this example does not reproduce the paper’s cohort
results. No real participant records are distributed. Choose the age
window, covariates, target population, and number of replicates for your
own analysis.

The main estimator uses the caller’s `future` plan when
`use_future = TRUE` (the default). The code above requests two separate
R workers and then restores the previous plan. Merely setting
`use_future = TRUE` does not start workers.

To inspect the published figures without rerunning the bootstrap, load
the saved result included with the package. The full reproduction script
is in
[notes/published-example](https://github.com/ShuntaroS/estimandYLL/tree/main/notes/published-example).

``` r

result <- readRDS(system.file("extdata", "yll-example.rds", package = "estimandYLL"))
```

### Intervals from a single model fit

`variance = "sandwich"` replaces the bootstrap with the empirical
sandwich variance of stacked estimating equations ([Zivich et al.,
2025](https://doi.org/10.48550/arXiv.2504.13291)). The model is fitted
once, so the same analysis returns intervals in a few seconds instead of
refitting 1,000 models. On this example the sandwich standard errors
were within 2% of the bootstrap standard errors.

``` r

result_sandwich <- estimate_yll(
  data = yll_toy,
  time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age", target_population = "unexposed",
  age_start = 40, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  variance = "sandwich"
)
```

`age_end` limits every integral. `age_interval` controls the spacing
between reported starting ages; the main estimator still integrates on a
one-year grid.

## What to specify

| Argument | Meaning |
|----|----|
| `target_population = "all"` | Average over all participants (ATE). |
| `target_population = "exposed"` | Average over those observed exposed (ATT). |
| `target_population = "unexposed"` | Average over those observed unexposed (ATU). |
| `reference_level`, `exposed_level` | Observed values corresponding to x = 0 and x = 1. |
| `variance = "bootstrap"` | Participant bootstrap with `B` replicates (default). |
| `variance = "sandwich"` | Empirical sandwich variance from a single model fit. |
| `ci_method = "normal"` | Standard error with a normal approximation (default). |
| `ci_method = "percentile"` | Quantiles of the bootstrap distribution; bootstrap only. |

The target population is required. Neither variance method changes the
fitted mortality model. Intervals are pointwise.

## Results and plots

`result$summary` contains `starting_age`, `erl_reference`,
`erl_exposed`, `yll`, `yll_se`, `ci_low` and `ci_high`. Lifetimes, their
difference and uncertainty are measured in years. The remaining elements
are `bootstrap_estimates`, `survival_curves`,
`bootstrap_survival_curves`, and `meta`.

``` r

library(ggplot2)
plot_yll(result, conf_band = TRUE) +
  labs(title = "YLL among participants without hypertension")
plot_survival(result, age_start = 50) + theme_bw()
# ggsave("yll.png", plot_yll(result), width = 6, height = 4, dpi = 300)
```

![YLL with 95% pointwise confidence intervals from 1,000 bootstrap
replicates](reference/figures/yll-example.png)

YLL with 95% pointwise confidence intervals from 1,000 bootstrap
replicates

Points show YLL and error bars show 95% confidence intervals, at
five-year starting-age intervals. At age 90 the integration window is
empty, so YLL and its interval are zero.

![Standardized survival from age 50 with 95% pointwise confidence
intervals](reference/figures/survival-example.png)

Standardized survival from age 50 with 95% pointwise confidence
intervals

Survival is restarted at one at age 50. Shaded areas show 95% pointwise
confidence intervals for each exposure scenario.

All plots are editable ggplot objects. The optional ggplot2 package is
needed only for plotting. Poisson and Royston-Parmar comparators are
available through
[`estimate_yll_poisson()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll_poisson.md)
and
[`estimate_yll_royston_parmar()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll_royston_parmar.md);
they standardize over the full sample. The latter requires the optional
rstpm2 package.

## Read the documentation

- [Getting
  started](https://shuntaros.github.io/estimandYLL/articles/getting-started.html)
- [Target populations and
  interpretation](https://shuntaros.github.io/estimandYLL/articles/target-populations.html)
- [Comparison
  models](https://shuntaros.github.io/estimandYLL/articles/comparison-models.html)
- [Migration from the development
  API](https://shuntaros.github.io/estimandYLL/articles/migration.html)
- [日本語の入門](https://shuntaros.github.io/estimandYLL/articles/getting-started-ja.html)
- [Function
  reference](https://shuntaros.github.io/estimandYLL/reference/index.html)

## Interpretation and limitations

The model uses observed person-period records after cohort entry. At
each starting age, individual survival curves are restarted at one and
then averaged over the target population’s baseline covariates. This
does not select only people observed alive at that age and is not a
principal-stratum effect.

A causal interpretation requires consistency, conditional
exchangeability, positivity, appropriate censoring and entry
assumptions, and an adequate model. The scenarios describe sustained
exposure, not a change initiated at the starting age. Age extrapolation
and sparse exposure-covariate combinations can make estimates sensitive
to the model. Read the interpretation article before applying the
package to an observational study.

## Development and validation

See [NEWS](https://shuntaros.github.io/estimandYLL/news/index.html) for
the changes in each version. [The validation
record](https://github.com/ShuntaroS/estimandYLL/blob/main/notes/validation/README.md)
gives the 100-replicate comparison with the previous development
version, and [the sandwich
validation](https://github.com/ShuntaroS/estimandYLL/blob/main/notes/sandwich-validation/README.md)
compares the sandwich variance with the bootstrap. Source scripts are
provided there.

## License

MIT © 2026 Shuntaro Sato.
