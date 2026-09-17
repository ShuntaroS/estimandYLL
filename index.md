# estimandYLL

[![R-CMD-check](https://github.com/ShuntaroS/estimandYLL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ShuntaroS/estimandYLL/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/ShuntaroS/estimandYLL/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/ShuntaroS/estimandYLL/actions/workflows/pkgdown.yaml)

Estimate **years of life lost (YLL)** between two binary exposure
scenarios, for a stated baseline population and age window.

**YLL = ERL under the reference level − ERL under the exposed level.**
ERL is expected residual lifetime restricted at the upper age. The
contrast has the same direction for every target population; the
interpretation changes with the population being studied.

## Install

Version 0.1.0 is distributed from GitHub. It has not been submitted to
CRAN.

``` r

install.packages("remotes")
remotes::install_github("ShuntaroS/estimandYLL")
```

## First estimate

``` r

library(estimandYLL)
data(yll_toy)
result <- estimate_yll(
  data = yll_toy[1:500, ],
  time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age", target_population = "unexposed",
  age_start = 50, age_end = 70, age_interval = 10,
  confounders_baseline = "sex", B = 0, use_future = FALSE
)
result$summary
```

This quick example uses 500 simulated people and skips confidence
intervals. For an analysis, choose the age window, covariates and target
population for your question, and use an adequate number of bootstrap
replicates, for example `B = 2000`. The data are synthetic; this example
does not reproduce the paper’s cohort results. No real participant
records are distributed.

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
| `ci_method = "normal"` | Bootstrap standard error with a normal approximation (default). |
| `ci_method = "percentile"` | Quantiles of the same participant bootstrap distribution. |

The target population is required. Both interval methods resample people
with replacement. Neither changes the fitted mortality model. Intervals
are pointwise.

## Results and plots

`result$summary` contains `starting_age`, `erl_reference`,
`erl_exposed`, `yll`, `yll_se`, `ci_low` and `ci_high`. Lifetimes, their
difference and uncertainty are measured in years. The remaining elements
are `bootstrap_estimates`, `survival_curves`,
`bootstrap_survival_curves`, and `meta`.

``` r

library(ggplot2)
plot_yll(result) + labs(title = "YLL among participants without hypertension")
plot_survival(result, age_start = 50) + theme_bw()
# ggsave("yll.png", plot_yll(result), width = 6, height = 4, dpi = 300)
```

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

See [NEWS](https://shuntaros.github.io/estimandYLL/NEWS.md) for the
breaking changes in 0.1.0 and [the validation
record](https://github.com/ShuntaroS/estimandYLL/blob/main/notes/validation/README.md)
for the 100-replicate comparison with the previous development version.
Source scripts are provided there.

## License

MIT © 2026 Shuntaro Sato.
