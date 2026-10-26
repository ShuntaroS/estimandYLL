# Validation of the sandwich variance in estimandYLL 0.2.0

`variance = "sandwich"` was compared with the participant bootstrap on the
packaged `yll_toy` data. Its standard errors were also checked against the
spread of estimates over repeated synthetic datasets. The method follows
Zivich et al. (2025), arXiv:2504.13291.

## Settings

Validation date: 2026-10-26. Runtime: R 4.6.1 on macOS (Apple silicon).
All results come from `run-validation.R`, run from the package root.

* Data: all 5,000 people in `yll_toy`. Repeated datasets use the generating
  process in `data-raw/yll_toy.R`, with 5,000 people each.
* Settings: smoking (reference `Never`, exposed `Current/Ever`) standardized to
  the exposed and to all participants; hypertension (reference `No`, exposed
  `Yes`) standardized to the unexposed.
* Covariates: sex, education_years, bmi. Starting ages 40, 50, ..., 90 years;
  upper age 90 years; left-rectangle integration.
* Bootstrap: 1,000 participant resamples, seed 20261026, 14 forked workers.
* Published example: the saved 1,000-replicate result in
  `inst/extdata/yll-example.rds` (hypertension, unexposed, four covariates,
  five-year starting ages).
* Simulation: 500 datasets.

## Results

Ranges cover starting ages 40 to 80 years. At age 90 the window is empty and
both standard errors are zero.

| Setting | Sandwich / bootstrap SE | Sandwich SE / SD of estimates over 500 datasets | Sandwich, seconds | Bootstrap, seconds (14 workers) |
|---|---:|---:|---:|---:|
| Smoking, exposed | 1.00 to 1.01 | 0.97 to 1.02 | 1.4 | 101.5 |
| Smoking, all | 1.00 to 1.01 | 0.97 to 1.02 | 2.6 | 188.5 |
| Hypertension, unexposed | 0.97 to 1.00 | 0.99 to 1.02 | 1.5 | 137.0 |

A 1,000-replicate bootstrap standard error has a Monte Carlo error of about
2%, and the ratio over 500 datasets has an error of about 3%. All ratios are
within these limits. At starting age 80 the simulation ratio was 0.97 to 0.99
in every setting. That window extends beyond the oldest attained age in the
data (84 years), so the estimates there rely on extrapolation.

For the published example, the sandwich standard errors were 0.98 to 1.02
times the saved bootstrap standard errors, and the point estimates were
identical. Full values are in `bootstrap-comparison.csv`,
`published-example-comparison.csv`, and `simulation-summary.csv`.

Coverage was not evaluated. The fitted model omits the smoking-by-hypertension
interaction in the generating process, so the value it targets would have to
be computed separately.

Sandwich times include the point estimate. Timing columns describe these runs
only and are not a benchmark.

## Unit tests

`tests/testthat/test-sandwich.R` compares the chunked covariance with a direct
influence-function calculation that uses `predict()` and numerical derivatives
(agreement within 1e-5). It also checks that splitting the target population
into chunks does not change the result, that point estimates equal those of
the bootstrap path, and that numeric exposures and models without covariates
are supported.
