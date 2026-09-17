# Migrating from the development API

Version 0.1.0 deliberately replaces the earlier development interface.
Old names and arguments are removed rather than retained as
compatibility wrappers.

| Earlier interface | Version 0.1.0 |
|----|----|
| `estimand_yll()` | [`estimate_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll.md) |
| [`estimate_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll.md) as an alias | The main function with explicit arguments. |
| Optional target population | Required `target_population`. |
| `intervention = "full_contrast"` | Always compare the two specified exposure levels; omit this argument. |
| `partial_change` and `change_*` | Removed. Population policy effects are outside this release. |
| `measure` / `life_year_change` | Removed. YLL always equals reference ERL minus exposed ERL. |
| `method` for the main function’s interval | `ci_method`. |
| `conditional_yll(method = "poisson")` | [`estimate_yll_poisson()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll_poisson.md) |
| `conditional_yll(method = "flexible_parametric")` | [`estimate_yll_royston_parmar()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll_royston_parmar.md) |
| `plot_yll_estimate(res)` / `plot_yll(res)` | `plot_yll(result)` |
| `plot_conditional_survival(res, age_start)` | `plot_survival(result, age_start)` |
| `plot_marginal_survival(res)` | `plot_survival(result, age_start = result$meta$age_start)` |

The former internal ATE/ATT/ATC wrapper functions are removed. Use
`"all"`, `"exposed"`, or `"unexposed"`; the manuscript uses ATU for the
last population. No population-specific sign reversal is applied.

## Result columns

`summary$yll` replaces the redundant `estimate` column. `le_reference`
and `le_exposed` become `erl_reference` and `erl_exposed` to identify
restricted expected residual lifetime. `summary` now includes `yll_se`.
Settings common to all rows, such as target population and interval
method, are in `meta`.

`detailed_results` is removed. Replicate estimates are available
directly in `bootstrap_estimates`. All starting-age survival curves are
in `survival_curves` and `bootstrap_survival_curves`, with descriptive
survival column names and an `iteration` column on bootstrap tables. The
separate marginal-curve copies are removed because they duplicated
curves at the first starting age.

## Validation and failure reporting

Selected columns must be complete, IDs must be unique, and the main
estimator requires integer ages. `B` is zero or at least two. A
zero-width age window is allowed and has zero ERL and YLL. `B = 0` is an
ordinary point-estimation request and does not generate the earlier
missing-interval warning.

Bootstrap failures are now reported in warnings and in
`meta$bootstrap_failures`, rather than being silently dropped by the
comparison estimators. Intervals use successful replicates; fewer than
two successes give missing intervals. Random state and the caller’s
future plan are preserved.
