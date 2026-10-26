# estimandYLL 0.2.0

* Add `variance = "sandwich"` to `estimate_yll()`. The model is fitted once,
  and intervals use the empirical sandwich variance of the pooled logistic
  score stacked with the standardization estimating equations (Zivich et al.,
  2025, arXiv:2504.13291). It accounts for both the hazard model and the target
  population's covariate distribution. Spline knots are treated as fixed.
* On `yll_toy`, sandwich standard errors were 0.97 to 1.02 times those of a
  1,000-replicate bootstrap and of the spread over 500 simulated datasets. One
  call took 1 to 3 seconds; see `notes/sandwich-validation`.
* With `variance = "sandwich"`, `survival_curves` gains `survival_reference_se`
  and `survival_exposed_se`, and `plot_survival()` draws bands from them. Only
  normal-approximation intervals are available.
* The variance is accumulated over groups of people, so memory use for
  predictions does not grow with the size of the target population.
* The default remains `variance = "bootstrap"`, so 0.1.0 results are reproduced
  unchanged. `meta$variance` records the method for every estimator.

# estimandYLL 0.1.0

* Archived on Zenodo as doi:10.5281/zenodo.23653506.

* Replace the development API with `estimate_yll()`, two explicitly named
  comparator functions, and `plot_yll()` / `plot_survival()`.
* Require the baseline target population for the main estimator. Always return
  reference-minus-exposed restricted expected residual lifetime; remove policy
  interventions, display measures, and estimand-specific wrappers.
* Keep the existing mortality models, annual pooled-logistic person-period
  expansion, individual conditioning, integration, and participant bootstrap.
* Return one common five-element result with descriptive ERL and survival names,
  bootstrap estimates, and explicit failure records. Flag finite but impossible
  bootstrap lifetimes without silently dropping their values.
* Keep normal and percentile intervals for the main estimator. All plotting
  functions return editable ggplot2 objects.
* Validate incomplete data, duplicate IDs and invalid settings before fitting;
  preserve the caller's random state and parallel/progress configuration.
* Replace obsolete documentation with an English tutorial, interpretation and
  comparator articles, a migration guide, and a Japanese introduction.
* Show public YLL examples at five-year starting-age intervals with 95% confidence
  intervals from 1,000 bootstrap replicates. Include the complete saved synthetic
  result and reproduction scripts, including a two-worker future check.
