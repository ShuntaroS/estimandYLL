# Changelog

## estimandYLL 0.1.0

- Replace the development API with
  [`estimate_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimate_yll.md),
  two explicitly named comparator functions, and
  [`plot_yll()`](https://shuntaros.github.io/estimandYLL/reference/plot_yll.md)
  /
  [`plot_survival()`](https://shuntaros.github.io/estimandYLL/reference/plot_survival.md).
- Require the baseline target population for the main estimator. Always
  return reference-minus-exposed restricted expected residual lifetime;
  remove policy interventions, display measures, and estimand-specific
  wrappers.
- Keep the existing mortality models, annual pooled-logistic
  person-period expansion, individual conditioning, integration, and
  participant bootstrap.
- Return one common five-element result with descriptive ERL and
  survival names, bootstrap estimates, and explicit failure records.
  Flag finite but impossible bootstrap lifetimes without silently
  dropping their values.
- Keep normal and percentile intervals for the main estimator. All
  plotting functions return editable ggplot2 objects.
- Validate incomplete data, duplicate IDs and invalid settings before
  fitting; preserve the caller’s random state and parallel/progress
  configuration.
- Replace obsolete documentation with an English tutorial,
  interpretation and comparator articles, a migration guide, and a
  Japanese introduction.
