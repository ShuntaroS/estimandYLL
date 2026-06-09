# Interventional effects with estimandYLL

This vignette explains how `estimandYLL` represents interventional
effects. The package emphasizes that the estimand is not only an
exposure contrast; it also requires a **target population**.

``` text
Estimand = target_population + intervention + measure
```

## Target population

`target_population` states whose counterfactual life expectancies are
averaged:

- `"all"`: everyone in the study population.
- `"exposed"`: individuals observed at the exposed level.
- `"unexposed"`: individuals observed at the reference level.

The hazard model is fit in the full dataset. The target population
changes only the averaging step after counterfactual survival curves
have been built.

## Interventions

`intervention = "full_contrast"` compares two static worlds:

``` math
\mathrm{LE}^{a=\mathrm{reference}}
\quad \text{versus} \quad
\mathrm{LE}^{a=\mathrm{exposed}}.
```

`intervention = "partial_change"` compares the natural observed exposure
distribution with a policy-like partial change. For example:

- observed-exposed people: 30% move to reference, 70% stay exposed;
- observed-unexposed people: stay unexposed.

Internally, this is a binary stochastic intervention. For each subject
the engine forms

``` math
S_i^*(t) =
(1 - p_i) S_i^{a=\mathrm{reference}}(t)
+ p_i S_i^{a=\mathrm{exposed}}(t),
```

then averages $`S_i^*(t)`$ in the stated target population.

## Measure

`measure = "yll"` reports

``` math
\mathrm{LE}_{reference} - \mathrm{LE}_{exposed}.
```

This is the usual years-of-life-lost orientation. A harmful exposure
tends to produce a positive value.

`measure = "life_year_change"` reports the clinically oriented change in
life expectancy for the stated intervention direction. For an ATC-like
question in which observed-unexposed people become exposed, a harmful
exposure therefore gives a negative value.

## Worked example

Imagine a programme where 20% of observed-unexposed individuals become
hypertensive, and observed-hypertensive individuals stay hypertensive.
The target population is explicitly the observed-unexposed group.

``` r

library(estimandYLL)
data(yll_toy)

res <- estimand_yll(
  data = yll_toy,
  target_population = "unexposed",
  intervention = "partial_change",
  change_from = "reference",
  change_to = "exposed",
  change_probability = 0.20,
  measure = "life_year_change",
  B = 50, show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)

res$summary
```

## Identification assumptions

The usual g-formula assumptions are required:

- consistency of the counterfactual outcomes under the hypothesised
  regimes;
- conditional exchangeability given the measured covariates;
- positivity for the chosen intervention and target population;
- a correctly specified discrete-time hazard model.

Partial-change policies can be more plausible than setting everyone to
one exposure state, especially when a full intervention would be
unrealistic or would strain positivity.
