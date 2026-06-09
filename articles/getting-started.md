# Getting started with estimandYLL

`estimandYLL` estimates **years of life lost (YLL)** and counterfactual
life expectancy under the parametric g-formula. The main interface is
[`estimand_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimand_yll.md).

The estimand is stated through three arguments:

1.  `target_population`: who to average over.
2.  `intervention`: which exposure worlds or policy change to compare.
3.  `measure`: how to display the result.

``` r

library(estimandYLL)
data(yll_toy)
str(yll_toy)
#> 'data.frame':    5000 obs. of  9 variables:
#>  $ id             : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ period         : num  15 15 15 15 15 ...
#>  $ event          : int  0 0 1 0 0 0 0 0 0 1 ...
#>  $ smoke_binary   : Factor w/ 2 levels "Never","Current/Ever": 1 1 1 2 2 2 2 1 1 2 ...
#>  $ hypertension   : Factor w/ 2 levels "No","Yes": 1 2 1 1 1 1 1 1 2 2 ...
#>  $ sex            : Factor w/ 2 levels "Female","Male": 2 2 2 1 1 2 1 2 1 1 ...
#>  $ education_years: num  14 15 16 11 12 15 18 14 19 14 ...
#>  $ bmi            : num  22 23 22 24.7 26.6 27.6 26.3 25.6 22.7 18.7 ...
#>  $ age            : num  54 49 64 58 47 61 51 57 43 65 ...
```

## Full contrasts

`intervention = "full_contrast"` compares the all-reference and
all-exposed counterfactual worlds. The target population determines
whether the result is ATE-like, ATT-like, or ATC-like.

``` r

res_all <- estimand_yll(
  data = yll_toy,
  target_population = "all",
  intervention = "full_contrast",
  measure = "yll",
  B = 50,
  method = "normal",
  show_progress = FALSE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  use_future = FALSE
)

res_all$summary
```

Observed-exposed and observed-unexposed target populations are specified
directly:

``` r

res_exposed <- estimand_yll(
  data = yll_toy,
  target_population = "exposed",
  intervention = "full_contrast",
  measure = "yll",
  B = 50, method = "normal", show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)

res_unexposed <- estimand_yll(
  data = yll_toy,
  target_population = "unexposed",
  intervention = "full_contrast",
  measure = "life_year_change",
  B = 50, method = "normal", show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)
```

With `measure = "life_year_change"`, an ATC-like harmful change is
reported as a negative value.

## Partial-change policies

`intervention = "partial_change"` compares the natural observed exposure
distribution with a policy-like change in which only a specified
fraction changes exposure state. In statistical terminology this is a
stochastic intervention, but the API uses clinical policy language.

Example: 30% of observed-exposed individuals move to the reference
level, averaged over observed-exposed individuals.

``` r

res_quit <- estimand_yll(
  data = yll_toy,
  target_population = "exposed",
  intervention = "partial_change",
  change_from = "exposed",
  change_to = "reference",
  change_probability = 0.30,
  measure = "life_year_change",
  B = 50, show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "smoke_binary",
  reference_level = "Never", exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi")
)
```

## Visualisation

``` r

plot_yll_estimate(res_all)
plot_conditional_survival(res_all, age_start = 60)
plot_marginal_survival(res_all)
```

## Result object

`summary` contains one row per starting age and includes the selected
`estimate`, confidence interval, `target_population`, `intervention`,
and `measure`. `detailed_results` additionally keeps the original `yll`,
life-year change, and arm-specific life expectancies for checking.
