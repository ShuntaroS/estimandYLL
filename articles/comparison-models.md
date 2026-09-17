# Poisson and Royston-Parmar comparison models

The comparison functions standardize both exposure scenarios over the
**full analytic sample**. They do not contrast the two observed groups’
separate covariate distributions. Compare them with
`estimate_yll(target_population = "all")` to isolate differences due to
the mortality model and estimation method. A comparison with the exposed
or unexposed target also changes the population.

## Poisson model

Follow-up is split by attained age using a Lexis expansion. The model
includes exposure, a natural spline of the interval midpoint, their
interaction, baseline covariates, and a log person-time offset.
Predicted interval death counts are converted to death probabilities as
one minus their negative exponential.

``` r

data(yll_toy)
poisson <- estimate_yll_poisson(
  yll_toy[1:500, ], time_var = "period", event_var = "event",
  exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age", confounders_baseline = "sex",
  age_start = 50, age_end = 70, age_interval = 10, B = 0
)
poisson$summary
#> # A tibble: 3 × 7
#>   starting_age erl_reference erl_exposed   yll yll_se ci_low ci_high
#>          <dbl>         <dbl>       <dbl> <dbl>  <dbl>  <dbl>   <dbl>
#> 1           50         18.6        18.2  0.441     NA     NA      NA
#> 2           60          9.45        8.92 0.531     NA     NA      NA
#> 3           70          0           0    0         NA     NA      NA
```

## Royston-Parmar model

This optional method uses
[`rstpm2::stpm2()`](https://rdrr.io/pkg/rstpm2/man/gsm.html) with
attained age and delayed entry. The default has four degrees of freedom
for the baseline log cumulative hazard and proportional exposure and
covariate effects. `rp_model = "stratified"` fits separate models by
exposure group; it changes the model specification.

``` r

if (requireNamespace("rstpm2", quietly = TRUE)) {
  rp <- estimate_yll_royston_parmar(
    yll_toy[1:500, ], time_var = "period", event_var = "event",
    exposure_var = "hypertension", reference_level = "No", exposed_level = "Yes",
    age_at_entry_var = "age", confounders_baseline = "sex",
    age_start = 50, age_end = 70, age_interval = 10, B = 0
  )
  rp$summary
}
#> # A tibble: 3 × 7
#>   starting_age erl_reference erl_exposed   yll yll_se ci_low ci_high
#>          <dbl>         <dbl>       <dbl> <dbl>  <dbl>  <dbl>   <dbl>
#> 1           50         18.7        17.8  0.928     NA     NA      NA
#> 2           60          9.44        9.02 0.421     NA     NA      NA
#> 3           70          0           0    0         NA     NA      NA
```

The manuscript did not report results for this model because its
bootstrap variance was unstable in that application. Its availability
here is not an assurance of stability on other data. Inspect model
warnings, replicate estimates and `meta$bootstrap_failures` before
interpreting an interval. Finite but impossible bootstrap ERLs are also
flagged in `meta$bootstrap_unreliable`. Their raw values are retained
for diagnosis; intervals based on such a distribution are unreliable.

## Shared output and uncertainty

Both functions return the same five elements as the main estimator and
use the same direction, reference ERL minus exposed ERL. Each person’s
curve is conditioned at the starting age before averaging. Both use
normal-approximation participant bootstrap intervals; they do not expose
the main function’s `ci_method` option. Their default prediction
interval is one year and can be changed with `prediction_interval`. The
reporting interval remains `age_interval`.
