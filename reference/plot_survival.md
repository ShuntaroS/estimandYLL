# Plot standardized survival from a specified starting age

Each curve is the mean of individual survival curves restarted at one at
`age_start`. Both scenarios use the same target population.

## Usage

``` r
plot_survival(
  result,
  age_start,
  conf_band = TRUE,
  reference_label = NULL,
  exposed_label = NULL
)
```

## Arguments

- result:

  A result from one of the package's estimation functions.

- age_start:

  A starting age present in `result$survival_curves`.

- conf_band:

  Include available pointwise confidence intervals.

- reference_label, exposed_label:

  Optional legend labels. By default, the exposure values recorded in
  the result are displayed.

## Value

An editable ggplot object. Confidence bands are pointwise bootstrap
intervals using the result's confidence level and interval method,
clipped to zero and one. They are not simultaneous bands.
