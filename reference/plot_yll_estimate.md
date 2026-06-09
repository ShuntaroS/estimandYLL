# Plot the selected YLL estimate across starting ages

Draws the main `estimate` column. For `measure = "yll"`, the y-axis is
years of life lost. For `measure = "life_year_change"`, the y-axis is
change in life expectancy.

## Usage

``` r
plot_yll_estimate(res, conf_band = TRUE)
```

## Arguments

- res:

  A result object returned by
  [`estimand_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimand_yll.md).

- conf_band:

  Logical. If `TRUE`, draw a confidence band when available.

## Value

A `ggplot` object.
