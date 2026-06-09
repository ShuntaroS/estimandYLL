# Plot YLL across starting ages

Draws the selected estimate as a function of the starting age. For new
[`estimand_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimand_yll.md)
results this uses `res$summary$estimate`; for legacy result objects it
falls back to `res$summary$yll`.

## Usage

``` r
plot_yll(res, conf_band = TRUE)
```

## Arguments

- res:

  A result object returned by
  [`estimand_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimand_yll.md).

- conf_band:

  Logical. If `TRUE` (default) and `ci_low`/`ci_high` are present in
  `res$summary`, draws a confidence band.

## Value

A `ggplot` object.
