# Plot years of life lost across starting ages

Plot years of life lost across starting ages

## Usage

``` r
plot_yll(result, conf_band = TRUE)
```

## Arguments

- result:

  A result from one of the package's estimation functions.

- conf_band:

  Include available pointwise confidence intervals.

## Value

A ggplot object, editable with `+ labs()` and `+ theme()` and saveable
with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  # The saved synthetic-data example uses B = 1000 and age_interval = 5.
  result <- readRDS(system.file("extdata", "yll-example.rds", package = "estimandYLL"))
  plot_yll(result, conf_band = TRUE) +
    ggplot2::labs(title = "YLL with 95% confidence intervals")
}
```
