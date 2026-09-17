# Run after generate-example.R, from the repository root.
# Rscript notes/published-example/publish-figures.R
library(estimandYLL)

result <- readRDS("local/published-example/result.rds")

# Check the actual fitted result before copying it into the public documentation.
stopifnot(
  result$meta$B == 1000,
  result$meta$n == 5000,
  result$meta$target_population == "unexposed",
  result$meta$age_interval == 5,
  result$meta$conf_level == 0.95,
  result$meta$ci_method == "normal",
  result$meta$use_future,
  result$meta$bootstrap_successful == 1000,
  nrow(result$meta$bootstrap_failures) == 0,
  nrow(result$meta$bootstrap_unreliable) == 0,
  identical(result$summary$starting_age, seq(40, 90, by = 5)),
  length(unique(result$bootstrap_estimates$iteration)) == 1000,
  all(is.finite(result$summary$ci_low)),
  all(is.finite(result$summary$ci_high))
)

yll_plot <- plot_yll(result, conf_band = TRUE) +
  ggplot2::labs(title = "YLL among participants without hypertension")
survival_plot <- plot_survival(result, age_start = 50) + ggplot2::theme_bw()

# Verify that the rendered intervals use the saved numerical limits.
plot_layers <- ggplot2::ggplot_build(yll_plot)$data
error_bars <- plot_layers[[length(plot_layers)]]
stopifnot(
  isTRUE(all.equal(error_bars$ymin, result$summary$ci_low)),
  isTRUE(all.equal(error_bars$ymax, result$summary$ci_high))
)

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
file.copy("local/published-example/result.rds", "inst/extdata/yll-example.rds",
          overwrite = TRUE)
ggplot2::ggsave("man/figures/yll-example.png", yll_plot,
               width = 6, height = 4, dpi = 180)
ggplot2::ggsave("man/figures/survival-example.png", survival_plot,
               width = 6, height = 4, dpi = 180)
message("Saved the complete result and both CI figures.")
