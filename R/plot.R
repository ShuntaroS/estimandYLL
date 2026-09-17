#' Plot years of life lost across starting ages
#'
#' @param result A result from one of the package's estimation functions.
#' @param conf_band Include available pointwise confidence intervals.
#' @return A ggplot object, editable with `+ labs()` and `+ theme()` and
#'   saveable with `ggplot2::ggsave()`.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   # The saved synthetic-data example uses B = 1000 and age_interval = 5.
#'   result <- readRDS(system.file("extdata", "yll-example.rds", package = "estimandYLL"))
#'   plot_yll(result, conf_band = TRUE) +
#'     ggplot2::labs(title = "YLL with 95% confidence intervals")
#' }
#' @export
plot_yll <- function(result, conf_band = TRUE) {
  yll_check_plot_input(result, conf_band)
  estimates <- result$summary
  plot <- ggplot2::ggplot(estimates, ggplot2::aes(x = .data$starting_age, y = .data$yll)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linetype = "dashed") +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Starting age (years)", y = "Years of life lost (years)") +
    ggplot2::theme_minimal()
  if (nrow(estimates) > 1L) plot <- plot + ggplot2::geom_line()
  if (conf_band && any(is.finite(estimates$ci_low))) {
    plot <- plot + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high), width = 0.5, na.rm = TRUE
    )
  }
  plot
}

#' Plot standardized survival from a specified starting age
#'
#' Each curve is the mean of individual survival curves restarted at one at
#' `age_start`. Both scenarios use the same target population.
#' @inheritParams plot_yll
#' @param age_start A starting age present in `result$survival_curves`.
#' @param reference_label,exposed_label Optional legend labels. By default,
#'   the exposure values recorded in the result are displayed.
#' @return An editable ggplot object. Confidence bands are pointwise bootstrap
#'   intervals using the result's confidence level and interval method, clipped
#'   to zero and one. They are not simultaneous bands.
#' @export
plot_survival <- function(result, age_start, conf_band = TRUE,
                           reference_label = NULL, exposed_label = NULL) {
  yll_check_plot_input(result, conf_band)
  yll_check_number(age_start, "age_start")
  curves <- result$survival_curves
  curves <- curves[curves$starting_age == age_start, , drop = FALSE]
  if (!nrow(curves)) stop("`age_start` is not a reported starting age in this result.", call. = FALSE)
  if (is.null(reference_label)) reference_label <- paste0("Reference: ", result$meta$reference_level)
  if (is.null(exposed_label)) exposed_label <- paste0("Exposed: ", result$meta$exposed_level)
  labels <- c(reference_label, exposed_label)
  columns <- c("survival_reference", "survival_exposed")
  bootstrap <- result$bootstrap_survival_curves
  bootstrap <- bootstrap[bootstrap$starting_age == age_start, , drop = FALSE]
  alpha <- (1 - result$meta$conf_level) / 2
  plot_data <- vector("list", 2L)
  for (scenario in seq_along(columns)) {
    column <- columns[scenario]
    scenario_data <- tibble::tibble(
      age = curves$age, survival = curves[[column]],
      scenario = labels[scenario], ci_low = NA_real_, ci_high = NA_real_
    )
    if (conf_band && nrow(bootstrap)) {
      for (row in seq_len(nrow(scenario_data))) {
        values <- bootstrap[[column]][bootstrap$age == scenario_data$age[row]]
        if (length(values) < 2L) next
        if (result$meta$ci_method == "normal") {
          margin <- stats::qnorm(1 - alpha) * stats::sd(values)
          bounds <- scenario_data$survival[row] + c(-margin, margin)
        } else {
          bounds <- stats::quantile(values, c(alpha, 1 - alpha), names = FALSE)
        }
        scenario_data$ci_low[row] <- max(0, bounds[1])
        scenario_data$ci_high[row] <- min(1, bounds[2])
      }
    }
    plot_data[[scenario]] <- scenario_data
  }
  plot_data <- bind_rows(plot_data)
  plot_data$scenario <- factor(plot_data$scenario, levels = labels)
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$age, y = .data$survival, colour = .data$scenario, fill = .data$scenario
  )) + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Age (years)", y = paste0("Survival probability from age ", age_start),
                  colour = "Exposure scenario", fill = "Exposure scenario") +
    ggplot2::theme_minimal()
  if (conf_band && any(is.finite(plot_data$ci_low))) {
    plot <- plot + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$ci_low, ymax = .data$ci_high),
      alpha = 0.2, colour = NA, na.rm = TRUE
    )
  }
  if (nrow(curves) > 1L) plot <- plot + ggplot2::geom_line(linewidth = 0.8)
  else plot <- plot + ggplot2::geom_point()
  plot
}

yll_check_plot_input <- function(result, conf_band) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install 'ggplot2' to use the plotting functions.", call. = FALSE)
  }
  if (!is.list(result) || !all(c("summary", "survival_curves", "meta") %in% names(result))) {
    stop("`result` must be returned by an estimandYLL estimation function.", call. = FALSE)
  }
  if (!is.logical(conf_band) || length(conf_band) != 1L || is.na(conf_band)) {
    stop("`conf_band` must be TRUE or FALSE.", call. = FALSE)
  }
}
