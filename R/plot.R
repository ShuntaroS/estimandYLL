# Build the conditional-survival table the plotting helper consumes.
#
# Reads `conditional_survival_point` (the per-a_start conditional curves from
# the engine) and optionally `conditional_survival_boot` (same across bootstrap
# iterations). Both are tibbles with an `age_start` column so filtering is
# straightforward. When no bootstrap curves are available we return just the
# point estimate columns.
#' @noRd
yll_conditional_curves_with_ci <- function(conditional_survival_point,
                                           conditional_survival_boot,
                                           age_start,
                                           conf_level = 0.95) {
  if (is.null(conditional_survival_point)) {
    stop("Result has no conditional_survival_point. Re-run estimation to enable plotting.", call. = FALSE)
  }

  point_cond <- conditional_survival_point[conditional_survival_point$age_start == age_start, , drop = FALSE]
  if (nrow(point_cond) == 0) {
    stop("age_start = ", age_start, " is not present in the conditional_survival_point.", call. = FALSE)
  }

  out <- tibble(
    age_temp     = point_cond$age_temp,
    surv_cond_g0 = point_cond$surv0,
    surv_cond_g1 = point_cond$surv1
  )

  if (!is.null(conditional_survival_boot) && nrow(conditional_survival_boot) > 0) {
    boot_sub <- conditional_survival_boot[conditional_survival_boot$age_start == age_start, , drop = FALSE]
    if (nrow(boot_sub) > 0) {
      alpha <- (1 - conf_level) / 2
      ci_df <- boot_sub |>
        group_by(age_temp) |>
        summarise(
          surv_cond_g0_low  = quantile(surv0, alpha,    na.rm = TRUE),
          surv_cond_g0_high = quantile(surv0, 1 - alpha, na.rm = TRUE),
          surv_cond_g1_low  = quantile(surv1, alpha,    na.rm = TRUE),
          surv_cond_g1_high = quantile(surv1, 1 - alpha, na.rm = TRUE),
          .groups = "drop"
        )
      out <- left_join(out, ci_df, by = "age_temp")
    }
  }

  out
}

# Build the marginal-survival table for `plot_marginal_survival()`.
#
# Mirrors `yll_conditional_curves_with_ci()` but on the *unconditional*
# population-marginal curves S(t) — i.e. before conditioning on survival to
# any chosen starting age. Pointwise CI bounds at each age are again the
# lower/upper `(1 - conf_level)/2` quantiles across bootstrap iterations.
#' @noRd
yll_marginal_curves_with_ci <- function(marginal_survival_point,
                                        marginal_survival_boot,
                                        conf_level = 0.95) {
  if (is.null(marginal_survival_point)) {
    stop("Result has no marginal_survival_point. Re-run estimation to enable plotting.", call. = FALSE)
  }

  out <- tibble(
    age_temp = marginal_survival_point$age_temp,
    surv0    = marginal_survival_point$surv0,
    surv1    = marginal_survival_point$surv1
  )

  if (!is.null(marginal_survival_boot) && nrow(marginal_survival_boot) > 0) {
    alpha <- (1 - conf_level) / 2
    ci_df <- marginal_survival_boot |>
      group_by(age_temp) |>
      summarise(
        surv0_low  = quantile(surv0, alpha,    na.rm = TRUE),
        surv0_high = quantile(surv0, 1 - alpha, na.rm = TRUE),
        surv1_low  = quantile(surv1, alpha,    na.rm = TRUE),
        surv1_high = quantile(surv1, 1 - alpha, na.rm = TRUE),
        .groups = "drop"
      )
    out <- left_join(out, ci_df, by = "age_temp")
  }

  out
}

#' Plot the conditional survival probability from a starting age
#'
#' Draws the conditional survival curves \eqn{S(y \mid a_{\text{start}})} for the
#' two intervention arms, optionally with a pointwise bootstrap confidence
#' band.
#'
#' @param res A result object returned by [estimand_yll()].
#' @param age_start Numeric. The starting age \eqn{a_{\text{start}}} used to
#'   condition the survival curve.
#' @param conf_band Logical. If `TRUE` (default) and bootstrap curves are
#'   available, draws a pointwise confidence band around each curve.
#' @param reference_label,exposed_label Character labels for the two
#'   intervention arms shown in the legend.
#' @param conf_level Confidence level for the pointwise band. Defaults to the
#'   value stored in `res$meta$conf_level`, falling back to `0.95`.
#'
#' @return A `ggplot` object.
#' @export
plot_conditional_survival <- function(res,
                                      age_start,
                                      conf_band = TRUE,
                                      reference_label = "Reference",
                                      exposed_label   = "Exposed",
                                      conf_level      = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_conditional_survival().", call. = FALSE)
  }

  if (is.null(conf_level)) {
    conf_level <- if (!is.null(res$meta$conf_level)) res$meta$conf_level else 0.95
  }

  cs <- yll_conditional_curves_with_ci(
    conditional_survival_point = res$conditional_survival_point,
    conditional_survival_boot  = res$conditional_survival_boot,
    age_start                  = age_start,
    conf_level                 = conf_level
  )

  has_ci <- all(c("surv_cond_g0_low", "surv_cond_g0_high",
                  "surv_cond_g1_low", "surv_cond_g1_high") %in% names(cs))

  long <- bind_rows(
    tibble(
      age_temp = cs$age_temp,
      surv     = cs$surv_cond_g0,
      ci_low   = if (has_ci) cs$surv_cond_g0_low  else NA_real_,
      ci_high  = if (has_ci) cs$surv_cond_g0_high else NA_real_,
      arm      = reference_label
    ),
    tibble(
      age_temp = cs$age_temp,
      surv     = cs$surv_cond_g1,
      ci_low   = if (has_ci) cs$surv_cond_g1_low  else NA_real_,
      ci_high  = if (has_ci) cs$surv_cond_g1_high else NA_real_,
      arm      = exposed_label
    )
  )
  long$arm <- factor(long$arm, levels = c(reference_label, exposed_label))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = age_temp, y = surv,
                                          colour = arm, fill = arm)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x      = "Age (years)",
      y      = sprintf("Conditional survival probability (from age %g)", age_start),
      colour = "Intervention",
      fill   = "Intervention"
    ) +
    ggplot2::theme_minimal()

  if (conf_band && has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}

#' Plot the selected YLL estimate across starting ages
#'
#' Draws the main `estimate` column returned by [estimand_yll()]. When
#' `measure = "life_year_change"`, positive values mean longer life expectancy
#' after the stated intervention and negative values mean shorter life
#' expectancy.
#'
#' @param res A result object returned by [estimand_yll()].
#' @param conf_band Logical. If `TRUE` (default) and `ci_low`/`ci_high` are
#'   present in `res$summary`, draws a confidence band.
#'
#' @return A `ggplot` object.
#' @export
plot_yll_estimate <- function(res, conf_band = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_yll_estimate().", call. = FALSE)
  }

  s <- res$summary
  if (!("estimate" %in% names(s))) {
    stop("Result has no `estimate` column. Use plot_yll() for legacy result objects.", call. = FALSE)
  }

  # English: `estimate` changes meaning with `measure`, so the y-axis label
  # must be read from metadata rather than hard-coded as YLL.
  # 日本語: estimate列はmeasureによって意味が変わるため、軸ラベルもmetaから決める。
  measure <- if (!is.null(res$meta$measure)) res$meta$measure else unique(s$measure)[[1]]
  y_label <- if (identical(measure, "life_year_change")) {
    "Life-year change (years)"
  } else {
    "Years of life lost (years)"
  }

  # English: The zero line is clinically useful: values above zero mean
  # longer life expectancy for life-year change, while values below zero mean
  # harm. For YLL, zero is the null exposure contrast.
  # 日本語: 0線は解釈の基準。life_year_changeでは正が余命増加、負が余命減少を表す。
  p <- ggplot2::ggplot(s, ggplot2::aes(x = starting_age, y = estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Starting age (years)",
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (conf_band && all(c("ci_low", "ci_high") %in% names(s)) &&
      any(is.finite(s$ci_low))) {
    # English: CI columns have already been oriented to the chosen measure in
    # `yll_add_measure_columns()`, so plotting can use them directly.
    # 日本語: 信頼区間はmeasureに合わせて符号調整済みなので、そのまま描画する。
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}

#' Plot YLL across starting ages
#'
#' Draws estimated years of life lost (YLL) as a function of the starting age,
#' with an optional confidence band when bootstrap CIs are available in the
#' result object.
#'
#' @param res A result object returned by [estimand_yll()].
#' @param conf_band Logical. If `TRUE` (default) and `ci_low`/`ci_high` are
#'   present in `res$summary`, draws a confidence band.
#'
#' @return A `ggplot` object.
#' @export
plot_yll <- function(res, conf_band = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_yll().", call. = FALSE)
  }

  s <- res$summary
  # English: Keep `plot_yll()` backward compatible. New results have
  # `estimate`; old internal/legacy results only have `yll`.
  # 日本語: 新APIではestimateを描き、旧結果では従来通りyllを描く。
  value_var <- if ("estimate" %in% names(s)) "estimate" else "yll"
  y_label <- if ("estimate" %in% names(s) && identical(res$meta$measure, "life_year_change")) {
    "Life-year change (years)"
  } else {
    "Years of life lost (years)"
  }
  s$.plot_value <- s[[value_var]]

  p <- ggplot2::ggplot(s, ggplot2::aes(x = starting_age, y = .plot_value)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Starting age (years)",
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (conf_band && all(c("ci_low", "ci_high") %in% names(s)) &&
      any(is.finite(s$ci_low))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}

#' Plot marginal survival curves for the two intervention arms
#'
#' Draws the population-level marginal survival curves \eqn{S(t)} for the two
#' intervention arms, optionally with pointwise bootstrap confidence bands.
#'
#' @inheritParams plot_conditional_survival
#'
#' @return A `ggplot` object.
#' @export
plot_marginal_survival <- function(res,
                                   conf_band = TRUE,
                                   reference_label = "Reference",
                                   exposed_label   = "Exposed",
                                   conf_level      = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_marginal_survival().", call. = FALSE)
  }

  if (is.null(conf_level)) {
    conf_level <- if (!is.null(res$meta$conf_level)) res$meta$conf_level else 0.95
  }

  ms <- yll_marginal_curves_with_ci(
    marginal_survival_point = res$marginal_survival_point,
    marginal_survival_boot  = res$marginal_survival_boot,
    conf_level              = conf_level
  )

  has_ci <- all(c("surv0_low", "surv0_high", "surv1_low", "surv1_high") %in% names(ms))

  long <- bind_rows(
    tibble(
      age_temp = ms$age_temp,
      surv     = ms$surv0,
      ci_low   = if (has_ci) ms$surv0_low  else NA_real_,
      ci_high  = if (has_ci) ms$surv0_high else NA_real_,
      arm      = reference_label
    ),
    tibble(
      age_temp = ms$age_temp,
      surv     = ms$surv1,
      ci_low   = if (has_ci) ms$surv1_low  else NA_real_,
      ci_high  = if (has_ci) ms$surv1_high else NA_real_,
      arm      = exposed_label
    )
  )
  long$arm <- factor(long$arm, levels = c(reference_label, exposed_label))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = age_temp, y = surv,
                                          colour = arm, fill = arm)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x      = "Age (years)",
      y      = "Marginal survival probability",
      colour = "Intervention",
      fill   = "Intervention"
    ) +
    ggplot2::theme_minimal()

  if (conf_band && has_ci) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_low, ymax = ci_high),
      alpha = 0.2, colour = NA
    )
  }

  p
}
