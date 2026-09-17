# Build the prediction grid used by the counterfactual step. For every
# subject we generate one row per integer age in [age_temp_start,
# age_temp_end]; the hazard model will then be evaluated at each of these
# rows under both exposure scenarios.
#
# We deliberately *do not* expand the full lifetime: shrinking the grid to
# [age_start, age_end] is mathematically harmless because the conditional
# survival S(y | a_start) cancels every hazard at ages below a_start, and it
# (a) cuts compute proportionally and (b) avoids extrapolating the hazard
# model below the earliest requested starting age. Requested ages can still
# lie outside the observed age support.
#
# `expo_original` preserves the observed exposure status (used downstream
# for selecting the standardization population); `expo` is set to
# the original value here and will be overwritten by the engine when it
# evaluates the two arms.
#' @noRd
yll_expand_counterfactual_data <- function(
    data,
    id_var,
    exposure_var,
    age_temp_start,
    age_temp_end
) {
  if (age_temp_end < age_temp_start) {
    stop("`age_temp_end` must be >= `age_temp_start`.", call. = FALSE)
  }
  n_ages <- as.integer(age_temp_end - age_temp_start + 1L)

  out <- data |>
    tidyr::uncount(weights = n_ages, .remove = FALSE) |>
    dplyr::group_by(.data[[id_var]]) |>
    dplyr::mutate(age_temp = age_temp_start + dplyr::row_number() - 1L) |>
    dplyr::ungroup() |>
    dplyr::rename(expo_original = dplyr::all_of(exposure_var))

  out[["expo"]] <- out$expo_original
  out |> dplyr::arrange(.data[[id_var]], age_temp)
}

# Convert per-interval hazards into per-subject survival curves under the two
# exposure scenarios.
#
# Discrete-time identity:
#   S(t) = prod_{u < t} (1 - h(u))
#
# In code: `cumprod(1 - hazard)` gives S(t+1) at row t (survival *after* the
# interval ends). To recover S(t) — survival *at the start* of the interval,
# which is what the integration step expects — we lag by one and seed the
# first interval at S = 1.
#' @noRd
yll_compute_individual_survival_curves <- function(df, id_var,
                                                   hazard0_var = "hazard0",
                                                   hazard1_var = "hazard1") {
  df |>
    arrange(.data[[id_var]], age_temp) |>
    group_by(.data[[id_var]]) |>
    mutate(
      .surv0_after_interval = cumprod(1 - .data[[hazard0_var]]),
      .surv1_after_interval = cumprod(1 - .data[[hazard1_var]]),
      surv0 = lag(.surv0_after_interval, default = 1),
      surv1 = lag(.surv1_after_interval, default = 1)
    ) |>
    ungroup()
}

# Numerically integrate a conditional survival curve to get residual life
# expectancy (LE = area under the curve).
#
# Two discretisation rules are supported:
#
#   * "left_rectangle" (default):
#       LE ≈ sum_i (a_{i+1} - a_i) * S(a_i)
#     i.e. the survival at the *start* of each interval represents the whole
#     interval. This is the natural companion to the discrete-time hazard
#     model, where S(a) is interpreted as survival up to the start of
#     interval [a, a+1).
#
#   * "trapezoidal":
#       LE ≈ sum_i (a_{i+1} - a_i) * (S(a_i) + S(a_{i+1})) / 2
#     i.e. linearly interpolate the survival curve within each interval.
#     This is slightly more accurate when the curve is smooth but mixes a
#     continuous-time intuition into a discrete-time model.
#
# A length-<=1 `age` vector is treated as zero area, which conveniently makes
# the function safe to call on degenerate cases without special-casing them
# at the call site.
#' @noRd
yll_integrate_le <- function(age, surv,
                             integration = c("left_rectangle", "trapezoidal")) {
  integration <- match.arg(integration)

  if (length(age) <= 1) {
    return(0)
  }

  ord <- order(age)
  age <- age[ord]
  surv <- surv[ord]

  d_age <- diff(age)

  if (identical(integration, "left_rectangle")) {
    return(sum(d_age * head(surv, -1)))
  }

  sum(d_age * (head(surv, -1) + tail(surv, -1)) / 2)
}
