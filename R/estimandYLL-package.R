#' estimandYLL: Years of Life Lost via the G-Formula
#'
#' Estimate differences in restricted expected residual lifetime under two
#' binary exposure scenarios, standardized to a stated baseline population.
#' See [estimate_yll()] for the main estimator and interpretation.
#'
#' @keywords internal
#' @importFrom dplyr arrange bind_rows group_by left_join mutate select summarise ungroup lag
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr uncount
#' @importFrom purrr map_dfr
#' @importFrom rlang .data :=
#' @importFrom methods selectMethod
#' @importFrom stats as.formula ave binomial glm predict poisson qnorm quantile sd setNames
#' @importFrom utils head tail txtProgressBar setTxtProgressBar
#' @importFrom survival Surv survSplit
#' @importFrom Epi Ns
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor with_progress
"_PACKAGE"

# Suppress R CMD check NOTEs about non-standard evaluation in dplyr/ggplot2.
utils::globalVariables(c(
  ".data",
  ".surv0_after_interval", ".surv1_after_interval",
  ".surv_after_interval", ".surv_at_age",
  "age_temp", "age_start", "starting_age",
  ".age_in", ".age_out", ".age_mid", ".death", ".event", ".expo",
  ".prediction_id", ".prediction_row",
  "tstart",
  "expo", "expo_original",
  "surv", "surv0", "surv1", "surv0_ind", "surv1_ind",
  "hazard0", "hazard1", "lex.dur", "lex.Xst",
  "ci_low", "ci_high", "ci_method",

  "le_reference", "le_exposed",
  "yll", "le_m0", "le_m1",
  "se_yll", "se_le_m0", "se_le_m1",
  "arm"
))
