#' estimandYLL: Years of Life Lost via the G-Formula
#'
#' Estimates Years of Life Lost (YLL) and counterfactual life expectancy
#' under user-specified interventions on a binary exposure using the
#' parametric g-formula on the attained-age time scale. The main interface,
#' [estimand_yll()], makes the target population, intervention, and result
#' measure explicit.
#'
#' @keywords internal
#' @importFrom dplyr arrange bind_rows filter group_by inner_join left_join
#'   mutate rename rename_with select summarise transmute ungroup lag
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
#' @importFrom progressr handlers progressor with_progress
"_PACKAGE"

# Suppress R CMD check NOTEs about non-standard evaluation in dplyr/ggplot2.
utils::globalVariables(c(
  ".data",
  ".surv0_after_interval", ".surv1_after_interval",
  ".surv_after_interval", ".surv_at_age",
  ".p_exposed", ".surv_mix",
  "age_temp", "age_start", "starting_age",
  ".age_in", ".age_out", ".age_mid", ".death", ".event", ".expo",
  ".prediction_id", ".prediction_row",
  "tstart",
  "expo", "expo_original",
  "surv", "surv0", "surv1", "surv0_ind", "surv1_ind",
  "hazard0", "hazard1", "lex.dur", "lex.Xst",
  "surv_cond_g0", "surv_cond_g1",
  "ci_low", "ci_high", "ci_method",
  "estimate", "measure", "target_population", "intervention",
  "life_year_change", "le_before", "le_after", ".plot_value",
  "le_reference", "le_exposed",
  "yll", "le_m0", "le_m1",
  "se_yll", "se_le_m0", "se_le_m1",
  "arm"
))
