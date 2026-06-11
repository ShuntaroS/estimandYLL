#' Estimate years of life lost from a clearly stated estimand
#'
#' `estimand_yll()` is the main user-facing interface. It asks users to state
#' three pieces of the estimand: the target population, the intervention, and
#' the measure used to display the result.
#'
#' @param data A `data.frame` with one row per individual, containing the
#'   variables specified by `id_var`, `time_var`, `event_var`, `exposure_var`,
#'   `age_at_entry_var`, and any `confounders_baseline`.
#' @param B Number of bootstrap iterations. Set to `0` to skip bootstrapping.
#' @param seed Integer seed for reproducibility.
#' @param conf_level Confidence level for the intervals.
#' @param method Either `"normal"` or `"percentile"`.
#' @param show_progress Show a progress bar during bootstrap.
#' @param id_var,time_var,event_var,exposure_var Column names in `data`.
#' @param reference_level,exposed_level Values labelling the two exposure
#'   levels. If `NULL` they are inferred from the data.
#' @param age_at_entry_var Column name for age at study entry.
#' @param age_start,age_end,age_interval Starting ages for which results are
#'   reported, given as `seq(age_start, age_end, age_interval)`.
#' @param confounders_baseline Character vector of baseline confounder column
#'   names.
#' @param integration Discretisation rule for the area under the conditional
#'   survival curve.
#' @param use_future If `TRUE`, parallelise the bootstrap with
#'   [future.apply::future_lapply()].
#' @param target_population Character. Who should the counterfactual life
#'   expectancies be averaged over? Use `"all"`, `"exposed"`, or
#'   `"unexposed"`. A custom function `f(data)` returning a logical vector may
#'   also be supplied for advanced use.
#' @param intervention Character. `"full_contrast"` compares the all-reference
#'   and all-exposed worlds. `"partial_change"` compares the natural observed
#'   exposure distribution with a policy-like partial change.
#' @param change_from,change_to For `intervention = "partial_change"`, the
#'   exposure state that some people move from and to. Use `"reference"` or
#'   `"exposed"` (the actual values of `reference_level` / `exposed_level` are
#'   also accepted).
#' @param change_probability For `intervention = "partial_change"`, the
#'   probability that a member of `change_from` moves to `change_to`.
#' @param measure Character. `"yll"` reports the usual `LE_reference -
#'   LE_exposed`. `"life_year_change"` reports the clinically oriented change
#'   in life expectancy for the stated intervention direction.
#'
#' @return A list with `summary`, `detailed_results`, `meta`, and survival
#'   curves. New API result columns include `estimate`, `measure`,
#'   `target_population`, `intervention`, `le_before`, `le_after`, and
#'   `life_year_change`.
#'
#' @export
estimand_yll <- function(
    data,
    id_var = "id",
    time_var,
    event_var,
    exposure_var,
    reference_level = NULL,
    exposed_level = NULL,
    age_at_entry_var,
    target_population = c("all", "exposed", "unexposed"),
    intervention = c("full_contrast", "partial_change"),
    change_from = NULL,
    change_to = NULL,
    change_probability = NULL,
    measure = c("yll", "life_year_change"),
    age_start = 50,
    age_end = 100,
    age_interval = 5,
    confounders_baseline = NULL,
    B = 1000,
    seed = 1,
    conf_level = 0.95,
    method = c("normal", "percentile"),
    integration = c("left_rectangle", "trapezoidal"),
    show_progress = TRUE,
    use_future = TRUE
) {
  # English: Match early so later helper functions can assume a single,
  # validated choice instead of repeatedly handling character vectors.
  # 日本語: ここで引数を1つに確定し、以降の補助関数では検証済みの値だけを扱う。
  intervention <- match.arg(intervention)
  measure <- match.arg(measure)
  method <- match.arg(method)
  integration <- match.arg(integration)

  # English: Fail before fitting the hazard model. Missing column errors from
  # inside model formulas are hard for clinical users to interpret.
  # 日本語: モデル式の中で落ちると原因が分かりにくいため、必要列は最初に確認する。
  yll_check_required_columns(
    data,
    c(id_var, time_var, event_var, exposure_var, age_at_entry_var, confounders_baseline)
  )
  yll_warn_entry_age_in_confounders(age_at_entry_var, confounders_baseline)

  # English: Resolve the binary exposure labels once. The rest of the package
  # works with the abstract pair "reference" and "exposed" so that the public
  # API can accept clinical labels such as "No"/"Yes" or "Never"/"Current".
  # 日本語: 実データ上のラベルを一度だけ確定し、内部では reference/exposed として扱う。
  exposure_levels <- yll_resolve_binary_levels(
    data[[exposure_var]],
    reference_level = reference_level,
    exposed_level = exposed_level
  )
  reference_level <- exposure_levels$reference_level
  exposed_level <- exposure_levels$exposed_level

  # English: This is the central estimand decision: who contributes to the
  # final marginal average. It does not subset the model-fitting data.
  # 日本語: ここで「誰で平均するか」を決める。ハザードモデルの学習データは削らない。
  target_spec <- yll_resolve_target_population_argument(
    target_population = target_population,
    reference_level = reference_level,
    exposed_level = exposed_level
  )

  # English: Build the two counterfactual worlds, then let the existing engine
  # do the modelling. 日本語: 介入の指定だけをここで作り、推定本体は既存エンジンに任せる。
  intervention_spec <- yll_build_intervention_spec(
    intervention = intervention,
    change_from = change_from,
    change_to = change_to,
    change_probability = change_probability,
    reference_level = reference_level,
    exposed_level = exposed_level
  )

  # English: The old generic engine already estimates survival under two
  # intervention arms and averages within a target population. The new API is
  # therefore a translation layer from clinical estimand language to engine
  # inputs, not a second implementation of the estimator.
  # 日本語: 既存エンジンは2つの介入下の生存曲線と標的集団平均を計算できるため、
  # 新APIは臨床的な指定をエンジン入力へ翻訳する層として実装する。
  res <- estimate_yll_gformula_intervention(
    data = data,
    intervention_reference = intervention_spec$intervention_reference,
    intervention_exposed = intervention_spec$intervention_exposed,
    target_population = target_spec$rule,
    B = B,
    seed = seed,
    conf_level = conf_level,
    method = method,
    show_progress = show_progress,
    id_var = id_var,
    time_var = time_var,
    event_var = event_var,
    exposure_var = exposure_var,
    reference_level = reference_level,
    exposed_level = exposed_level,
    age_at_entry_var = age_at_entry_var,
    age_start = age_start,
    age_end = age_end,
    age_interval = age_interval,
    confounders_baseline = confounders_baseline,
    integration = integration,
    use_future = use_future
  )

  # English: The engine always returns the mathematical YLL orientation
  # (`LE_reference - LE_exposed`). This final step adds the requested display
  # measure without changing the underlying survival estimates.
  # 日本語: エンジンの基本出力は常にYLL向きなので、最後に表示用のmeasureを追加する。
  yll_apply_measure(
    res = res,
    measure = measure,
    target_population_label = target_spec$label,
    intervention_label = intervention,
    intervention_direction = intervention_spec$direction
  )
}

#' @rdname estimand_yll
#' @export
estimate_yll <- function(...) {
  estimand_yll(...)
}

# Convert a user-facing target population label into the row-selection rule
# used by the engine. 日本語: target populationはestimandの中心なので、
# "誰で平均するか" をここで明示的に解決する。
#' @noRd
yll_resolve_target_population_argument <- function(target_population,
                                                   reference_level,
                                                   exposed_level) {
  if (is.function(target_population)) {
    # English: Advanced users may define their own target population, e.g. a
    # clinical subgroup. The function is evaluated later on counterfactual
    # data containing `expo_original`.
    # 日本語: 上級者向けに任意の標的集団関数も許す。評価時にはexpo_originalが使える。
    return(list(label = "custom", rule = target_population))
  }

  target_population <- match.arg(target_population, c("all", "exposed", "unexposed"))

  if (identical(target_population, "all")) {
    # English: NULL means "do not filter" in the lower-level engine.
    # 日本語: 内部エンジンではNULLが「全体で平均」を意味する。
    return(list(label = "all", rule = NULL))
  }
  if (identical(target_population, "exposed")) {
    # English: ATT-like target population. We select by the observed exposure,
    # not by the post-intervention exposure assignment.
    # 日本語: ATT的な標的集団。介入後ではなく観察時の曝露状態で選ぶ。
    return(list(
      label = "exposed",
      rule = function(data) data$expo_original == exposed_level
    ))
  }

  # English: ATC-like target population. Again, membership is defined by the
  # observed exposure status preserved as `expo_original`.
  # 日本語: ATC的な標的集団。expo_originalに保存された観察時非曝露で選ぶ。
  list(
    label = "unexposed",
    rule = function(data) data$expo_original == reference_level
  )
}

#' @noRd
yll_normalize_change_state <- function(x, reference_level, exposed_level, arg_name) {
  # English: Partial-change policies are easier to read with words
  # ("reference"/"exposed"), but accepting the actual data labels makes the
  # function forgiving in scripts and notebooks.
  # 日本語: partial_changeではreference/exposedが読みやすいが、実データのラベルも許す。
  if (is.null(x)) {
    stop("`", arg_name, "` is required when `intervention = \"partial_change\"`.", call. = FALSE)
  }

  x_chr <- as.character(x)
  if (identical(x_chr, "reference") || identical(x_chr, as.character(reference_level))) {
    return("reference")
  }
  if (identical(x_chr, "exposed") || identical(x_chr, as.character(exposed_level))) {
    return("exposed")
  }

  stop(
    "`", arg_name, "` must be \"reference\", \"exposed\", or one of the observed exposure levels.",
    call. = FALSE
  )
}

#' @noRd
yll_check_change_probability <- function(change_probability) {
  # English: This is a probability of changing state among people currently in
  # `change_from`; it is not a prevalence target. Values must stay in [0, 1].
  # 日本語: これはchange_from集団のうち状態が変わる確率であり、介入後有病率ではない。
  if (is.null(change_probability) || length(change_probability) != 1L ||
      is.na(change_probability) || change_probability < 0 || change_probability > 1) {
    stop("`change_probability` must be a single number between 0 and 1.", call. = FALSE)
  }
  as.numeric(change_probability)
}

# Build the intervention pair passed to the generic engine.
# English: "partial_change" replaces the old stochastic API with a clinical
# policy language. 日本語: stochasticという言葉を表に出さず、一部が変わる介入として指定する。
#' @noRd
yll_build_intervention_spec <- function(intervention,
                                        change_from,
                                        change_to,
                                        change_probability,
                                        reference_level,
                                        exposed_level) {
  if (identical(intervention, "full_contrast")) {
    # English: Deterministic contrast. The engine interprets exposure-level
    # labels as degenerate interventions: everyone assigned reference vs
    # everyone assigned exposed.
    # 日本語: 決定論的な比較。曝露ラベルを渡すと「全員reference」「全員exposed」
    # として内部で解釈される。
    return(list(
      intervention_reference = reference_level,
      intervention_exposed = exposed_level,
      direction = "full_contrast"
    ))
  }

  from <- yll_normalize_change_state(change_from, reference_level, exposed_level, "change_from")
  to <- yll_normalize_change_state(change_to, reference_level, exposed_level, "change_to")
  if (identical(from, to)) {
    stop("`change_from` and `change_to` must be different.", call. = FALSE)
  }
  p_change <- yll_check_change_probability(change_probability)

  # English: The reference arm for a partial-change policy is the natural
  # observed exposure distribution. This is encoded as P(exposed|observed
  # reference)=0 and P(exposed|observed exposed)=1.
  # 日本語: partial_changeの比較元は自然経過。観察時非曝露は非曝露、観察時曝露は曝露
  # のままという退化した確率介入として表す。
  natural <- yll_make_binary_stochastic_intervention(
    prob_exposed_if_unexposed = 0,
    prob_exposed_if_exposed = 1,
    reference_level = reference_level,
    exposed_level = exposed_level
  )

  if (identical(from, "reference") && identical(to, "exposed")) {
    # English: Example: 20% of observed-unexposed people become exposed; those
    # already exposed remain exposed. This is ATC-like when averaged over
    # `target_population = "unexposed"`.
    # 日本語: 例: 非曝露者の20%が曝露になる。target_population="unexposed"ならATC的。
    changed <- yll_make_binary_stochastic_intervention(
      prob_exposed_if_unexposed = p_change,
      prob_exposed_if_exposed = 1,
      reference_level = reference_level,
      exposed_level = exposed_level
    )
    direction <- "reference_to_exposed"
  } else {
    # English: Example: 30% of observed-exposed people move to reference; those
    # already unexposed remain unexposed. This is ATT-like when averaged over
    # `target_population = "exposed"`.
    # 日本語: 例: 曝露者の30%が非曝露になる。target_population="exposed"ならATT的。
    changed <- yll_make_binary_stochastic_intervention(
      prob_exposed_if_unexposed = 0,
      prob_exposed_if_exposed = 1 - p_change,
      reference_level = reference_level,
      exposed_level = exposed_level
    )
    direction <- "exposed_to_reference"
  }

  list(
    intervention_reference = natural,
    intervention_exposed = changed,
    direction = direction
  )
}

# Add the clinically oriented measure to the result object while preserving
# the original YLL columns. 日本語: 既存のYLL列は残し、主結果としてestimateを追加する。
#' @noRd
yll_apply_measure <- function(res,
                              measure,
                              target_population_label,
                              intervention_label,
                              intervention_direction) {
  # English: The multiplier determines whether `life_year_change` is equal to
  # YLL or -YLL. We keep this as a separate helper so the sign convention is
  # visible and testable.
  # 日本語: life_year_changeがYLLと同じ向きか反対向きかをここで決める。
  # 符号規則を独立関数にして、読みやすくテストしやすくする。
  multiplier <- yll_life_year_change_multiplier(
    target_population_label = target_population_label,
    intervention_label = intervention_label,
    intervention_direction = intervention_direction
  )

  # English: Add the same measure columns to both detailed and summary tables.
  # `detailed_results` keeps every diagnostic column; `summary` is what users
  # usually print in manuscripts or slides.
  # 日本語: detailed_resultsとsummaryの両方にmeasure列を付ける。論文表ではsummaryが主役。
  res$detailed_results <- yll_add_measure_columns(
    res$detailed_results,
    measure = measure,
    target_population_label = target_population_label,
    intervention_label = intervention_label,
    life_year_change_multiplier = multiplier
  )

  res$summary <- yll_add_measure_columns(
    res$summary,
    measure = measure,
    target_population_label = target_population_label,
    intervention_label = intervention_label,
    life_year_change_multiplier = multiplier
  )

  # English: Store the estimand labels in metadata so plots and downstream
  # tables can report what was estimated without reconstructing the call.
  # 日本語: 推定したtarget/intervention/measureをメタデータに残し、図表で再利用できるようにする。
  res$meta$target_population <- target_population_label
  res$meta$intervention <- intervention_label
  res$meta$intervention_direction <- intervention_direction
  res$meta$measure <- measure

  class(res) <- unique(c("yll_estimate", class(res)))
  res
}

#' @noRd
yll_life_year_change_multiplier <- function(target_population_label,
                                            intervention_label,
                                            intervention_direction) {
  if (identical(intervention_label, "partial_change")) {
    # English: Partial-change results are returned by the engine as
    # natural - changed. Clinically, "life-year change" should read as
    # changed - natural, so we flip the sign.
    # 日本語: エンジンは「自然経過 - 介入後」で返すが、臨床的な変化量は
    # 「介入後 - 自然経過」なので符号を反転する。
    return(-1)
  }

  # For a full contrast, the clinical direction depends on who is targeted:
  # exposed people are imagined moving to reference; unexposed/all are read as
  # moving toward exposed. 日本語: ATCでは「非曝露者が曝露になる」ため符号を反転する。
  if (identical(target_population_label, "exposed")) {
    return(1)
  }
  -1
}

#' @noRd
yll_add_measure_columns <- function(df,
                                    measure,
                                    target_population_label,
                                    intervention_label,
                                    life_year_change_multiplier) {
  if (!("yll" %in% names(df))) {
    # English: Some malformed or future result tables may not carry YLL.
    # Returning unchanged keeps this helper defensive and avoids a secondary
    # error that would hide the original issue.
    # 日本語: yll列が無い表では余計なエラーを出さず、そのまま返す。
    return(df)
  }

  # English: Always keep both orientations. `estimate` is the user-selected
  # primary column; `yll` and `life_year_change` remain available for checking.
  # The measure choice must be evaluated *before* the mutate: inside mutate
  # the name `measure` refers to the just-created column (one value per row),
  # not to the function argument.
  # 日本語: 選択された主結果はestimateだが、yllとlife_year_changeの両方を残す。
  # mutate内では`measure`が新設の列を指してしまうため、判定は先に済ませる。
  measure_label <- measure
  measure_is_yll <- identical(measure, "yll")
  out <- df |>
    mutate(
      target_population = target_population_label,
      intervention = intervention_label,
      measure = measure_label,
      life_year_change = life_year_change_multiplier * .data[["yll"]],
      estimate = if (measure_is_yll) .data[["yll"]] else .data[["life_year_change"]]
    )

  if (all(c("le_reference", "le_exposed") %in% names(out))) {
    # English: Translate model-arm labels into before/after labels. This is
    # especially helpful for clinical readers of partial-change policies.
    # 日本語: モデル上のreference/exposedを、臨床的に読みやすいbefore/afterへ写す。
    out <- yll_add_before_after_le(out, life_year_change_multiplier)
  }

  if (all(c("ci_low", "ci_high") %in% names(out)) &&
      identical(measure, "life_year_change") &&
      identical(life_year_change_multiplier, -1)) {
    # English: When multiplying an interval by -1, lower and upper bounds swap.
    # This is easy to miss and would otherwise produce inverted CIs.
    # 日本語: -1を掛けると信頼区間の下限・上限が入れ替わるため、必ず順序も反転する。
    old_low <- out$ci_low
    old_high <- out$ci_high
    out$ci_low <- -old_high
    out$ci_high <- -old_low
  }

  out
}

#' @noRd
yll_add_before_after_le <- function(df, life_year_change_multiplier) {
  if (identical(life_year_change_multiplier, 1)) {
    # English: For exposed -> reference readings, the "before" state is the
    # exposed arm and the "after" state is the reference arm.
    # 日本語: 曝露者が非曝露へ変わる読み方では、before=曝露, after=非曝露。
    return(df |>
      mutate(
        le_before = .data[["le_exposed"]],
        le_after = .data[["le_reference"]]
      ))
  }

  # English: For reference -> exposed readings and partial-change natural-vs-
  # changed comparisons, the original reference/natural arm is "before" and
  # the changed/exposed arm is "after".
  # 日本語: 非曝露から曝露への読み方やpartial_changeでは、before=自然/非曝露側,
  # after=介入後/曝露側として表示する。
  df |>
    mutate(
      le_before = .data[["le_reference"]],
      le_after = .data[["le_exposed"]]
    )
}
