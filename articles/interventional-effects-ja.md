# estimandYLL における interventional effects

この vignette では、`estimandYLL` が interventional effects をどのように
表現するかを説明します。この package では、estimand
は単なる曝露コントラスト ではなく、**target population**
も必要であることを重視します。

``` text
Estimand = target_population + intervention + measure
```

## Target population

`target_population` は、誰の反事実的な平均余命を平均するかを指定します。

- `"all"`: 研究対象集団の全員。
- `"exposed"`: 観察上 exposed level にいた人。
- `"unexposed"`: 観察上 reference level にいた人。

ハザードモデルは全データで fit します。`target_population`
が変えるのは、 反事実的な survival curve
を作った後の平均化のステップだけです。

## Interventions

`intervention = "full_contrast"` は、2つの static worlds を比較します。

``` math
\mathrm{LE}^{a=\mathrm{reference}}
\quad \text{versus} \quad
\mathrm{LE}^{a=\mathrm{exposed}}.
```

`intervention = "partial_change"` は、自然な観察曝露分布と、 policy-like
な partial change を比較します。例えば次のような介入です。

- 観察上の曝露者: 30%が reference に移り、70%は exposed のまま。
- 観察上の非曝露者: 非曝露のまま。

内部的には、これは二値曝露に対する stochastic intervention です。
各対象者について、推定エンジンは次の survival curve を作ります。

``` math
S_i^*(t) =
(1 - p_i) S_i^{a=\mathrm{reference}}(t)
+ p_i S_i^{a=\mathrm{exposed}}(t),
```

そのうえで、指定された target population 内で $`S_i^*(t)`$
を平均します。

## Measure

`measure = "yll"` は次を報告します。

``` math
\mathrm{LE}_{reference} - \mathrm{LE}_{exposed}.
```

これは通常の years-of-life-lost
の向きです。有害曝露では正の値になりやすい 表示です。

`measure = "life_year_change"`
は、指定した介入方向に沿った臨床的に解釈しやすい
平均余命の変化を報告します。例えば、観察上の非曝露者が曝露状態になるという
ATC-like な問いでは、有害曝露の影響は負の値になります。

## Worked example

観察上の非高血圧者の20%が高血圧になり、観察上の高血圧者は高血圧のまま、
という programme を考えます。target population
は明示的に観察上の非高血圧者です。

``` r

library(estimandYLL)
data(yll_toy)

res <- estimand_yll(
  data = yll_toy,
  target_population = "unexposed",
  intervention = "partial_change",
  change_from = "reference",
  change_to = "exposed",
  change_probability = 0.20,
  measure = "life_year_change",
  B = 50, show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)

res$summary
```

## Identification assumptions

通常の g-formula と同じ仮定が必要です。

- 仮定した regime のもとでの counterfactual outcomes に対する
  consistency。
- 測定された共変量で条件づけた conditional exchangeability。
- 選択した intervention と target population に対する positivity。
- discrete-time hazard model の正しい指定。

Partial-change policy
は、全員を1つの曝露状態に固定する介入よりも現実的な場合が あります。特に
full intervention が非現実的であったり、positivity の問題を強める
場合には有用です。
