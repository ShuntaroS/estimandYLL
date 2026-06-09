# estimandYLL 入門

`estimandYLL` は、パラメトリック g-formula に基づいて **years of life
lost (YLL)** と反事実的な平均余命を推定する package です。 主な入口は
[`estimand_yll()`](https://shuntaros.github.io/estimandYLL/reference/estimand_yll.md)
です。

この package では、estimand を以下の3つの引数で明示します。

1.  `target_population`: 誰で平均するか。
2.  `intervention`: どの曝露世界、またはどの policy change
    を比較するか。
3.  `measure`: 結果をどの向きで表示するか。

``` r

library(estimandYLL)
data(yll_toy)
str(yll_toy)
#> 'data.frame':    5000 obs. of  9 variables:
#>  $ id             : int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ period         : num  15 15 15 15 15 ...
#>  $ event          : int  0 0 1 0 0 0 0 0 0 1 ...
#>  $ smoke_binary   : Factor w/ 2 levels "Never","Current/Ever": 1 1 1 2 2 2 2 1 1 2 ...
#>  $ hypertension   : Factor w/ 2 levels "No","Yes": 1 2 1 1 1 1 1 1 2 2 ...
#>  $ sex            : Factor w/ 2 levels "Female","Male": 2 2 2 1 1 2 1 2 1 1 ...
#>  $ education_years: num  14 15 16 11 12 15 18 14 19 14 ...
#>  $ bmi            : num  22 23 22 24.7 26.6 27.6 26.3 25.6 22.7 18.7 ...
#>  $ age            : num  54 49 64 58 47 61 51 57 43 65 ...
```

## Full contrasts

`intervention = "full_contrast"` は、全員が reference level
である反事実世界と、 全員が exposed level
である反事実世界を比較します。`target_population` によって、 結果が
ATE-like、ATT-like、ATC-like のどれに対応するかが決まります。

``` r

res_all <- estimand_yll(
  data = yll_toy,
  target_population = "all",
  intervention = "full_contrast",
  measure = "yll",
  B = 50,
  method = "normal",
  show_progress = FALSE,
  id_var = "id",
  time_var = "period",
  event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No",
  exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50,
  age_end = 90,
  age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary"),
  use_future = FALSE
)

res_all$summary
```

観察上の曝露者、観察上の非曝露者を target population
として明示できます。

``` r

res_exposed <- estimand_yll(
  data = yll_toy,
  target_population = "exposed",
  intervention = "full_contrast",
  measure = "yll",
  B = 50, method = "normal", show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)

res_unexposed <- estimand_yll(
  data = yll_toy,
  target_population = "unexposed",
  intervention = "full_contrast",
  measure = "life_year_change",
  B = 50, method = "normal", show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "hypertension",
  reference_level = "No", exposed_level = "Yes",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi", "smoke_binary")
)
```

`measure = "life_year_change"` を使うと、ATC-like な有害方向の変化は
負の値として表示されます。

## Partial-change policies

`intervention = "partial_change"` は、自然な観察曝露分布と、
一定割合の対象者だけ曝露状態が変わる policy-like な介入を比較します。
統計学的には stochastic intervention に相当しますが、API では臨床家にも
読みやすい policy language を使っています。

例: 観察上の曝露者の30%が reference level に移る policy を、
観察上の曝露者で平均します。

``` r

res_quit <- estimand_yll(
  data = yll_toy,
  target_population = "exposed",
  intervention = "partial_change",
  change_from = "exposed",
  change_to = "reference",
  change_probability = 0.30,
  measure = "life_year_change",
  B = 50, show_progress = FALSE, use_future = FALSE,
  id_var = "id", time_var = "period", event_var = "event",
  exposure_var = "smoke_binary",
  reference_level = "Never", exposed_level = "Current/Ever",
  age_at_entry_var = "age",
  age_start = 50, age_end = 90, age_interval = 5,
  confounders_baseline = c("sex", "education_years", "bmi")
)
```

## 可視化

``` r

plot_yll_estimate(res_all)
plot_conditional_survival(res_all, age_start = 60)
plot_marginal_survival(res_all)
```

## 結果オブジェクト

`summary` には開始年齢ごとの結果が入り、選択された
`estimate`、信頼区間、 `target_population`、`intervention`、`measure`
が含まれます。 `detailed_results` には、元の `yll`、life-year change、各
arm の平均余命も 確認用に保持されます。
