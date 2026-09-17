# Validation of estimandYLL 0.1.0

The refactor preserves the previous estimator's numerical results on the full
5,000-person `yll_toy` dataset. It also exposes a pre-existing numerical
instability in the Royston-Parmar bootstrap; that model's intervals in this
validation are not usable as uncertainty estimates.

## Baseline and settings

Baseline source: GitHub commit `aa9973ec0b76ee409aae462699f13350fb1a266a`.
Validation date: 2026-09-17. Runtime: R 4.6.1 on macOS (Apple silicon).
Old and new estimators were run in separate R sessions. The baseline result
builders were instrumented only to retain bootstrap estimates that the old API
did not expose. The estimation algorithms were not modified for the comparison.

* Data: all 5,000 people in the unchanged packaged synthetic dataset.
* Exposure: hypertension, reference `No`, exposed `Yes`.
* Covariates: sex, education_years, bmi, smoke_binary.
* Starting ages: 40, 50, 60, 70, 80, 90 years; upper age: 90 years.
* Integration: left rectangles on a one-year prediction grid.
* Bootstrap: 100 participant resamples, seed 20260917, sequential execution.
* Confidence level: 95%. Normal and percentile intervals for the main method;
  normal intervals for both comparison models.
* Royston-Parmar model: combined model, four baseline spline degrees of freedom.

Percentile intervals were computed from the stored resamples rather than
running another bootstrap. Unit tests separately exercise the public
`ci_method = "percentile"` path.

## Results

| Setting | Successful replicates, old / new | Maximum point-estimate difference | Maximum CI endpoint difference |
|---|---:|---:|---:|
| Pooled logistic, all | 100 / 100 | 0 years | 0 years |
| Pooled logistic, exposed | 100 / 100 | 0 years | 0 years |
| Pooled logistic, unexposed | 100 / 100 | 0 years | 8.89e-16 years |
| Poisson | 100 / 100 | 0 years | 8.89e-16 years |
| Royston-Parmar | 93 / 93 | 0 years | 0.015625 years |

The main-method percentile interval endpoints match exactly. Point survival
curves match exactly for every setting. All point-estimate ERLs and YLLs match.
Strict equality was not required; the near-zero differences satisfy the intended
level of agreement. Timing columns describe these runs only: competing local
processes make them unsuitable as a performance benchmark.

The Royston-Parmar bootstrap YLL values themselves match exactly, including
unstable values as large as 9.63e15 years in absolute magnitude. CI endpoints
reach approximately +/-2.55e15 years at starting age 40. Against that scale,
the 0.015625-year endpoint difference reflects floating-point aggregation,
not a meaningful improvement or deterioration. The intervals are numerically
unstable in **both** versions and should not be interpreted.

Replicates 6, 7, 22, 55, 62, 64, and 66 failed in both versions. The current
version records the error `function cannot be evaluated at initial parameters`.
The old estimator discarded the reasons. Its failed iteration numbers were
recovered from missing replicate IDs. No additional bootstrap replicates were
silently substituted for failures.

Warning counts in the CSV files describe the numerical comparison runs, before
the additional impossible-ERL diagnostics were added. Those diagnostics do not
change estimates, intervals, or the successful replicate set.

The refactored package reports failed replicates and flags finite ERLs outside
the possible age window. It retains the raw unstable estimates for diagnosis
rather than changing the model or silently removing extreme values.

## Files

* [Estimate and interval comparisons](estimate-comparison.csv)
* [Maximum differences by setting and metric](difference-summary.csv)
* [Point survival-curve comparisons](survival-comparison.csv)
* [Replicate success counts and elapsed times](run-status.csv)
* [Failed replicate numbers and reasons](bootstrap-failures.csv)
* [Run old and new estimators](run-comparison.R)
* [Summarize stored runs](summarize-comparison.R)

The full RDS objects and console logs are local build artifacts under
`local/validation/`, excluded from Git and the source package. The CSV summaries
and scripts above are version controlled. They contain synthetic results only.

## Reproduce

From the package repository root, prepare the fixed baseline source:

```sh
mkdir -p local/legacy
git archive aa9973ec0b76ee409aae462699f13350fb1a266a | tar -x -C local/legacy
Rscript notes/validation/run-comparison.R legacy
Rscript notes/validation/run-comparison.R current
Rscript notes/validation/summarize-comparison.R
```

The runner skips completed RDS outputs so an interrupted validation can resume.
Use a new output directory or remove only those validation RDS files before
requesting a fresh full comparison. The legacy and current sessions must remain
separate; do not load both package namespaces in one R session.

## Additional checks

The automated tests cover target-population weighting, sign reversal when the
exposure labels are swapped, zero-width windows, conditioning before averaging,
invariance to the earliest requested age, interval construction, unique IDs for
resampled copies, missing data and invalid arguments, explicit failure records,
and editable ggplot objects with and without intervals.

An installed-package check with two future workers confirms reproducibility
across repeated runs and between the two-worker and sequential future plans.
The ordinary sequential bootstrap can use a different random stream.
Large known-model simulation tests run in development; small functional tests
remain enabled on CRAN.
