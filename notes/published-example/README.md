# Published example with confidence intervals

The README, English and Japanese introductory articles, and `plot_yll()` help
use a saved result from the full 5,000-person synthetic `yll_toy` dataset.
No observed participant data are distributed.

## Settings

* Method: pooled logistic g-formula.
* Target population: people observed without hypertension (`unexposed`).
* Exposure: hypertension; reference `No`, exposed `Yes`.
* Covariates: sex, education_years, bmi, smoke_binary.
* Starting ages: 40 through 90 years in five-year steps; upper age 90.
* Integration: annual left rectangles.
* Bootstrap: 1,000 participant resamples, seed 20260917.
* Execution: `future::multisession` with two workers and `use_future = TRUE`.
* Intervals: 95% pointwise normal intervals from the bootstrap standard error.

Run completed on 2026-09-17 with R 4.6.1 on macOS: all 1,000 replicates
succeeded, with no failed replicates or impossible-ERL flags. The separate
two-worker reproducibility check also passed. The published error-bar endpoints
were verified against the saved summary table.

The complete fitted result, including all bootstrap estimates and survival
curves, is saved in `inst/extdata/yll-example.rds`. Documentation builds read
that file rather than repeating the bootstrap. All displayed plots are drawn
by the package's plotting functions from that result.

The age-90 estimate and interval are zero because the integration window has
zero width. This is an illustration using synthetic data, not a reproduction
of the manuscript's cohort findings.

## Reproduce

Install the current package and its optional `future` and `ggplot2` dependencies,
then run the following from the repository root:

```sh
Rscript notes/published-example/generate-example.R
Rscript notes/published-example/publish-figures.R
Rscript notes/published-example/verify-future.R
```

The first script skips computation if `local/published-example/result.rds`
already exists. To deliberately repeat the fit, move that file elsewhere first.
The second script checks the settings, successful replicate count, failure
records, numerical reliability flags, and plotted CI endpoints before copying
the result and saving the figures.

The last script performs a separate four-replicate check with 1,000 people.
It compares two-worker execution, a repeated two-worker run, and the sequential
future backend. It requires identical replicate estimates, survival curves,
summaries and intervals, with all four replicates successful in each run.

The [summary table](summary.csv) contains each starting age's ERLs, YLL, standard
error and interval. Session details and run logs are retained locally under
`local/`, outside the source package.

## Relationship to the 100-replicate comparison

The previous [old-versus-new validation](../validation/README.md) remains a
separate 100-replicate comparison. This 1,000-replicate run supplies the public
figures. Its point estimates agree at shared starting ages; its confidence
limits need not be identical because the bootstrap count and execution backend
differ. The Poisson and Royston-Parmar comparisons are not rerun here.
