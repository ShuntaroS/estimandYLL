## Release

First planned CRAN submission, version 0.1.0. The package has not been submitted.

## Local test environment

* macOS, Apple silicon, R 4.6.1
* R CMD check --as-cran, including the PDF manual
* Full development test suite and installed-package future tests with two workers

Post-publication check on 2026-09-17: 0 ERRORs, 0 WARNINGs, 2 NOTEs.
The notes concern a new submission and an old system HTML Tidy executable.
The earlier website URL note was resolved by publishing the pkgdown site.
The local HTML Tidy note concerns the available validator, not an R code or
documentation build failure. Examples, tests, vignettes and the PDF manual pass.

## Statistical validation

The full 5,000-person synthetic dataset was evaluated against commit aa9973e
in separate R sessions, with 100 participant bootstrap replicates per setting.
Point estimates and point survival curves match for all five settings. Main
method and Poisson intervals agree to floating-point precision. Normal and
percentile intervals were checked for the main method.

The Royston-Parmar comparator is numerically unstable on this validation dataset
in both versions: the same seven of 100 replicates fail, and some other replicates
produce impossible lifetimes. Its raw results are preserved and explicitly
flagged rather than silently discarded. These intervals must not be interpreted.
Full comparison tables and reproduction scripts are in notes/validation/ on
GitHub, excluded from the CRAN source tarball.

## Runtime and optional dependencies

Examples use small synthetic subsets and no bootstrap. Long known-model tests
are skipped on CRAN; small tests exercise all three estimators when optional
rstpm2 is available. ggplot2 and rstpm2 are optional dependencies. No network
access is required to estimate, test, or build the vignettes. Parallel workers
are controlled by the caller; examples do not create a worker pool.
