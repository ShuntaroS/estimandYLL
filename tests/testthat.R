# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.

if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(estimandYLL)

  test_check("estimandYLL")
} else {
  message("Package 'testthat' is not installed; skipping tests.")
}
