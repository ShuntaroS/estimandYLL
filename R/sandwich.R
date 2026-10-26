# Empirical sandwich variance for the pooled logistic estimator.
#
# The pooled logistic score, summed within person, and the standardized
# survival probabilities form a block-triangular system of estimating
# equations (Zivich et al., arXiv:2504.13291). The variance therefore follows
# from person-level influence contributions, and the model is fitted once.
#
# For person i, the influence contribution to a standardized survival
# probability S(a | a0) is
#
#   I(i in target) {S_i(a | a0) - S(a | a0)} / n_target  +  q_i G,
#
# where q_i = U_i' V is the person's score times the model covariance and G is
# the derivative of S(a | a0) with respect to the coefficients. ERL and YLL are
# linear in these probabilities, so one covariance matrix serves every output.
# Its pieces are accumulated over chunks of the target population, which keeps
# memory proportional to `chunk_rows`. Knot positions are treated as fixed.
#' @noRd
yll_sandwich_logistic <- function(fit, person_ids, data, id_var, exposure_var,
                                  reference_level, exposed_level,
                                  confounders_baseline, target_population,
                                  age_start, age_end, age_interval, integration,
                                  chunk_rows = 2e5) {
  keep <- !is.na(stats::coef(fit))
  beta <- stats::coef(fit)[keep]
  bread <- stats::vcov(fit)[keep, keep, drop = FALSE]
  if (any(!is.finite(bread))) {
    stop("The hazard model covariance is not finite, so the sandwich variance ",
         "cannot be computed. Inspect the model or use variance = \"bootstrap\".", call. = FALSE)
  }

  # Person-level score: person-period contributions summed within person.
  person <- match(person_ids, data[[id_var]])
  model_design <- stats::model.matrix(fit)[, keep, drop = FALSE]
  score <- rowsum(model_design * (fit$y - fit$fitted.values), person, reorder = TRUE)
  if (nrow(score) != nrow(data)) {
    stop("Every person must contribute at least one person-period record.", call. = FALSE)
  }
  influence_beta <- score %*% bread
  omega <- crossprod(influence_beta)

  target <- switch(target_population,
    all = seq_len(nrow(data)),
    exposed = which(data[[exposure_var]] == exposed_level),
    unexposed = which(data[[exposure_var]] == reference_level)
  )
  n_target <- length(target)
  ages <- seq.int(age_start, age_end)
  n_ages <- length(ages)
  starting_ages <- seq(age_start, age_end, by = age_interval)
  # One point per (starting age, age) pair; `first` and `column` index ages.
  points <- bind_rows(lapply(starting_ages, function(starting_age) {
    tibble::tibble(starting_age = starting_age, age = seq.int(starting_age, age_end))
  }))
  points$first <- as.integer(round(points$starting_age - age_start)) + 1L
  points$column <- as.integer(round(points$age - age_start)) + 1L
  n_points <- nrow(points)
  # Survival at the starting age is one for everyone, so those points have no
  # variance. Leaving them out of the gradient keeps it exactly zero.
  moving <- which(points$column > points$first)
  by_column <- split(moving, points$column[moving])
  by_first <- split(moving, points$first[moving])

  terms_rhs <- stats::delete.response(stats::terms(fit))
  design <- function(grid, level) {
    grid$expo <- level
    frame <- stats::model.frame(terms_rhs, grid, xlev = fit$xlevels)
    stats::model.matrix(terms_rhs, frame, contrasts.arg = fit$contrasts)[, keep, drop = FALSE]
  }

  n_coef <- length(beta)
  n_columns <- 2L * n_points
  # Cross-products are accumulated around the first chunk's means, which
  # avoids cancellation when survival probabilities are close to one.
  shift <- NULL
  sum_points <- numeric(n_columns)
  cross_points <- matrix(0, n_columns, n_columns)
  cross_influence <- matrix(0, n_columns, n_coef)
  gradient <- matrix(0, n_coef, n_columns)
  cumulative <- 1 * outer(seq_len(n_ages - 1L), seq_len(n_ages - 1L), "<=")
  people_per_chunk <- max(1L, floor(chunk_rows / n_ages))
  chunks <- split(target, ceiling(seq_along(target) / people_per_chunk))

  for (rows in chunks) {
    n_rows <- length(rows)
    grid <- data[rep(rows, each = n_ages), confounders_baseline, drop = FALSE]
    grid$age_temp <- rep(ages, times = n_rows)
    chunk_points <- matrix(0, n_rows, n_columns)
    levels <- list(reference_level, exposed_level)
    for (scenario in 1:2) {
      offset <- (scenario - 1L) * n_points
      z <- design(grid, levels[[scenario]])
      hazard <- stats::plogis(drop(z %*% beta))
      # Log survival to the start of each age, counted from age_start.
      log_step <- matrix(log1p(-hazard), n_rows, n_ages, byrow = TRUE)
      log_surv <- matrix(0, n_rows, n_ages)
      if (n_ages > 1L) log_surv[, -1L] <- log_step[, -n_ages, drop = FALSE] %*% cumulative
      surv <- exp(log_surv[, points$column, drop = FALSE] - log_surv[, points$first, drop = FALSE])
      chunk_points[, offset + seq_len(n_points)] <- surv

      # d S(a | a0) / d beta = -S(a | a0) {C(a) - C(a0)}, where C(a) sums
      # h_u z_u over ages u before a.
      weighted <- z * hazard
      age_rows <- (seq_len(n_rows) - 1L) * n_ages
      accumulated <- matrix(0, n_rows, n_coef)
      for (column in seq_len(n_ages)) {
        if (column > 1L) accumulated <- accumulated + weighted[age_rows + column - 1L, , drop = FALSE]
        to_age <- by_column[[as.character(column)]]
        if (!is.null(to_age)) {
          gradient[, offset + to_age] <- gradient[, offset + to_age] -
            crossprod(accumulated, surv[, to_age, drop = FALSE])
        }
        from_age <- by_first[[as.character(column)]]
        if (!is.null(from_age)) {
          gradient[, offset + from_age] <- gradient[, offset + from_age] +
            crossprod(accumulated, surv[, from_age, drop = FALSE])
        }
      }
    }
    if (is.null(shift)) shift <- colMeans(chunk_points)
    chunk_points <- chunk_points - rep(shift, each = n_rows)
    sum_points <- sum_points + colSums(chunk_points)
    cross_points <- cross_points + crossprod(chunk_points)
    cross_influence <- cross_influence + crossprod(chunk_points, influence_beta[rows, , drop = FALSE])
  }

  offset_mean <- sum_points / n_target
  estimate <- shift + offset_mean
  gradient <- gradient / n_target
  centered <- (cross_points - n_target * tcrossprod(offset_mean)) / n_target^2
  mixed <- (cross_influence - tcrossprod(offset_mean, colSums(influence_beta[target, , drop = FALSE]))) %*%
    gradient / n_target
  covariance <- centered + mixed + t(mixed) + crossprod(gradient, omega %*% gradient)

  # ERL is a weighted sum of the survival points, as in yll_integrate_le().
  weights <- yll_integration_weights(points, integration)
  erl_map <- matrix(0, length(starting_ages), n_points)
  erl_map[cbind(match(points$starting_age, starting_ages), seq_len(n_points))] <- weights
  empty <- matrix(0, length(starting_ages), n_points)
  maps <- list(reference = cbind(erl_map, empty), exposed = cbind(empty, erl_map))
  maps$yll <- maps$reference - maps$exposed
  standard_error <- function(map) sqrt(pmax(0, rowSums((map %*% covariance) * map)))
  point_se <- sqrt(pmax(0, diag(covariance)))

  list(
    estimates = tibble::tibble(
      starting_age = starting_ages,
      erl_reference = drop(maps$reference %*% estimate),
      erl_exposed = drop(maps$exposed %*% estimate),
      yll = drop(maps$yll %*% estimate)
    ),
    standard_errors = tibble::tibble(
      starting_age = starting_ages,
      erl_reference_se = standard_error(maps$reference),
      erl_exposed_se = standard_error(maps$exposed),
      yll_se = standard_error(maps$yll)
    ),
    curves = tibble::tibble(
      starting_age = points$starting_age,
      age = points$age,
      survival_reference_se = point_se[seq_len(n_points)],
      survival_exposed_se = point_se[n_points + seq_len(n_points)]
    )
  )
}

# Integration weights on the one-year grid, matching yll_integrate_le().
#' @noRd
yll_integration_weights <- function(points, integration) {
  weights <- rep(1, nrow(points))
  for (rows in split(seq_len(nrow(points)), points$starting_age)) {
    last <- rows[length(rows)]
    if (length(rows) == 1L) {
      weights[rows] <- 0
    } else if (integration == "left_rectangle") {
      weights[last] <- 0
    } else {
      weights[c(rows[1L], last)] <- 0.5
    }
  }
  weights
}

# Normal-approximation intervals from sandwich standard errors.
#' @noRd
yll_sandwich_intervals <- function(estimates, sandwich, conf_level) {
  standard_errors <- sandwich$standard_errors
  estimates$yll_se <- standard_errors$yll_se[match(estimates$starting_age, standard_errors$starting_age)]
  margin <- stats::qnorm(1 - (1 - conf_level) / 2) * estimates$yll_se
  estimates$ci_low <- estimates$yll - margin
  estimates$ci_high <- estimates$yll + margin
  estimates
}
