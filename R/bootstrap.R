# Resampling copies must have distinct IDs, even when the same person is drawn twice.
# Retain the established sampling order to make comparisons with older versions easy.
yll_resample_people <- function(data, id_var) {
  ids <- data[[id_var]]
  sampled_ids <- sample(ids, length(ids), replace = TRUE)
  sample_data <- data[match(sampled_ids, ids), , drop = FALSE]
  copy_number <- ave(seq_along(sampled_ids), sampled_ids, FUN = seq_along)
  sample_data[[id_var]] <- paste0(sampled_ids, "_rep", copy_number)
  sample_data
}

yll_run_bootstrap <- function(B, estimate_once, data, id_var,
                              use_future, show_progress) {
  if (B == 0L) return(list())
  one_sample <- function(iteration) {
    # Keep failures visible instead of silently narrowing the bootstrap distribution.
    tryCatch({
      sampled_data <- yll_resample_people(data, id_var)
      result <- estimate_once(sampled_data)
      result$estimates$iteration <- iteration
      result$curves$iteration <- iteration
      result
    }, error = function(error) {
      list(failure = tibble::tibble(iteration = iteration, reason = conditionMessage(error)))
    })
  }
  if (use_future) {
    if (show_progress) {
      return(progressr::with_progress({
        progress <- progressr::progressor(steps = B)
        future.apply::future_lapply(seq_len(B), function(iteration) {
          result <- one_sample(iteration)
          progress()
          result
        }, future.seed = TRUE)
      }))
    }
    return(future.apply::future_lapply(seq_len(B), one_sample, future.seed = TRUE))
  }
  if (show_progress) {
    progress <- utils::txtProgressBar(min = 0, max = B, style = 3)
    on.exit(close(progress), add = TRUE)
    return(lapply(seq_len(B), function(iteration) {
      result <- one_sample(iteration)
      utils::setTxtProgressBar(progress, iteration)
      result
    }))
  }
  lapply(seq_len(B), one_sample)
}

# Both interval methods use the same nonparametric participant bootstrap.
yll_confidence_intervals <- function(estimates, bootstrap_estimates, conf_level, ci_method) {
  estimates$yll_se <- NA_real_
  estimates$ci_low <- NA_real_
  estimates$ci_high <- NA_real_
  alpha <- (1 - conf_level) / 2
  for (row in seq_len(nrow(estimates))) {
    values <- bootstrap_estimates$yll[
      bootstrap_estimates$starting_age == estimates$starting_age[row]
    ]
    if (length(values) < 2L) next
    estimates$yll_se[row] <- stats::sd(values)
    if (ci_method == "normal") {
      margin <- stats::qnorm(1 - alpha) * estimates$yll_se[row]
      estimates$ci_low[row] <- estimates$yll[row] - margin
      estimates$ci_high[row] <- estimates$yll[row] + margin
    } else {
      bounds <- stats::quantile(values, c(alpha, 1 - alpha), names = FALSE)
      estimates$ci_low[row] <- bounds[1]
      estimates$ci_high[row] <- bounds[2]
    }
  }
  estimates
}

yll_build_result <- function(point, replicates, meta) {
  failures <- bind_rows(lapply(replicates, function(result) result$failure))
  if (nrow(failures) == 0L) {
    failures <- tibble::tibble(iteration = integer(), reason = character())
  }
  succeeded <- vapply(replicates, function(result) is.null(result$failure), logical(1))
  successful <- replicates[succeeded]
  bootstrap_estimates <- point$estimates[0, ]
  bootstrap_estimates$iteration <- integer()
  bootstrap_curves <- point$curves[0, ]
  bootstrap_curves$iteration <- integer()
  if (length(successful) > 0L) {
    bootstrap_estimates <- bind_rows(lapply(successful, function(result) result$estimates))
    bootstrap_curves <- bind_rows(lapply(successful, function(result) result$curves))
  }
  meta$bootstrap_successful <- sum(succeeded)
  meta$bootstrap_failures <- failures
  # A model can return finite numbers even when its survival predictions are
  # invalid. Keep the raw results for diagnosis, but never present them silently.
  window <- meta$age_end - bootstrap_estimates$starting_age
  outside_window <- bootstrap_estimates$erl_reference < -1e-8 |
    bootstrap_estimates$erl_exposed < -1e-8 |
    bootstrap_estimates$erl_reference > window + 1e-8 |
    bootstrap_estimates$erl_exposed > window + 1e-8
  meta$bootstrap_unreliable <- bootstrap_estimates[outside_window, ]
  if (nrow(meta$bootstrap_unreliable) > 0L) {
    warning("Bootstrap ERL estimates outside the possible age window were returned. ",
            "The model is numerically unstable; the resulting confidence intervals are unreliable. ",
            "Raw estimates are retained in meta$bootstrap_unreliable for diagnosis.", call. = FALSE)
  }
  if (nrow(failures) > 0L) {
    warning(nrow(failures), " of ", meta$B, " bootstrap replicates failed. ",
            "Intervals use the successful replicates and may be unreliable. ",
            "See result$meta$bootstrap_failures for iteration numbers and reasons. ",
            "First failure: ", failures$reason[1], call. = FALSE)
  }
  if (meta$B > 0L && sum(succeeded) < 2L) {
    warning("Fewer than two bootstrap replicates succeeded; confidence intervals are NA.", call. = FALSE)
  }
  list(
    summary = yll_confidence_intervals(point$estimates, bootstrap_estimates,
                                       meta$conf_level, meta$ci_method),
    bootstrap_estimates = bootstrap_estimates,
    survival_curves = point$curves,
    bootstrap_survival_curves = bootstrap_curves,
    meta = meta
  )
}
