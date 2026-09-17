# Validate inputs before fitting, so errors refer to the user's data and arguments.
yll_check_number <- function(value, name, lower = -Inf, upper = Inf, integer = FALSE) {
  valid <- is.numeric(value) && length(value) == 1L && is.finite(value)
  if (valid) valid <- value >= lower && value <= upper
  if (valid && integer) valid <- value == floor(value)
  if (!valid) stop("`", name, "` must be a finite ",
                   if (integer) "integer" else "number", " in [", lower, ", ", upper, "].", call. = FALSE)
}

yll_check_inputs <- function(data, id_var, time_var, event_var, exposure_var,
                             age_at_entry_var, confounders_baseline,
                             age_start, age_end, age_interval, B, seed,
                             conf_level, show_progress, use_future,
                             model, prediction_interval, rp_df) {
  if (!is.data.frame(data) || nrow(data) < 2L || anyDuplicated(names(data))) {
    stop("`data` must be a data frame with at least two people and unique column names.", call. = FALSE)
  }
  column_arguments <- list(id_var = id_var, time_var = time_var, event_var = event_var,
                           exposure_var = exposure_var, age_at_entry_var = age_at_entry_var)
  for (name in names(column_arguments)) {
    value <- column_arguments[[name]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("`", name, "` must be one column name.", call. = FALSE)
    }
  }
  if (!is.null(confounders_baseline) &&
      (!is.character(confounders_baseline) || anyNA(confounders_baseline))) {
    stop("`confounders_baseline` must be a character vector of column names.", call. = FALSE)
  }
  columns <- unique(c(unlist(column_arguments), confounders_baseline))
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns)) stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  if (anyNA(data[, columns, drop = FALSE])) {
    stop("Selected variables contain missing values. Prepare a complete-case dataset explicitly.", call. = FALSE)
  }
  for (name in columns) {
    if (is.numeric(data[[name]]) && any(!is.finite(data[[name]]))) {
      stop("Column `", name, "` contains non-finite values.", call. = FALSE)
    }
  }
  if (anyDuplicated(data[[id_var]])) stop("`id_var` must identify one unique person per row.", call. = FALSE)
  if (!is.numeric(data[[time_var]]) || any(data[[time_var]] <= 0)) {
    stop("Follow-up duration must be positive and numeric, in years.", call. = FALSE)
  }
  if (!is.numeric(data[[age_at_entry_var]]) || any(data[[age_at_entry_var]] < 0)) {
    stop("Age at entry must be nonnegative and numeric, in years.", call. = FALSE)
  }
  if (!(is.numeric(data[[event_var]]) || is.logical(data[[event_var]])) ||
      any(!data[[event_var]] %in% c(0, 1))) {
    stop("The event indicator must be coded 0/1.", call. = FALSE)
  }
  if (!any(data[[event_var]] == 1)) stop("At least one observed event is required.", call. = FALSE)
  integer_ages <- model == "pooled_logistic"
  yll_check_number(age_start, "age_start", lower = 0, integer = integer_ages)
  yll_check_number(age_end, "age_end", lower = age_start, integer = integer_ages)
  yll_check_number(age_interval, "age_interval", lower = .Machine$double.eps, integer = integer_ages)
  yll_check_number(B, "B", lower = 0, integer = TRUE)
  if (B == 1) stop("`B` must be 0 or at least 2.", call. = FALSE)
  yll_check_number(seed, "seed", lower = 0, upper = .Machine$integer.max, integer = TRUE)
  yll_check_number(conf_level, "conf_level", lower = 0, upper = 1)
  if (conf_level == 0 || conf_level == 1) stop("`conf_level` must be strictly between 0 and 1.", call. = FALSE)
  yll_check_number(prediction_interval, "prediction_interval", lower = .Machine$double.eps)
  yll_check_number(rp_df, "rp_df", lower = 1, integer = TRUE)
  for (value in list(show_progress, use_future)) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop("Progress and parallel options must be TRUE or FALSE.", call. = FALSE)
    }
  }
  # Formula construction and intermediate columns need unambiguous variable names.
  reserved <- c("tstart", "tgroup", "age_temp", "expo", "expo_original",
                ".age_in", ".age_out", ".age_mid", ".expo", ".event",
                ".death", ".prediction_id", ".prediction_row", "lex.dur", "lex.Xst")
  if (any(columns %in% reserved)) stop("Selected column names conflict with internal model columns: ",
                                     paste(intersect(columns, reserved), collapse = ", "), call. = FALSE)
  if (any(make.names(columns) != columns)) {
    stop("Selected columns must have syntactically valid R names; rename them before estimation.", call. = FALSE)
  }
  if (age_at_entry_var %in% confounders_baseline) {
    warning("Age at entry is included among baseline covariates, in addition to the attained-age spline. ",
            "Check that this is intentional.", call. = FALSE)
  }
  data[, columns, drop = FALSE]
}

yll_exposure_levels <- function(exposure, reference_level, exposed_level) {
  observed <- unique(exposure)
  if (length(observed) != 2L) stop("Exactly two observed exposure levels are required.", call. = FALSE)
  ordered <- if (is.factor(exposure)) levels(droplevels(exposure)) else observed
  inferred <- is.null(reference_level) || is.null(exposed_level)
  if (is.null(reference_level) && is.null(exposed_level)) {
    reference_level <- ordered[1]
    exposed_level <- ordered[2]
  } else if (is.null(reference_level)) {
    reference_level <- ordered[!ordered %in% exposed_level]
  } else if (is.null(exposed_level)) {
    exposed_level <- ordered[!ordered %in% reference_level]
  }
  if (length(reference_level) != 1L || length(exposed_level) != 1L ||
      anyNA(c(reference_level, exposed_level)) ||
      !reference_level %in% observed || !exposed_level %in% observed ||
      as.character(reference_level) == as.character(exposed_level)) {
    stop("Reference and exposed levels must be distinct observed exposure values.", call. = FALSE)
  }
  if (inferred) message("Exposure levels: reference = '", reference_level,
                        "', exposed = '", exposed_level, "'. Specify both explicitly to set the direction.")
  list(reference = reference_level, exposed = exposed_level)
}
