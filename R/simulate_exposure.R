#' Simulate exposure profiles
#'
#' Generates one or more exposure profiles in chronological order, from the
#' earliest observation to the most recent observation. The first profile value
#' corresponds internally to the maximum retrospective lag, whereas the final
#' profile value corresponds internally to lag 0.
#'
#' All user-facing pattern arguments are expressed in chronological time rather
#' than retrospective lag. Therefore, `time = 0` always identifies the earliest
#' observation and `time = lag_max` identifies the most recent observation.
#'
#' @param lag_max Non-negative integer maximum lag. The generated profile has
#'   length `lag_max + 1` and chronological names from `time_0` to
#'   `time_<lag_max>`.
#' @param n Positive integer number of profiles to simulate. Default is 1.
#' @param mode Character. `"profile"` simulates only the background process;
#'   `"pattern"` also applies the selected chronological pattern.
#' @param pattern Character pattern used when `mode = "pattern"`:
#'   `"random_times"`, `"alternating_times"`, or `"block_times"`.
#' @param fixed_value Finite numeric scalar assigned to selected chronological
#'   times for `"random_times"` and `"block_times"`.
#' @param n_times Positive integer number of chronological times selected by
#'   `"random_times"`.
#' @param time_range Optional integer vector identifying eligible chronological
#'   times. Values must lie between 0 and `lag_max`. If `NULL`, every
#'   chronological time is eligible.
#' @param alternating_values Finite numeric vector with at least two values,
#'   recycled chronologically over `time_range` for `"alternating_times"`.
#' @param block_times Integer vector of length 2 defining the first and final
#'   chronological times of the inclusive block used by `"block_times"`.
#' @param background Named list controlling the background process. Supported
#'   distributions are `"normal"`, `"empirical"`, `"ar1"`, and `"fixed"`.
#' @param bounds Optional finite numeric vector of length 2. Bounds are applied
#'   to the simulated background and again after applying the pattern.
#' @param seed Optional finite integer seed. When supplied, the complete set of
#'   simulations is reproducible.
#' @param cumulative Logical. If `TRUE`, applies `cumsum()` in chronological
#'   order after the pattern and bounds have been applied. In this case,
#'   `bounds` constrain the exposure increments, not the final cumulative values.
#'
#' @return If `n = 1`, a list with `profile` and `meta`. If `n > 1`, a list with
#'   `profiles` and `meta`. Every profile is ordered chronologically and named
#'   `time_0`, `time_1`, ..., `time_<lag_max>`.
#'
#' @details
#' The returned order is always:
#'
#' `time_0` = earliest observation = internal lag `lag_max`
#'
#' `time_<lag_max>` = most recent observation = internal lag 0
#'
#' The backward lag mapping is stored in `meta$retrospective_lag` and is not
#' exposed through the pattern arguments. Users specify all patterns in
#' chronological time.
#'
#' @export
simulate_exposure <- function(
    lag_max,
    n = 1,
    mode = c("profile", "pattern"),
    pattern = c("random_times", "alternating_times", "block_times"),
    fixed_value = NULL,
    n_times = NULL,
    time_range = NULL,
    alternating_values = NULL,
    block_times = NULL,
    background = list(dist = "normal", mean = 0, sd = 1),
    bounds = NULL,
    seed = NULL,
    cumulative = FALSE
) {

  # =========================================================
  # HELPERS
  # =========================================================

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  is_whole_number <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      abs(x - round(x)) < sqrt(.Machine$double.eps)
  }

  validate_finite_scalar <- function(x, argument) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
      stop("`", argument, "` must be one finite numeric value.")
    }
    as.numeric(x)
  }

  validate_time_vector <- function(x, argument, lag_max_i) {
    if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)) ||
        any(abs(x - round(x)) >= sqrt(.Machine$double.eps))) {
      stop("`", argument, "` must contain finite integer times.")
    }

    x <- as.integer(round(x))

    if (any(x < 0L) || any(x > lag_max_i)) {
      stop(
        "`", argument, "` must contain chronological times between 0 and ",
        "`lag_max`."
      )
    }

    if (anyDuplicated(x)) {
      stop("`", argument, "` must contain unique chronological times.")
    }

    x
  }

  # =========================================================
  # BASIC VALIDATION
  # =========================================================

  if (!is_whole_number(lag_max) || lag_max < 0) {
    stop("`lag_max` must be a non-negative integer.")
  }
  lag_max_i <- as.integer(round(lag_max))

  if (!is_whole_number(n) || n <= 0) {
    stop("`n` must be a positive integer.")
  }
  n <- as.integer(round(n))

  if (!is.null(seed)) {
    if (!is_whole_number(seed)) {
      stop("`seed` must be NULL or one finite integer.")
    }
    set.seed(as.integer(round(seed)))
  }

  if (!is.logical(cumulative) || length(cumulative) != 1L ||
      is.na(cumulative)) {
    stop("`cumulative` must be TRUE or FALSE.")
  }

  mode_i <- match.arg(mode)

  if (mode_i == "pattern") {
    pattern_i <- match.arg(pattern)
  } else {
    pattern_i <- NA_character_
    if (!missing(pattern)) {
      message("`pattern` is ignored when `mode = 'profile'`.")
    }
  }

  if (!is.list(background)) {
    stop("`background` must be a list.")
  }

  if (!is.null(bounds)) {
    if (!is.numeric(bounds) || length(bounds) != 2L || anyNA(bounds) ||
        any(!is.finite(bounds))) {
      stop("`bounds` must be a finite numeric vector of length 2.")
    }

    bounds <- as.numeric(bounds)

    if (bounds[1L] > bounds[2L]) {
      stop("`bounds[1]` must be less than or equal to `bounds[2]`.")
    }
  }

  # =========================================================
  # CHRONOLOGICAL INDEX
  # =========================================================

  profile_length <- lag_max_i + 1L
  chronological_time <- seq.int(0L, lag_max_i)
  retrospective_lag <- rev(chronological_time)
  profile_names <- paste0("time_", chronological_time)

  if (is.null(time_range)) {
    time_range_i <- chronological_time
  } else {
    time_range_i <- validate_time_vector(
      time_range,
      "time_range",
      lag_max_i
    )
    time_range_i <- sort(time_range_i)
  }

  # =========================================================
  # BACKGROUND VALIDATION AND GENERATION
  # =========================================================

  dist <- background$dist %||% "normal"

  if (!is.character(dist) || length(dist) != 1L || is.na(dist) ||
      !nzchar(dist)) {
    stop("`background$dist` must be one non-empty character value.")
  }

  dist <- match.arg(
    tolower(dist),
    choices = c("normal", "empirical", "ar1", "fixed")
  )

  generate_background <- function() {
    if (dist == "normal") {
      mean_value <- validate_finite_scalar(
        background$mean %||% 0,
        "background$mean"
      )
      sd_value <- validate_finite_scalar(
        background$sd %||% 1,
        "background$sd"
      )

      if (sd_value < 0) {
        stop("`background$sd` must be non-negative.")
      }

      return(stats::rnorm(profile_length, mean_value, sd_value))
    }

    if (dist == "empirical") {
      values <- background$values

      if (!is.numeric(values) || !length(values)) {
        stop(
          "For `background$dist = 'empirical'`, `background$values` ",
          "must be a non-empty numeric vector."
        )
      }

      values <- values[is.finite(values)]

      if (!length(values)) {
        stop("`background$values` must contain at least one finite value.")
      }

      return(sample(values, size = profile_length, replace = TRUE))
    }

    if (dist == "ar1") {
      mean_value <- validate_finite_scalar(
        background$mean %||% 0,
        "background$mean"
      )
      sd_value <- validate_finite_scalar(
        background$sd %||% 1,
        "background$sd"
      )
      phi <- validate_finite_scalar(
        background$phi %||% 0.7,
        "background$phi"
      )

      if (sd_value < 0) {
        stop("`background$sd` must be non-negative.")
      }

      if (abs(phi) >= 1) {
        stop("For `background$dist = 'ar1'`, `background$phi` must lie in (-1, 1).")
      }

      x <- numeric(profile_length)
      x[1L] <- stats::rnorm(1L, mean_value, sd_value)

      innovation_sd <- if (sd_value == 0) {
        0
      } else {
        sd_value * sqrt(1 - phi^2)
      }

      if (profile_length > 1L) {
        for (i in 2:profile_length) {
          x[i] <- mean_value +
            phi * (x[i - 1L] - mean_value) +
            stats::rnorm(1L, 0, innovation_sd)
        }
      }

      return(x)
    }

    fixed_background <- validate_finite_scalar(
      background$value,
      "background$value"
    )

    rep(fixed_background, profile_length)
  }

  # =========================================================
  # PATTERN VALIDATION
  # =========================================================

  if (mode_i == "pattern" && pattern_i %in% c("random_times", "block_times")) {
    fixed_value <- validate_finite_scalar(fixed_value, "fixed_value")

    if (!is.null(bounds) &&
        (fixed_value < bounds[1L] || fixed_value > bounds[2L])) {
      warning(
        "`fixed_value` lies outside `bounds` and will be truncated.",
        call. = FALSE
      )
    }
  }

  if (mode_i == "pattern" && pattern_i == "random_times") {
    if (!is_whole_number(n_times) || n_times <= 0) {
      stop("`n_times` must be a positive integer for `pattern = 'random_times'`.")
    }

    n_times <- as.integer(round(n_times))

    if (n_times > length(time_range_i)) {
      stop(
        "`n_times` cannot exceed the number of chronological times available ",
        "in `time_range`."
      )
    }
  }

  if (mode_i == "pattern" && pattern_i == "alternating_times") {
    if (!is.numeric(alternating_values) || length(alternating_values) < 2L ||
        anyNA(alternating_values) || any(!is.finite(alternating_values))) {
      stop(
        "`alternating_values` must contain at least two finite numeric values ",
        "for `pattern = 'alternating_times'`."
      )
    }

    alternating_values <- as.numeric(alternating_values)

    if (!is.null(bounds) &&
        any(alternating_values < bounds[1L] |
            alternating_values > bounds[2L])) {
      warning(
        "Some `alternating_values` lie outside `bounds` and will be truncated.",
        call. = FALSE
      )
    }
  }

  if (mode_i == "pattern" && pattern_i == "block_times") {
    if (is.null(block_times) || length(block_times) != 2L) {
      stop(
        "`block_times` must contain exactly two chronological times for ",
        "`pattern = 'block_times'`."
      )
    }

    block_times <- validate_time_vector(
      block_times,
      "block_times",
      lag_max_i
    )

    block_times <- sort(block_times)
  }

  # =========================================================
  # INTERNAL PROFILE GENERATOR
  # =========================================================

  simulate_one <- function() {
    x <- generate_background()

    if (!is.null(bounds)) {
      x <- pmin(pmax(x, bounds[1L]), bounds[2L])
    }

    selected_times <- integer(0)

    if (mode_i == "pattern") {
      if (pattern_i == "random_times") {
        selected_times <- sort(
          sample(
            time_range_i,
            size = n_times,
            replace = FALSE
          )
        )

        x[selected_times + 1L] <- fixed_value
      }

      if (pattern_i == "alternating_times") {
        selected_times <- sort(time_range_i)

        chronological_values <- rep(
          alternating_values,
          length.out = length(selected_times)
        )

        x[selected_times + 1L] <- chronological_values
      }

      if (pattern_i == "block_times") {
        selected_times <- seq.int(
          block_times[1L],
          block_times[2L]
        )

        x[selected_times + 1L] <- fixed_value
      }
    }

    if (!is.null(bounds)) {
      x <- pmin(pmax(x, bounds[1L]), bounds[2L])
    }

    if (cumulative) {
      x <- cumsum(x)
    }

    names(x) <- profile_names

    list(
      profile = x,
      meta = list(
        lag_max = lag_max_i,
        profile_length = profile_length,
        profile_order = "chronological",
        chronological_time = chronological_time,
        retrospective_lag = retrospective_lag,
        mode = mode_i,
        pattern = pattern_i,
        background = background,
        bounds = bounds,
        seed = seed,
        time_range = time_range_i,
        selected_times = selected_times,
        selected_positions = selected_times + 1L,
        selected_internal_lags = lag_max_i - selected_times,
        cumulative = cumulative,
        bounds_apply_to = if (cumulative) {
          "increments_before_cumulative_sum"
        } else {
          "profile_values"
        }
      )
    )
  }

  # =========================================================
  # RUN SIMULATIONS
  # =========================================================

  simulations <- replicate(
    n,
    simulate_one(),
    simplify = FALSE
  )

  if (n == 1L) {
    return(simulations[[1L]])
  }

  list(
    profiles = lapply(simulations, `[[`, "profile"),
    meta = lapply(simulations, `[[`, "meta")
  )
}
