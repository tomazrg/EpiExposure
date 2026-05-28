#' Simulate an exposure-lag profile (full or patterned)
#'
#' Creates a numeric vector' @param n_lags Integer. Number of random lags (required for "random_lags").#' Creates a numeric vector of length (lag_max + 1) representing an exposure
#' @param lag_range Integer vector. Candidate lags (subset of 0:lag_max).
#'   If NULL, defaults to 0:lag_max.
#' @param alternating_values Numeric vector. Values to cycle through
#'   (required for "alternating_lags", length >= 2).
#' @param block_lags Integer length-2. c(lag_start, lag_end) (required for "block_lags").
#' @param background List controlling background generation. Supported:
#'   dist = "normal" | "empirical" | "ar1" | "fixed"
#'   - normal: mean, sd
#'   - empirical: values
#'   - ar1: mean, sd, phi
#'   - fixed: value
#' @param bounds Optional numeric length-2 c(min, max) to clamp values.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return list(profile = numeric, meta = list)
#' @export
simulate_exposure_profile <- function(
    lag_max,
    mode = c("profile", "pattern"),
    pattern = c("random_lags", "alternating_lags", "block_lags"),
    fixed_value = NULL,
    n_lags = NULL,
    lag_range = NULL,
    alternating_values = NULL,
    block_lags = NULL,
    background = list(dist = "normal", mean = 0, sd = 1),
    bounds = NULL,
    seed = NULL
) {
  # ---- helpers (local; avoids global operators) ----
  nullcoalesce <- function(a, b) if (!is.null(a)) a else b
  is_whole_number <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x) && abs(x - round(x)) < .Machine$double.eps^0.5

  # ---- validate core args ----
  if (!is_whole_number(lag_max) || lag_max < 0) {
    stop("lag_max must be a non-negative integer.")
  }
  lag_max <- as.integer(round(lag_max))

  mode <- match.arg(mode)
  pattern <- match.arg(pattern)

  if (!is.null(seed)) {
    if (!is_whole_number(seed)) stop("seed must be an integer if provided.")
    set.seed(as.integer(round(seed)))
  }

  N <- lag_max + 1L
  lags_all <- 0:lag_max

  # ---- validate bounds ----
  if (!is.null(bounds)) {
    if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) {
      stop("bounds must be numeric length-2 c(min, max).")
    }
    if (bounds[1] > bounds[2]) stop("bounds[1] must be <= bounds[2].")
  }

  # ---- validate background ----
  if (!is.list(background)) stop("background must be a list.")
  dist <- nullcoalesce(background$dist, "normal")
  if (!is.character(dist) || length(dist) != 1L) stop("background$dist must be a single string.")

  gen_background <- function() {
    if (dist == "normal") {
      mu <- nullcoalesce(background$mean, 0)
      sd <- nullcoalesce(background$sd, 1)
      if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) stop("background$mean must be a finite number.")
      if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd < 0) stop("background$sd must be a finite number >= 0.")
      stats::rnorm(N, mean = mu, sd = sd)

    } else if (dist == "empirical") {
      vals <- background$values
      if (is.null(vals) || !is.numeric(vals) || length(vals) == 0L) {
        stop("background$values (numeric) is required for dist='empirical'.")
      }
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0L) stop("background$values has no finite values.")
      sample(vals, size = N, replace = TRUE)

    } else if (dist == "ar1") {
      mu <- nullcoalesce(background$mean, 0)
      sd <- nullcoalesce(background$sd, 1)
      phi <- nullcoalesce(background$phi, 0.7)

      if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu)) stop("background$mean must be a finite number.")
      if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd < 0) stop("background$sd must be a finite number >= 0.")
      if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || abs(phi) >= 1) {
        stop("background$phi must be a finite number with |phi| < 1 for dist='ar1'.")
      }

      x <- numeric(N)
      x[1] <- stats::rnorm(1, mean = mu, sd = sd)

      # innovation sd to maintain marginal sd approximately = sd
      eps_sd <- if (sd == 0) 0 else sd * sqrt(1 - phi^2)

      for (i in 2:N) {
        x[i] <- mu + phi * (x[i - 1] - mu) + stats::rnorm(1, mean = 0, sd = eps_sd)
      }
      x

    } else if (dist == "fixed") {
      val <- background$value
      if (is.null(val) || !is.numeric(val) || length(val) != 1L || !is.finite(val)) {
        stop("background$value (finite numeric scalar) is required for dist='fixed'.")
      }
      rep(val, N)

    } else {
      stop("Unknown background dist: ", dist,
           ". Supported: 'normal', 'empirical', 'ar1', 'fixed'.")
    }
  }

  x <- gen_background()

  # clamp if requested (before pattern)
  if (!is.null(bounds)) {
    x <- pmin(pmax(x, bounds[1]), bounds[2])
  }

  fixed_lags <- integer(0)

  # ---- pattern overlay ----
  if (mode == "pattern") {

    if (is.null(lag_range)) {
      lag_range <- lags_all
    } else {
      if (!is.numeric(lag_range)) stop("lag_range must be numeric/integer.")
      lag_range <- as.integer(round(lag_range))
      lag_range <- lag_range[is.finite(lag_range)]
    }

    lag_range <- intersect(unique(lag_range), lags_all)
    if (length(lag_range) == 0L) stop("lag_range has no valid lags within 0:lag_max.")

    if (pattern == "random_lags") {
      if (is.null(n_lags) || !is_whole_number(n_lags) || n_lags <= 0) {
        stop("n_lags must be a positive integer for pattern='random_lags'.")
      }
      if (is.null(fixed_value) || !is.numeric(fixed_value) || length(fixed_value) != 1L || !is.finite(fixed_value)) {
        stop("fixed_value (finite numeric scalar) must be provided for pattern='random_lags'.")
      }

      n_lags <- as.integer(round(n_lags))
      n_lags <- min(n_lags, length(lag_range))
      fixed_lags <- sort(sample(lag_range, size = n_lags, replace = FALSE))
      x[fixed_lags + 1L] <- fixed_value

    } else if (pattern == "alternating_lags") {
      if (is.null(alternating_values) || !is.numeric(alternating_values) || length(alternating_values) < 2L) {
        stop("alternating_values (numeric, length >= 2) must be provided for pattern='alternating_lags'.")
      }
      alternating_values <- alternating_values[is.finite(alternating_values)]
      if (length(alternating_values) < 2L) stop("alternating_values must contain at least 2 finite values.")

      fixed_lags <- sort(lag_range)
      alt <- rep(alternating_values, length.out = length(fixed_lags))
      x[fixed_lags + 1L] <- alt

    } else if (pattern == "block_lags") {
      if (is.null(block_lags) || !is.numeric(block_lags) || length(block_lags) != 2L) {
        stop("block_lags must be numeric/integer length-2 c(lag_start, lag_end) for pattern='block_lags'.")
      }
      if (is.null(fixed_value) || !is.numeric(fixed_value) || length(fixed_value) != 1L || !is.finite(fixed_value)) {
        stop("fixed_value (finite numeric scalar) must be provided for pattern='block_lags'.")
      }

      a <- as.integer(round(min(block_lags)))
      b <- as.integer(round(max(block_lags)))
      a <- max(0L, a)
      b <- min(lag_max, b)

      fixed_lags <- a:b
      x[fixed_lags + 1L] <- fixed_value
    }
  }

  # clamp after overlay
  if (!is.null(bounds)) {
    x <- pmin(pmax(x, bounds[1]), bounds[2])
  }

  # name by lag for readability
  names(x) <- paste0("lag_", lags_all)

  list(
    profile = x,
    meta = list(
      lag_max = lag_max,
      mode = mode,
      pattern = if (mode == "pattern") pattern else NA_character_,
      background = background,
      bounds = bounds,
      seed = seed,
      fixed_lags = fixed_lags
    )
  )
}

#' Simulate exposure history across lags (alias)
#'
#' User-facing short alias to \code{simulate_exposure_profile()}.
#'
#' @inheritParams simulate_exposure_profile
#' @return list(profile = numeric, meta = list)
#' @export
simulate_exposure <- function(...) {
  simulate_exposure_profile(...)
}

#' history across lags 0:lag_max. You can generate only a background profile
#' (mode = "profile") or overlay a controlled pattern (mode = "pattern").
#'
#' @param lag_max Integer. Maximum lag (profile length = lag_max + 1).
#' @param mode Character. "profile" (only background) or "pattern" (background + overlay).
#' @param pattern Character. Pattern type if mode="pattern":
#'   "random_lags", "alternating_lags", "block_lags".
#' @param fixed_value Numeric. Value to impose in pattern-controlled lags
#'   (required for "random_lags" and "block_lags").
