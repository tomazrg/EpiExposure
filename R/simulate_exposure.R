#' Simulate exposure history across lags (full or patterned)
#'
#' Creates one or multiple exposure-lag profiles of length (lag_max + 1)
#'
#' @param lag_max Integer. Maximum lag (profile length = lag_max + 1)
#' @param n Integer. Number of simulations (default = 1)
#' @param mode Character. "profile" or "pattern"
#' @param pattern Character. Pattern type
#' @param fixed_value Numeric. Value for pattern-controlled lags
#' @param n_lags Integer for "random_lags"
#' @param lag_range Integer vector of lags
#' @param alternating_values Numeric vector (>=2)
#' @param block_lags Integer length-2
#' @param background List controlling background
#' @param bounds Optional numeric length-2
#' @param seed Optional seed
#' @param cumulative Logical. If TRUE, returns cumulative (cumsum) profile
#'
#' @return
#' - If n = 1: list(profile, meta)
#' - If n > 1: list(profiles, meta)
#' @export
simulate_exposure <- function(
    lag_max,
    n = 1,
    mode = c("profile", "pattern"),
    pattern = c("random_lags", "alternating_lags", "block_lags"),
    fixed_value = NULL,
    n_lags = NULL,
    lag_range = NULL,
    alternating_values = NULL,
    block_lags = NULL,
    background = list(dist = "normal", mean = 0, sd = 1),
    bounds = NULL,
    seed = NULL,
    cumulative = FALSE
) {

  # ---- helpers ----
  nullcoalesce <- function(a, b) if (!is.null(a)) a else b
  is_whole_number <- function(x) is.numeric(x) && length(x) == 1 &&
    is.finite(x) && abs(x - round(x)) < .Machine$double.eps^0.5

  # ---- validate n ----
  if (!is_whole_number(n) || n <= 0) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(round(n))

  # ---- internal function ----
  simulate_one <- function() {

    if (!is_whole_number(lag_max) || lag_max < 0) {
      stop("lag_max must be a non-negative integer.")
    }
    lag_max_i <- as.integer(round(lag_max))

    mode_i <- match.arg(mode)
    pattern_i <- match.arg(pattern)

    if (!is.null(seed)) {
      if (!is_whole_number(seed)) stop("seed must be an integer if provided.")
      set.seed(as.integer(round(seed)))
    }

    N <- lag_max_i + 1L
    lags_all <- 0:lag_max_i

    if (!is.null(bounds)) {
      if (!is.numeric(bounds) || length(bounds) != 2L || any(!is.finite(bounds))) {
        stop("bounds must be numeric length-2.")
      }
      if (bounds[1] > bounds[2]) stop("bounds[1] must be <= bounds[2].")
    }

    if (!is.list(background)) stop("background must be a list.")
    dist <- nullcoalesce(background$dist, "normal")

    # ---- background ----
    gen_background <- function() {
      if (dist == "normal") {
        mu <- nullcoalesce(background$mean, 0)
        sd <- nullcoalesce(background$sd, 1)
        stats::rnorm(N, mu, sd)

      } else if (dist == "empirical") {
        vals <- background$values
        vals <- vals[is.finite(vals)]
        sample(vals, size = N, replace = TRUE)

      } else if (dist == "ar1") {
        mu <- nullcoalesce(background$mean, 0)
        sd <- nullcoalesce(background$sd, 1)
        phi <- nullcoalesce(background$phi, 0.7)

        x <- numeric(N)
        x[1] <- stats::rnorm(1, mu, sd)

        eps_sd <- if (sd == 0) 0 else sd * sqrt(1 - phi^2)

        for (i in 2:N) {
          x[i] <- mu + phi * (x[i - 1] - mu) +
            stats::rnorm(1, 0, eps_sd)
        }
        x

      } else if (dist == "fixed") {
        rep(background$value, N)

      } else {
        stop("Unknown background dist.")
      }
    }

    x <- gen_background()

    if (!is.null(bounds)) {
      x <- pmin(pmax(x, bounds[1]), bounds[2])
    }

    fixed_lags <- integer(0)

    # ---- pattern ----
    if (mode_i == "pattern") {

      lag_range_i <- if (is.null(lag_range)) lags_all else intersect(lag_range, lags_all)

      if (pattern_i == "random_lags") {

        n_lags_i <- min(n_lags, length(lag_range_i))
        fixed_lags <- sort(sample(lag_range_i, n_lags_i))
        x[fixed_lags + 1L] <- fixed_value

      } else if (pattern_i == "alternating_lags") {

        fixed_lags <- sort(lag_range_i)
        alt <- rep(alternating_values, length.out = length(fixed_lags))
        x[fixed_lags + 1L] <- alt

      } else if (pattern_i == "block_lags") {

        a <- max(0L, min(block_lags))
        b <- min(lag_max_i, max(block_lags))
        fixed_lags <- a:b
        x[fixed_lags + 1L] <- fixed_value
      }
    }

    # ---- bounds AFTER pattern ----
    if (!is.null(bounds)) {
      x <- pmin(pmax(x, bounds[1]), bounds[2])
    }

    # If the variable needs to be cumulative across the lags.
    if (isTRUE(cumulative)) {
      x <- cumsum(x)
    }

    names(x) <- paste0("lag_", lags_all)

    list(
      profile = x,
      meta = list(
        lag_max = lag_max_i,
        mode = mode_i,
        pattern = if (mode_i == "pattern") pattern_i else NA_character_,
        background = background,
        bounds = bounds,
        seed = seed,
        fixed_lags = fixed_lags,
        cumulative = cumulative
      )
    )
  }

  # ---- run simulations ----
  sims <- replicate(n, simulate_one(), simplify = FALSE)

  # ---- backward compatibility ----
  if (n == 1) {
    return(sims[[1]])
  }

  profiles <- lapply(sims, function(s) s$profile)
  meta <- lapply(sims, function(s) s$meta)

  return(list(
    profiles = profiles,
    meta = meta
  ))
}
