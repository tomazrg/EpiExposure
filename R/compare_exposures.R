#' Compare exposure-lag profiles (global or lagwise, with optional grouping)
#'
#' Compares exposure-lag profiles between two or more scenarios, epidemics,
#' or simulated exposure histories. The function supports both global
#' (aggregated) comparisons and lag-wise comparisons, enabling detailed
#' temporal analysis of exposure differences.
#'
#' It also supports grouping of exposures (e.g., multiple epidemics),
#' aggregating profiles per lag using a summary function (median or mean),
#' which is particularly useful when working with large observational datasets.
#'
#' @param exposures Named list of numeric vectors, where each element
#'   represents an exposure-lag profile. Each vector must have the same length.
#'   Alternatively, a single numeric vector is accepted and converted to a list.
#' @param exposure1 Legacy numeric vector (optional).
#' @param exposure2 Legacy numeric vector (optional).
#'   Used for backward compatibility if `exposures` is not provided.
#' @param lags Optional numeric vector indicating lag positions
#'   (same length as exposure profiles). Defaults to `0:(n - 1)`.
#' @param reverse Logical. If `TRUE`, each exposure profile is reversed before
#'   analysis. Use this when the input data are ordered from oldest to most
#'   recent (e.g., days before present). If `FALSE`, the function assumes the
#'   first element corresponds to lag 0 (most recent) and the last element to
#'   lag max (oldest).
#' @param group Optional vector defining group membership for each exposure.
#'   Must have the same length as `exposures`. If provided, profiles are first
#'   aggregated within each group and comparisons are performed between groups
#'   instead of individual exposures.
#' @param agg_fun Aggregation function used when `group` is provided.
#'   Options are `"median"` (default, recommended) or `"mean"`.
#' @param mode Type of comparison:
#'   \describe{
#'     \item{"global"}{Aggregated comparison using distances, correlation,
#'     and summary statistics across the full lag profile (default).}
#'     \item{"lagwise"}{Lag-by-lag comparison, returning a long-format data
#'     frame suitable for visualization and temporal analysis.}
#'   }
#' @param q Quantile threshold (between 0 and 1) used to compute overlap of
#'   high-exposure regions (default = 0.8).
#' @param eps Small constant used to avoid division by zero (default = 1e-12).
#'
#' @return A data.frame.
#'
#' If `mode = "global"`, returns one row per pairwise comparison with metrics:
#' \itemize{
#'   \item L1, L2, Linf distances
#'   \item Differences in sum, mean, and standard deviation
#'   \item Correlation between profiles
#'   \item Peak lag timing and shift
#'   \item Center of mass and shift
#'   \item Overlap of high-exposure regions
#' }
#'
#' If `mode = "lagwise"`, returns one row per lag per pair:
#' \itemize{
#'   \item lag
#'   \item value1, value2
#'   \item diff, abs_diff, ratio
#' }
#'
#' @details
#' This function provides two complementary comparison strategies:
#'
#' \strong{Global comparison} summarizes the overall similarity between
#' exposure profiles using distance metrics, correlation, and timing features.
#'
#' \strong{Lag-wise comparison} retains full temporal resolution, allowing
#' identification of when differences between exposure profiles occur,
#' which is particularly useful for plotting and interpretation.
#'
#' When working with observational data (e.g., multiple epidemics), providing
#' a `group` argument is strongly recommended to avoid excessive pairwise
#' comparisons and to obtain representative exposure profiles via aggregation.
#'
#' @export
compare_exposures <- function(
    exposures = NULL,
    exposure1 = NULL,
    exposure2 = NULL,
    lags = NULL,
    reverse = FALSE,
    group = NULL,
    agg_fun = c("median", "mean"),
    mode = c("global", "lagwise"),
    q = 0.8,
    eps = 1e-12
) {

  mode <- match.arg(mode)
  agg_fun <- match.arg(agg_fun)

  # ------------------------------------------------------------
  # BACKWARD COMPATIBILITY
  # ------------------------------------------------------------
  if (is.null(exposures)) {

    if (is.null(exposure1) || is.null(exposure2)) {
      stop("Provide either `exposures` OR both `exposure1` and `exposure2`.")
    }

    exposures <- list(
      exposure1 = exposure1,
      exposure2 = exposure2
    )

  } else if (is.numeric(exposures)) {

    exposures <- list(exposure = exposures)
  }

  if (!is.list(exposures) || is.null(names(exposures))) {
    stop("`exposures` must be a named list.")
  }

  lengths <- sapply(exposures, length)

  if (length(unique(lengths)) != 1) {
    stop("All exposure profiles must have the same length.")
  }

  n <- lengths[1]

  if (n < 2) {
    stop("Exposure profiles must have length >= 2.")
  }

  for (nm in names(exposures)) {
    x <- exposures[[nm]]
    if (!is.numeric(x)) stop(paste("Exposure", nm, "must be numeric."))
    if (any(!is.finite(x))) stop(paste("Exposure", nm, "must contain finite values."))
  }

  if (is.null(lags)) {
    lags <- 0:(n - 1)
  }

  if (length(lags) != n) {
    stop("lags must match exposure length.")
  }

  lags <- as.numeric(lags)

  # ------------------------------------------------------------
  # ORIENTATION
  # ------------------------------------------------------------
  if (reverse) {
    exposures <- lapply(exposures, rev)
    lags <- rev(lags)

    message(
      "Profiles reversed: interpreting input as oldest → most recent. ",
      "After reversal, first value = lag 0 (most recent)."
    )
  } else {
    message(
      "Assuming exposures are lag-ordered: first value = lag 0 (most recent)."
    )
  }

  # ------------------------------------------------------------
  # GROUP AGGREGATION
  # ------------------------------------------------------------
  if (!is.null(group)) {

    if (length(group) != length(exposures)) {
      stop("Length of `group` must match number of exposures.")
    }

    group_levels <- unique(group)

    mat <- do.call(rbind, exposures)

    agg_fun_use <- if (agg_fun == "median") stats::median else base::mean

    exposures <- lapply(group_levels, function(g) {
      rows <- group == g
      apply(mat[rows, , drop = FALSE], 2, agg_fun_use)
    })

    names(exposures) <- group_levels
  }

  combs <- combn(names(exposures), 2, simplify = FALSE)

  # ------------------------------------------------------------
  # LAGWISE MODE
  # ------------------------------------------------------------
  if (mode == "lagwise") {

    out_list <- lapply(combs, function(cb) {

      name1 <- cb[1]
      name2 <- cb[2]

      x1 <- exposures[[name1]]
      x2 <- exposures[[name2]]

      data.frame(
        exposure1 = name1,
        exposure2 = name2,
        lag = lags,
        value1 = x1,
        value2 = x2,
        diff = x1 - x2,
        abs_diff = abs(x1 - x2),
        ratio = x1 / (x2 + eps)
      )
    })

    out <- do.call(rbind, out_list)
    rownames(out) <- NULL
    return(out)
  }

  # ------------------------------------------------------------
  # GLOBAL MODE
  # ------------------------------------------------------------
  out_list <- lapply(combs, function(cb) {

    name1 <- cb[1]
    name2 <- cb[2]

    x1 <- exposures[[name1]]
    x2 <- exposures[[name2]]

    d <- x1 - x2
    absd <- abs(d)

    L1 <- sum(absd)
    L2 <- sqrt(sum(d^2))
    Linf <- max(absd)

    sum1 <- sum(x1); sum2 <- sum(x2)
    mean1 <- mean(x1); mean2 <- mean(x2)

    sd1 <- stats::sd(x1)
    sd2 <- stats::sd(x2)

    corr <- suppressWarnings(stats::cor(x1, x2))

    peak_lag1 <- lags[which.max(x1)]
    peak_lag2 <- lags[which.max(x2)]

    timing_shift_peak <- peak_lag1 - peak_lag2

    w1 <- pmax(x1, 0)
    w2 <- pmax(x2, 0)

    com1 <- sum(lags * w1) / (sum(w1) + eps)
    com2 <- sum(lags * w2) / (sum(w2) + eps)

    center_of_mass_shift <- com1 - com2

    thr1 <- stats::quantile(x1, probs = q, names = FALSE)
    thr2 <- stats::quantile(x2, probs = q, names = FALSE)

    hi1 <- x1 >= thr1
    hi2 <- x2 >= thr2

    overlap_above_q <- sum(hi1 & hi2) / (sum(hi1 | hi2) + eps)

    data.frame(
      exposure1 = name1,
      exposure2 = name2,
      n_lags = n,
      L1 = L1,
      L2 = L2,
      Linf = Linf,
      diff_sum = sum1 - sum2,
      diff_mean = mean1 - mean2,
      diff_sd = sd1 - sd2,
      corr = corr,
      peak_lag1 = peak_lag1,
      peak_lag2 = peak_lag2,
      timing_shift_peak = timing_shift_peak,
      com1 = com1,
      com2 = com2,
      center_of_mass_shift = center_of_mass_shift,
      overlap_above_q = overlap_above_q,
      q = q
    )
  })

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  return(out)
}
