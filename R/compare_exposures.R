#' Compare exposure-lag profiles (pairwise or multiple)
#'
#' @param exposures Named list of exposure profiles (preferred)
#' @param exposure1 Legacy single exposure
#' @param exposure2 Legacy second exposure
#' @param lags Optional integer vector (same length as exposures) exp - lags = c(7, 14, 21, 28) is equal to value 1 = lag 7, valor 2 = lag 14...
#' @param q Quantile threshold for high-exposure overlap
#' @param eps Small constant for stability
#'
#' @return data.frame
#' @export
compare_exposures <- function(
    exposures = NULL,
    exposure1 = NULL,
    exposure2 = NULL,
    lags = NULL,
    q = 0.8,
    eps = 1e-12
) {

  # ------------------------------------------------------------
  # ✅ BACKWARD COMPATIBILITY
  # ------------------------------------------------------------
  if (is.null(exposures)) {

    if (is.null(exposure1) || is.null(exposure2)) {
      stop("Provide either `exposures` (list) OR both `exposure1` and `exposure2`.")
    }

    exposures <- list(
      exposure1 = exposure1,
      exposure2 = exposure2
    )
  }

  if (!is.list(exposures) || is.null(names(exposures))) {
    stop("`exposures` must be a named list.")
  }

  n_exp <- length(exposures)
  if (n_exp < 2) {
    stop("At least two exposure profiles are required.")
  }

  # ------------------------------------------------------------
  # ✅ VALIDATIONS
  # ------------------------------------------------------------
  lengths <- sapply(exposures, length)
  if (length(unique(lengths)) != 1) {
    stop("All exposure profiles must have the same length.")
  }

  n <- lengths[1]

  if (n < 2) stop("Exposure profiles must have length >= 2.")

  for (nm in names(exposures)) {
    x <- exposures[[nm]]
    if (!is.numeric(x)) stop(paste("Exposure", nm, "must be numeric."))
    if (any(!is.finite(x))) stop(paste("Exposure", nm, "must contain finite values."))
  }

  if (!is.numeric(q) || length(q) != 1 || q <= 0 || q >= 1) {
    stop("q must be a single number between 0 and 1.")
  }

  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("eps must be a positive numeric scalar.")
  }

  if (is.null(lags)) {
    lags <- 0:(n - 1)
  }

  if (length(lags) != n) {
    stop("lags must have the same length as exposures.")
  }

  lags <- as.numeric(lags)

  # ------------------------------------------------------------
  # ✅ PAIRWISE COMPARISONS
  # ------------------------------------------------------------
  combs <- combn(names(exposures), 2, simplify = FALSE)

  out_list <- lapply(combs, function(cb) {

    name1 <- cb[1]
    name2 <- cb[2]

    exposure1 <- exposures[[name1]]
    exposure2 <- exposures[[name2]]

    d <- exposure1 - exposure2
    absd <- abs(d)

    # distances
    L1 <- sum(absd)
    L2 <- sqrt(sum(d^2))
    Linf <- max(absd)

    # summaries
    sum1 <- sum(exposure1); sum2 <- sum(exposure2)
    mean1 <- mean(exposure1); mean2 <- mean(exposure2)
    sd1 <- stats::sd(exposure1); sd2 <- stats::sd(exposure2)

    # correlation
    corr <- suppressWarnings(stats::cor(exposure1, exposure2))

    # peak timing
    peak_lag1 <- lags[which.max(exposure1)]
    peak_lag2 <- lags[which.max(exposure2)]
    timing_shift_peak <- peak_lag1 - peak_lag2

    # center of mass
    w1 <- pmax(exposure1, 0)
    w2 <- pmax(exposure2, 0)

    com1 <- sum(lags * w1) / (sum(w1) + eps)
    com2 <- sum(lags * w2) / (sum(w2) + eps)

    center_of_mass_shift <- com1 - com2

    # overlap
    thr1 <- stats::quantile(exposure1, probs = q, names = FALSE)
    thr2 <- stats::quantile(exposure2, probs = q, names = FALSE)

    hi1 <- exposure1 >= thr1
    hi2 <- exposure2 >= thr2

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
  out
}
