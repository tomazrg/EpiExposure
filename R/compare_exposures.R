#' Compare two exposure-lag profiles
#'
#' @param exposure1 Numeric vector (lag profile)
#' @param exposure2 Numeric vector (lag profile)
#' @param lags Optional integer vector of same length as exposures (default 0:(n-1))
#' @param q Quantile used to define "high exposure" regions for overlap (default 0.8)
#' @param eps Small constant to avoid division by zero (default 1e-12)
#'
#' @return data.frame with comparison metrics
#' @export
compare_exposures <- function(exposure1, exposure2, lags = NULL, q = 0.8, eps = 1e-12) {
  
  # ---- validations ----
  if (!is.numeric(exposure1) || !is.numeric(exposure2)) stop("Both exposures must be numeric vectors.")
  if (length(exposure1) != length(exposure2)) stop("Exposure profiles must have the same length.")
  if (length(exposure1) < 2) stop("Exposure profiles must have length >= 2.")
  if (any(!is.finite(exposure1)) || any(!is.finite(exposure2))) stop("Exposure profiles must contain only finite values.")
  if (!is.numeric(q) || length(q) != 1 || !is.finite(q) || q <= 0 || q >= 1) stop("q must be a single number between 0 and 1.")
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) stop("eps must be a positive numeric scalar.")
  
  n <- length(exposure1)
  if (is.null(lags)) lags <- 0:(n - 1)
  if (length(lags) != n) stop("lags must have the same length as exposures.")
  lags <- as.numeric(lags)
  
  d <- exposure1 - exposure2
  absd <- abs(d)
  
  # Distances
  L1 <- sum(absd)
  L2 <- sqrt(sum(d^2))
  Linf <- max(absd)
  
  # Basic summaries
  sum1 <- sum(exposure1); sum2 <- sum(exposure2)
  mean1 <- mean(exposure1); mean2 <- mean(exposure2)
  sd1 <- stats::sd(exposure1); sd2 <- stats::sd(exposure2)
  
  # Correlation (shape similarity)
  corr <- suppressWarnings(stats::cor(exposure1, exposure2))
  
  # Peak timing
  peak_lag1 <- lags[which.max(exposure1)]
  peak_lag2 <- lags[which.max(exposure2)]
  timing_shift_peak <- peak_lag1 - peak_lag2
  
  # Center of mass (timing of exposure mass)
  w1 <- pmax(exposure1, 0)
  w2 <- pmax(exposure2, 0)
  com1 <- sum(lags * w1) / (sum(w1) + eps)
  com2 <- sum(lags * w2) / (sum(w2) + eps)
  center_of_mass_shift <- com1 - com2
  
  # Overlap in "high exposure" region (Jaccard)
  thr1 <- stats::quantile(exposure1, probs = q, names = FALSE, type = 7)
  thr2 <- stats::quantile(exposure2, probs = q, names = FALSE, type = 7)
  hi1 <- exposure1 >= thr1
  hi2 <- exposure2 >= thr2
  overlap_above_q <- sum(hi1 & hi2) / (sum(hi1 | hi2) + eps)
  
  data.frame(
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
}