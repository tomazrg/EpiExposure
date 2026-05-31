#' Compute ECI decomposition by lag (ECI-by-lag)
#'
#' Returns ECI_raw and ECI_weighted, plus a per-lag decomposition of the
#' weighted ECI using numerical derivatives on the DLNM linear predictor.
#'
#' Idea:
#'   eta(x) = sum_k beta_k * B_k(x)   (DLNM linear predictor contribution)
#'   w_l(x) ≈ [eta(x + eps*e_l) - eta(x)] / eps
#'   contribution_l = x_l * w_l
#'
#' @param profile Numeric vector, length = lag_max + 1, ordered as lag 0..lag_max
#' @param fit Fitted model returned by fit_epidlnm() (stores epiexposure_spec attrs)
#' @param eps Small numeric perturbation for finite differences (default 1e-6)
#' @param center If TRUE, uses central difference; else forward difference
#' @param absolute If TRUE, report absolute contributions as well
#'
#' @return list with:
#'   - ECI_raw
#'   - ECI_weighted (eta)
#'   - by_lag (data.frame)
#' @export
compute_eci_by_lag <- function(profile, fit, eps = 1e-6, center = TRUE, absolute = TRUE) {
  
  # -----------------------
  # Validations
  # -----------------------
  if (!is.numeric(profile) || any(!is.finite(profile))) stop("profile must be a finite numeric vector.")
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) stop("eps must be a positive numeric scalar.")
  if (!is.logical(center) || length(center) != 1) stop("center must be TRUE/FALSE.")
  if (!is.logical(absolute) || length(absolute) != 1) stop("absolute must be TRUE/FALSE.")
  
  spec <- attr(fit, "epiexposure_spec")
  vars <- attr(fit, "epiexposure_vars")
  
  if (is.null(spec) || !is.list(spec)) stop("fit is missing epiexposure_spec attribute.")
  if (is.null(vars) || length(vars) != 1) stop("compute_eci_by_lag currently supports a single exposure variable fit.")
  
  var <- vars[1]
  spec_v <- spec[[var]]
  if (is.null(spec_v)) stop("Missing spec for variable: ", var)
  
  lag_max <- as.integer(spec_v$lag_max)
  if (length(profile) != lag_max + 1L) stop("Profile length must be lag_max + 1.")
  
  # -----------------------
  # Helper: extract beta_cb from supported engines
  # -----------------------
  extract_beta <- function(model) {
    if (inherits(model, "glmmTMB")) {
      b <- glmmTMB::fixef(model)$cond
      return(b)
    }
    if (inherits(model, "merMod")) {
      b <- lme4::fixef(model)
      return(b)
    }
    if (inherits(model, "glm") || inherits(model, "gam")) {
      b <- stats::coef(model)
      return(b)
    }
    if (inherits(model, "gls")) {
      b <- stats::coef(model)
      return(b)
    }
    if (inherits(model, "brmsfit")) {
      b <- brms::fixef(model)[, "Estimate"]
      names(b) <- rownames(brms::fixef(model))
      return(b)
    }
    if (inherits(model, "inla")) {
      b <- model$summary.fixed$mean
      return(b)
    }
    stop("Unsupported model class for ECI-by-lag decomposition.")
  }
  
  model <- fit
  beta <- extract_beta(model)
  
  cb_names <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)
  if (length(cb_names) == 0) stop("No cb_ coefficients found for variable: ", var)
  
  beta_cb <- beta[cb_names]
  
  # -----------------------
  # Helper: compute eta(profile) = sum(cb_row * beta_cb)
  # -----------------------
  eta_from_profile <- function(p) {
    cb <- dlnm::crossbasis(
      p,
      lag = lag_max,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )
    cb_row <- as.numeric(cb[lag_max + 1L, ])
    if (length(cb_row) != length(beta_cb)) {
      stop("Mismatch between crossbasis dimension and cb coefficients. Spec/model inconsistent.")
    }
    sum(cb_row * beta_cb)
  }
  
  # -----------------------
  # Compute ECI_raw and ECI_weighted (eta)
  # -----------------------
  ECI_raw <- sum(profile)
  eta0 <- eta_from_profile(profile)
  ECI_weighted <- eta0
  
  # -----------------------
  # Compute lag weights w_l via finite differences
  # -----------------------
  lags <- 0:lag_max
  w <- numeric(length(lags))
  
  if (center) {
    for (i in seq_along(lags)) {
      p_plus <- profile; p_minus <- profile
      p_plus[i]  <- p_plus[i]  + eps
      p_minus[i] <- p_minus[i] - eps
      w[i] <- (eta_from_profile(p_plus) - eta_from_profile(p_minus)) / (2 * eps)
    }
  } else {
    for (i in seq_along(lags)) {
      p_plus <- profile
      p_plus[i] <- p_plus[i] + eps
      w[i] <- (eta_from_profile(p_plus) - eta0) / eps
    }
  }
  
  contribution <- profile * w
  
  by_lag <- data.frame(
    var = var,
    lag = lags,
    exposure = profile,
    weight = w,
    contribution = contribution
  )
  
  if (absolute) {
    by_lag$abs_contribution <- abs(by_lag$contribution)
    by_lag$abs_weight <- abs(by_lag$weight)
  }
  
  # percent contribution (based on sum of abs contributions for stability)
  denom <- sum(abs(by_lag$contribution))
  by_lag$percent_contribution <- if (denom > 0) 100 * abs(by_lag$contribution) / denom else NA_real_
  
  list(
    ECI_raw = ECI_raw,
    ECI_weighted = ECI_weighted,
    by_lag = by_lag
  )
}

#' Identify critical lags from ECI-by-lag decomposition
#'
#' @param eci_by_lag Output from compute_eci_by_lag()
#' @param top_n Number of critical lags to return
#' @param metric Which metric to rank by: "abs_contribution" or "abs_weight" or "percent_contribution"
#'
#' @return data.frame of top critical lags
#' @export
identify_critical_lags_eci <- function(eci_by_lag, top_n = 10, metric = c("abs_contribution", "percent_contribution", "abs_weight")) {
  
  metric <- match.arg(metric)
  if (!is.list(eci_by_lag) || is.null(eci_by_lag$by_lag)) stop("eci_by_lag must be output from compute_eci_by_lag().")
  df <- eci_by_lag$by_lag
  
  if (!metric %in% names(df)) stop("Requested metric not found in by_lag output: ", metric)
  if (!is.numeric(top_n) || length(top_n) != 1 || top_n <= 0) stop("top_n must be a positive integer.")
  top_n <- as.integer(top_n)
  
  df <- df[order(df[[metric]], decreasing = TRUE), , drop = FALSE]
  head(df, top_n)
}

#' Plot ECI-by-lag decomposition
#'
#' @param eci_by_lag Output from compute_eci_by_lag()
#' @param what What to plot: "contribution" or "weight" or "percent_contribution"
#'
#' @return ggplot object
#' @export
plot_eci_by_lag <- function(eci_by_lag, what = c("contribution", "weight", "percent_contribution")) {
  
  what <- match.arg(what)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  
  df <- eci_by_lag$by_lag
  ylab <- switch(
    what,
    contribution = "ECI contribution (x_l × w_l)",
    weight = "Lag weight (d eta / d x_l)",
    percent_contribution = "Contribution (%)"
  )
  
  ggplot2::ggplot(df, ggplot2::aes(x = lag, y = .data[[what]])) +
    ggplot2::geom_col(fill = "#636363") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = paste0("ECI-by-lag: ", what),
      x = "Lag (days)",
      y = ylab
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
}