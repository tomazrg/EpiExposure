#' Compute Exposure Cumulative Impact (ECI)
#'
#' Computes:
#' 1) ECI_raw (sum of exposure values; unweighted)
#' 2) ECI_weighted using DLNM model coefficients
#'
#' @param profile Numeric vector of exposure values (lag profile)
#' @param fit Fitted model from fit_epidlnm() (required for weighted ECI)
#'
#' @return data.frame with ECI_raw and ECI_weighted
#' @export
compute_eci <- function(profile, fit = NULL) {

  # -----------------------
  # Validation
  # -----------------------
  if (!is.numeric(profile)) {
    stop("profile must be a numeric vector.")
  }

  if (any(!is.finite(profile))) {
    stop("profile must contain only finite values.")
  }

  # -----------------------
  # 1) ECI RAW
  # -----------------------
  eci_raw <- sum(profile)

  # -----------------------
  # 2) ECI WEIGHTED (DLNM-based)
  # -----------------------
  eci_weighted <- NA_real_

  if (!is.null(fit)) {

    spec <- attr(fit, "epiexposure_spec")
    vars <- attr(fit, "epiexposure_vars")

    if (is.null(spec)) {
      stop("fit does not contain epiexposure_spec.")
    }

    if (length(vars) != 1) {
      stop("compute_eci currently supports single exposure variable.")
    }

    var <- vars[1]
    spec_v <- spec[[var]]

    if (is.null(spec_v)) {
      stop("Missing spec for variable: ", var)
    }

    lag_max <- as.integer(spec_v$lag_max)

    if (length(profile) != lag_max + 1L) {
      stop("Profile length must be lag_max + 1.")
    }

    # -----------------------
    # Rebuild crossbasis
    # -----------------------
    cb <- dlnm::crossbasis(
      profile,
      lag = lag_max,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )

    cb_row <- as.numeric(cb[lag_max + 1L, ])

    # -----------------------
    # Extract coefficients by engine
    # -----------------------
    beta <- NULL
    model <- fit

    if (inherits(model, "glmmTMB")) {

      beta <- fixef(model)$cond

    } else if (inherits(model, "merMod")) {

      beta <- lme4::fixef(model)

    } else if (inherits(model, "glm") || inherits(model, "gam")) {

      beta <- stats::coef(model)

    } else if (inherits(model, "gls")) {

      beta <- stats::coef(model)

    } else if (inherits(model, "brmsfit")) {

      # brms: usar média posterior dos coeficientes
      beta <- brms::fixef(model)[, "Estimate"]
      names(beta) <- rownames(brms::fixef(model))

    } else if (inherits(model, "inla")) {

      # INLA: usar efeitos fixos
      beta <- model$summary.fixed$mean

    } else {

      stop("Unsupported model class for weighted ECI.")
    }

    # -----------------------
    # Select cb_* coefficients
    # -----------------------
    cb_names <- grep(paste0("^cb_", var, "_"), names(beta), value = TRUE)

    if (length(cb_names) != length(cb_row)) {
      stop("Mismatch between cb coefficients and constructed basis.")
    }

    beta_cb <- beta[cb_names]

    # -----------------------
    # Compute weighted ECI (eta)
    # -----------------------
    eci_weighted <- sum(cb_row * beta_cb)
  }

  # -----------------------
  # OUTPUT
  # -----------------------
  data.frame(
    ECI_raw = eci_raw,
    ECI_weighted = eci_weighted
  )
}
