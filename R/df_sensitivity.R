#' Compute sensitivity of DLNM effects (analytical or finite derivative)
#'
#' Calculates the derivative of the exposure–response relationship from
#' predicted values or effects. The function (via GAM smoothing) and finite differences, with optional#' predicted values or effects. The function supports both analytical
#' elasticity computation and critical point detection.
#'
#' @param df Data frame containing at least columns `value` (predictor)
#'   and `prediction` or `effect` (response).
#' @param x Character. Name of predictor variable (default = "value")
#' @param y Character. Name of response variable (default = "prediction")
#' @param scenario_var Optional character. Grouping variable (e.g., scenario or id).
#'   If NULL, computation is done on the full dataset.
#' @param k Integer. Basis dimension for GAM smoothing (default = 10)
#' @param method Character. Derivative method:
#'   \describe{
#'     \item{"analytical"}{Uses GAM and lpmatrix to compute smooth derivative (default).}
#'     \item{"finite"}{Uses finite differences (robust but noisier).}
#'   }
#' @param elasticity Logical. If TRUE, computes elasticity:
#'   \deqn{(dy/dx) * (x/y)}
#' @param critical Logical. If TRUE, detects critical points (sign change in derivative)
#' @param eps Small numeric value for numerical stability (default = 1e-6)
#'
#' @return A data.frame including:
#' \itemize{
#'   \item `sensitivity` (dy/dx)
#'   \item `elasticity` (optional)
#'   \item `critical` (TRUE/FALSE indicator of critical points)
#' }
#'
#' @details
#' The analytical method fits a smooth GAM and computes derivatives via
#' the linear predictor matrix (`lpmatrix`). The finite method computes
#' local differences between adjacent observations.
#'
#' Analytical derivatives are smoother and preferred for continuous
#' exposure–response curves, whereas finite differences are more robust
#' when the data are irregular or sparse.
#'
#' @export
df_sensitivity <- function(df,
                           x = "value",
                           y = "prediction",
                           scenario_var = "scenario",
                           k = 10,
                           method = c("analytical", "finite"),
                           elasticity = TRUE,
                           critical = TRUE,
                           eps = 1e-6) {

  method <- match.arg(method)

  # ------------------------------------------------------------
  # ✅ VALIDATIONS
  # ------------------------------------------------------------
  if (!x %in% names(df)) stop(paste("Column", x, "not found in df."))
  if (!y %in% names(df)) stop(paste("Column", y, "not found in df."))

  if (!is.null(scenario_var) && !scenario_var %in% names(df)) {
    stop(paste("scenario_var", scenario_var, "not found in df."))
  }

  if (!is.numeric(k) || k <= 0) stop("k must be a positive integer.")
  k <- as.integer(k)

  if (!is.numeric(eps) || eps <= 0) {
    stop("eps must be a positive numeric scalar.")
  }

  # ------------------------------------------------------------
  # ✅ Analytical derivative via GAM
  # ------------------------------------------------------------
  compute_derivative_analytical <- function(d) {

    d <- d[order(d[[x]]), ]

    dx_check <- diff(d[[x]])
    if (any(dx_check <= 0, na.rm = TRUE)) {
      warning(
        "Non-strictly increasing '", x,
        "' detected. Results may be unstable."
      )
    }

    form <- as.formula(paste0(y, " ~ s(", x, ", k=", k, ", bs='cs')"))
    fit  <- mgcv::gam(form, data = d)

    Xp <- predict(fit, newdata = d, type = "lpmatrix")

    d_eps <- d
    d_eps[[x]] <- d_eps[[x]] + eps

    Xp_eps <- predict(fit, newdata = d_eps, type = "lpmatrix")

    dXp <- (Xp_eps - Xp) / eps
    dy  <- as.vector(dXp %*% coef(fit))

    d$sensitivity <- dy
    d
  }

  # ------------------------------------------------------------
  # ✅ Finite difference
  # ------------------------------------------------------------
  compute_derivative_finite <- function(d) {

    d <- d[order(d[[x]]), ]

    dx <- diff(d[[x]])
    dy <- diff(d[[y]])

    sens <- rep(NA_real_, length(dx))
    valid <- abs(dx) > eps

    sens[valid] <- dy[valid] / dx[valid]

    d$sensitivity <- c(NA, sens)
    d
  }

  # ------------------------------------------------------------
  # ✅ Combined processing
  # ------------------------------------------------------------
  compute_all <- function(d) {

    d <- d[order(d[[x]]), ]

    if (method == "analytical") {
      d <- compute_derivative_analytical(d)
    } else {
      d <- compute_derivative_finite(d)
    }

    # Elasticity
    if (elasticity) {
      d$elasticity <- ifelse(
        abs(d[[y]]) > eps,
        d$sensitivity * (d[[x]] / d[[y]]),
        NA_real_
      )
    }

    # Critical points
    if (critical) {

      s <- d$sensitivity

      s_sign <- sign(s)
      s_sign[!is.finite(s_sign)] <- NA

      sc <- c(NA, diff(s_sign))

      crit_idx <- which(!is.na(sc) & sc != 0)

      d$critical <- FALSE
      d$critical[crit_idx] <- TRUE
    }

    d
  }

  # ------------------------------------------------------------
  # ✅ Apply per scenario
  # ------------------------------------------------------------
  if (!is.null(scenario_var)) {

    split_list <- split(df, df[[scenario_var]])
    res_list   <- lapply(split_list, compute_all)

    out <- do.call(rbind, res_list)

    out <- out[order(out[[scenario_var]], out[[x]]), ]

  } else {

    out <- compute_all(df)

    out <- out[order(out[[x]]), ]
  }

  rownames(out) <- NULL

  return(out)
}
