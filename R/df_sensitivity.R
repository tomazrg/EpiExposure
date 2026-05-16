#' Compute sensitivity of DLNM effects (analytical or finite derivative)
#'
#' Calculates the derivative of the exposure–response relationship,
#' optionally including elasticity and critical point detection.
#'
#' @param df Data frame with columns: value and prediction/effect,
#'   and optionally a scenario column
#' @param x Name of predictor variable (default = "value")
#' @param y Name of response variable (default = "prediction")
#' @param scenario_var Optional grouping variable (default = "scenario")
#' @param k Basis dimension for GAM smoothing (default = 10)
#' @param method Derivative method: "analytical" (default) or "finite"
#' @param elasticity Logical, compute elasticity (default = TRUE)
#' @param critical Logical, detect critical points (default = TRUE)
#' @param eps Small epsilon for numerical derivative (default = 1e-6)
#'
#' @return data.frame with sensitivity, optional elasticity and critical points
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

  stopifnot(x %in% names(df), y %in% names(df))

  if (!is.null(scenario_var)) {
    stopifnot(scenario_var %in% names(df))
  }

  # ------------------------------------------------------------
  # ✅ Analytical derivative via GAM (lpmatrix)
  # ------------------------------------------------------------
  compute_derivative_analytical <- function(d) {

    d <- d[order(d[[x]]), ]

    form <- as.formula(paste0(y, " ~ s(", x, ", k=", k, ", bs='cs')"))
    fit  <- mgcv::gam(form, data = d)

    # base prediction matrix
    Xp <- predict(fit, newdata = d, type = "lpmatrix")

    # perturbed matrix
    d_eps <- d
    d_eps[[x]] <- d_eps[[x]] + eps

    Xp_eps <- predict(fit, newdata = d_eps, type = "lpmatrix")

    # derivative of basis
    dXp <- (Xp_eps - Xp) / eps

    # derivative
    dy <- as.vector(dXp %*% coef(fit))

    d$sensitivity <- dy
    d
  }

  # ------------------------------------------------------------
  # ✅ Finite difference (simple fallback)
  # ------------------------------------------------------------
  compute_derivative_finite <- function(d) {

    d <- d[order(d[[x]]), ]

    dx <- diff(d[[x]])
    dy <- diff(d[[y]])

    d$sensitivity <- c(NA, dy / dx)
    d
  }

  # ------------------------------------------------------------
  # ✅ Apply derivative
  # ------------------------------------------------------------
  compute_all <- function(d) {

    if (method == "analytical") {
      d <- compute_derivative_analytical(d)
    } else {
      d <- compute_derivative_finite(d)
    }

    # --------------------------------------------------------
    # ✅ Elasticity
    # --------------------------------------------------------
    if (elasticity) {

      d$elasticity <- with(d,
                           ifelse(abs(d[[y]]) > 1e-8,
                                  sensitivity * (d[[x]] / d[[y]]),
                                  NA)
      )
    }

    # --------------------------------------------------------
    # ✅ Critical points detection
    # --------------------------------------------------------
    if (critical) {

      s <- d$sensitivity

      # sign change detection
      sc <- c(NA, diff(sign(s)))

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

    out <- split(df, df[[scenario_var]]) |>
      lapply(compute_all) |>
      do.call(rbind, _)

  } else {

    out <- compute_all(df)
  }

  rownames(out) <- NULL

  out
}
