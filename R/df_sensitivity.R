#' Compute sensitivity of exposure–response curves
#'
#' Calculates the derivative of the exposure–response relationship from
#' predicted values or effects.
#'
#' The function supports both analytical derivatives obtained through
#' GAM smoothing and finite differences, with optional elasticity
#' computation and critical-point detection.
#'
#' @param data Data frame containing at least columns `value`
#'   (predictor) and `prediction` or `effect` (response).
#' @param x Character. Name of predictor variable
#'   (default = `"value"`).
#' @param y Character. Name of response variable
#'   (default = `"prediction"`).
#' @param scenario_var Optional character. Grouping variable
#'   (e.g. scenario or id). If `NULL`, computation is performed on the
#'   entire dataset.
#' @param k Positive integer basis dimension used for GAM smoothing.
#'   Default = `10`.
#' @param method Character.
#'   \describe{
#'     \item{"analytical"}{
#'       Computes derivatives from a GAM fitted to the curve.
#'     }
#'     \item{"finite"}{
#'       Computes finite-difference derivatives.
#'     }
#'   }
#' @param elasticity Logical.
#'   If `TRUE`, computes elasticity:
#'   \deqn{
#'   (dy/dx)(x/y)
#'   }
#' @param critical Logical.
#'   If `TRUE`, detects critical points based on derivative sign changes.
#' @param eps Small positive numeric constant used for numerical stability.
#'
#' @return A data.frame including:
#' \itemize{
#'   \item `sensitivity`
#'   \item `elasticity` (optional)
#'   \item `critical` (optional)
#' }
#'
#' @details
#' The analytical method fits a smooth GAM and computes derivatives
#' using the linear-predictor matrix (`lpmatrix`).
#'
#' The finite method computes local differences between adjacent
#' observations.
#'
#' Analytical derivatives are smoother and generally preferred for
#' continuous exposure–response curves, whereas finite differences are
#' more robust when curves are sparse or irregular.
#'
#' @export
df_sensitivity <- function(
    data,
    x = "value",
    y = "prediction",
    scenario_var = "scenario",
    k = 10,
    method = c("analytical", "finite"),
    elasticity = TRUE,
    critical = TRUE,
    eps = 1e-6
) {

  method <- match.arg(method)

  # ============================================================
  # VALIDATION
  # ============================================================

  if (!is.data.frame(data)) {
    stop(
      "`data` must be a data.frame."
    )
  }

  if (nrow(data) < 2L) {
    stop(
      "At least two observations are required."
    )
  }

  if (!x %in% names(data)) {
    stop(
      paste("Column", x, "not found in `data`.")
    )
  }

  if (!y %in% names(data)) {
    stop(
      paste("Column", y, "not found in `data`.")
    )
  }

  if (
    any(!is.finite(data[[x]]))
  ) {
    stop(
      paste0(
        "`",
        x,
        "` contains non-finite values."
      )
    )
  }

  if (
    any(!is.finite(data[[y]]))
  ) {
    stop(
      paste0(
        "`",
        y,
        "` contains non-finite values."
      )
    )
  }

  if (
    !is.null(scenario_var) &&
    !scenario_var %in% names(data)
  ) {
    stop(
      paste(
        "scenario_var",
        scenario_var,
        "not found in `data`."
      )
    )
  }

  if (
    !is.numeric(k) ||
    length(k) != 1L ||
    !is.finite(k) ||
    k <= 0
  ) {
    stop(
      "`k` must be a positive integer."
    )
  }

  k <- as.integer(k)

  if (
    !is.numeric(eps) ||
    length(eps) != 1L ||
    !is.finite(eps) ||
    eps <= 0
  ) {
    stop(
      "`eps` must be a positive numeric scalar."
    )
  }

  # ============================================================
  # ANALYTICAL DERIVATIVE
  # ============================================================

  compute_derivative_analytical <- function(d) {

    d <- d[
      order(d[[x]]),
      ,
      drop = FALSE
    ]

    dx_check <- diff(d[[x]])

    if (
      any(dx_check <= 0, na.rm = TRUE)
    ) {
      warning(
        paste0(
          "Non-strictly increasing values detected in `",
          x,
          "`. Results may be unstable."
        ),
        call. = FALSE
      )
    }

    k_use <- min(
      k,
      max(
        3L,
        nrow(d) - 1L
      )
    )

    form <- stats::as.formula(
      paste0(
        y,
        " ~ s(",
        x,
        ", k = ",
        k_use,
        ", bs = 'cs')"
      )
    )

    fit <- mgcv::gam(
      form,
      data = d
    )

    Xp <- predict(
      fit,
      newdata = d,
      type = "lpmatrix"
    )

    d_eps <- d
    d_eps[[x]] <- d_eps[[x]] + eps

    Xp_eps <- predict(
      fit,
      newdata = d_eps,
      type = "lpmatrix"
    )

    dXp <- (Xp_eps - Xp) / eps

    d$sensitivity <- as.vector(
      dXp %*% stats::coef(fit)
    )

    d
  }

  # ============================================================
  # FINITE DIFFERENCE
  # ============================================================

  compute_derivative_finite <- function(d) {

    d <- d[
      order(d[[x]]),
      ,
      drop = FALSE
    ]

    dx <- diff(d[[x]])
    dy <- diff(d[[y]])

    sensitivity <- rep(
      NA_real_,
      length(dx)
    )

    valid <- abs(dx) > eps

    sensitivity[valid] <- dy[valid] / dx[valid]

    d$sensitivity <- c(
      NA_real_,
      sensitivity
    )

    d
  }

  # ============================================================
  # MAIN CALCULATION
  # ============================================================

  compute_all <- function(d) {

    d <- d[
      order(d[[x]]),
      ,
      drop = FALSE
    ]

    if (method == "analytical") {

      d <- compute_derivative_analytical(d)

    } else {

      d <- compute_derivative_finite(d)

    }

    if (elasticity) {

      d$elasticity <- ifelse(
        abs(d[[y]]) > eps,
        d$sensitivity *
          (d[[x]] / d[[y]]),
        NA_real_
      )
    }

    if (critical) {

      s <- d$sensitivity

      s_sign <- sign(s)

      s_sign[
        !is.finite(s_sign)
      ] <- NA

      sign_change <- c(
        NA,
        diff(s_sign)
      )

      critical_index <- which(
        !is.na(sign_change) &
          sign_change != 0
      )

      d$critical <- FALSE

      d$critical[
        critical_index
      ] <- TRUE
    }

    d
  }

  # ============================================================
  # APPLY BY SCENARIO
  # ============================================================

  if (!is.null(scenario_var)) {

    split_list <- split(
      data,
      data[[scenario_var]]
    )

    result_list <- lapply(
      split_list,
      compute_all
    )

    out <- do.call(
      rbind,
      result_list
    )

    out <- out[
      order(
        out[[scenario_var]],
        out[[x]]
      ),
      ,
      drop = FALSE
    ]

  } else {

    out <- compute_all(data)

    out <- out[
      order(out[[x]]),
      ,
      drop = FALSE
    ]
  }

  rownames(out) <- NULL

  out
}
