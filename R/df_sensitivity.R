#' Compute sensitivity of exposure-response curves
#'
#' Estimates the local derivative of an already calculated exposure-response
#' curve. The function can work directly with supplied curve values or first
#' smooth the curve with a GAM and differentiate the fitted smooth.
#'
#' `df_sensitivity()` is a **post-processing** function. It does not fit an
#' EpiExposure DLNM, reconstruct exposure histories, or redefine the fitted lag
#' window. Consequently, rows supplied to this function are curve-evaluation
#' points, not necessarily the `max_lag + 1` observations that constituted an
#' original exposure history.
#'
#' @param data Non-empty data frame containing the curve to differentiate.
#'   At minimum, it must contain the columns selected by `x` and `y`.
#'
#'   Within each group defined by `scenario_var`, `x` must contain unique finite
#'   numeric values. Duplicate `x` values are rejected rather than aggregated
#'   silently because a derivative requires one unambiguous curve value at each
#'   evaluation position.
#' @param x Character scalar naming the predictor/evaluation coordinate.
#'   Default is `"value"`.
#' @param y Character scalar naming the curve value to differentiate.
#'   Default is `"prediction"`.
#' @param scenario_var `NULL` or a character vector naming one or more grouping
#'   columns. The derivative is computed independently within every unique
#'   combination of these columns.
#'
#'   The historical default is `"scenario"`. Set `scenario_var = NULL` to
#'   differentiate the entire supplied data set as one curve. If several
#'   exposure variables, scenarios, or parameter draws are present, include all
#'   columns needed to identify one unique curve, for example
#'   `scenario_var = c("var", "scenario")` or
#'   `scenario_var = c("scenario", "sample")`.
#' @param k Positive integer GAM basis dimension requested when
#'   `method = "analytical"`. Default is `10`.
#'
#'   The effective value is reduced, when necessary, to the number of distinct
#'   `x` values in the current curve. It is never increased silently. A
#'   non-constant curve analyzed with `method = "analytical"` requires at least
#'   three distinct `x` values and `k >= 3`.
#' @param method Character. One of:
#'
#'   - `"analytical"`: historical EpiExposure name for a GAM-based smooth
#'     derivative. A univariate GAM is fitted with
#'     `mgcv::gam(..., method = "REML")`, using the auxiliary smooth basis
#'     selected by `smooth_basis`, then its linear-predictor matrix is
#'     differentiated numerically using a small within-range perturbation;
#'   - `"finite"`: differentiates the supplied curve directly using local
#'     finite-difference formulas. With at least three points, second-order
#'     three-point formulas are used and support irregular `x` spacing.
#'
#'   The label `"analytical"` is retained for backward compatibility, but this
#'   method is not a symbolic analytical derivative. It is a derivative of the
#'   fitted GAM smooth obtained by finite differencing its `lpmatrix`.
#' @param smooth_basis Character. Auxiliary `mgcv` smooth basis used only when
#'   `method = "analytical"`. One of:
#'
#'   - `"cs"` (default): shrinkage cubic regression spline;
#'   - `"cr"`: cubic regression spline without the additional shrinkage
#'     modification;
#'   - `"tp"`: thin plate regression spline;
#'   - `"ts"`: shrinkage thin plate regression spline.
#'
#'   This choice belongs only to the **post-processing GAM fitted inside
#'   `df_sensitivity()`**. It does **not** redefine, replace, inherit, or need to
#'   match the exposure-response or lag-response basis used by
#'   `define_exposures()`, `fit_epidlnm()`, or `find_bestfit()`. The DLNM has
#'   already been fitted before this function is called.
#'
#'   `"cs"` remains the default because it provides a low-rank univariate
#'   cubic smooth with shrinkage, which is a useful conservative default when
#'   differentiating an already estimated curve. The alternatives are exposed
#'   to support transparent sensitivity/robustness analyses of this secondary
#'   smoothing step.
#'
#'   `smooth_basis` is validated even when `method = "finite"` for API
#'   consistency, but it has no computational effect in that method and the
#'   output metadata records it as unused.
#' @param elasticity Logical scalar. If `TRUE`, compute local elasticity
#'
#'   \deqn{
#'     E(x) = \frac{dy}{dx}\frac{x}{y}.
#'   }
#'
#'   For `method = "analytical"`, the denominator is the GAM-fitted curve value,
#'   so the derivative and elasticity refer to the same fitted function. For
#'   `method = "finite"`, the denominator is the supplied `y`.
#'
#'   Elasticity is returned as `NA` when its denominator is numerically
#'   indistinguishable from zero.
#' @param critical Logical scalar. If `TRUE`, detect derivative sign changes.
#'   A sign change from positive to negative is classified as a local maximum;
#'   a change from negative to positive is classified as a local minimum.
#'
#'   Because the true zero crossing can lie between evaluated `x` values, the
#'   nearest supplied/evaluated row is flagged and the interpolated crossing is
#'   returned in `critical_x`.
#' @param eps Positive finite numeric scalar smaller than 1. Default is `1e-6`.
#'
#'   It is used as a relative numerical tolerance. For the GAM derivative, the
#'   perturbation is proportional to the observed `x` span and is constrained
#'   to remain inside that span. For elasticity and critical-point detection,
#'   scale-aware tolerances are derived from the corresponding curve/derivative
#'   magnitudes.
#'
#' @return A data frame containing all original input columns, ordered within
#'   group by increasing `x`, plus:
#'
#'   - `sensitivity`: estimated local derivative `dy/dx`;
#'   - `sensitivity_curve_value`: value of the curve whose derivative is being
#'     reported. This equals the GAM-fitted value for `method = "analytical"`
#'     and the supplied `y` for `method = "finite"`;
#'   - `elasticity`, when `elasticity = TRUE`;
#'   - `critical`, `critical_type`, and `critical_x`, when `critical = TRUE`.
#'
#'   `critical_type` is `"local_maximum"` or `"local_minimum"` for detected
#'   derivative sign changes and `NA` otherwise. `critical_x` is a linear
#'   interpolation of the derivative zero crossing between the bracketing
#'   evaluation points.
#'
#'   Output attributes document the derivative method, grouping variables,
#'   numerical contract, requested/effective GAM basis dimensions, and the
#'   auxiliary `smooth_basis` used by the GAM derivative.
#'
#' @details
#' ## GAM-based derivative
#'
#' For `method = "analytical"`, EpiExposure fits the descriptive curve
#'
#' \deqn{
#'   y = s(x) + \epsilon
#' }
#'
#' with REML smoothing-parameter estimation and the auxiliary `mgcv` basis
#' selected through `smooth_basis`.
#'
#' The supported bases are deliberately restricted to four standard
#' one-dimensional choices:
#'
#' - `"cs"`: shrinkage cubic regression spline;
#' - `"cr"`: cubic regression spline;
#' - `"tp"`: thin plate regression spline;
#' - `"ts"`: shrinkage thin plate regression spline.
#'
#' This auxiliary spline is **not the DLNM spline**. It is fitted only after an
#' exposure-response curve has already been produced. For example, a DLNM may
#' have been fitted with a natural spline exposure basis and the resulting
#' curve may still be differentiated here with `smooth_basis = "cs"` or
#' `"tp"`. There is no requirement that the two bases match because they serve
#' different statistical roles.
#'
#' The DLNM basis determines the epidemiological model itself. In contrast,
#' `smooth_basis` controls only how the already evaluated curve is optionally
#' smoothed before numerical differentiation. Users interested in robustness
#' can therefore compare `"cs"`, `"cr"`, `"tp"`, and `"ts"` without changing
#' the previously fitted DLNM.
#'
#' Let \eqn{X_p(x)} denote the GAM linear-predictor matrix and
#' \eqn{\hat\beta} its fitted coefficient vector. The fitted curve is
#'
#' \deqn{
#'   \hat f(x) = X_p(x)\hat\beta.
#' }
#'
#' The derivative mapping is approximated from two nearby prediction matrices:
#'
#' \deqn{
#'   X'_p(x) \approx
#'   \frac{X_p(x_+) - X_p(x_-)}{x_+ - x_-},
#' }
#'
#' and
#'
#' \deqn{
#'   \hat f'(x) = X'_p(x)\hat\beta.
#' }
#'
#' Interior points use a centered perturbation. Boundary points use a
#' one-sided perturbation that remains inside the observed `x` range; the
#' function does not intentionally extrapolate beyond the supplied curve.
#'
#' The external method name `"analytical"` is retained because it was part of
#' the previous EpiExposure API. Methodologically, however, this is a
#' GAM-smoothed numerical derivative, not symbolic differentiation.
#'
#' ## Why `"cs"` remains the default
#'
#' Derivatives can amplify small-scale wiggles in an estimated curve. A
#' shrinkage cubic regression spline is therefore retained as the default
#' auxiliary smoother because it offers a simple one-dimensional cubic basis
#' while allowing stronger penalization of weak structure. This is a
#' post-processing default, not a statement that `"cs"` is universally optimal.
#'
#' `smooth_basis = "cr"` is useful when the analyst wants the corresponding
#' cubic regression spline without the extra shrinkage modification.
#' `"tp"` provides the general thin plate regression spline commonly used by
#' `mgcv`, while `"ts"` adds shrinkage to that family.
#'
#' Because the estimated derivative can depend on the secondary smoothing
#' choice, reporting or checking sensitivity across these bases can be useful
#' when derivative-based scientific conclusions depend on fine features of the
#' curve.
#'
#' ## Direct finite differences
#'
#' With `method = "finite"`, no smoother is fitted. For three or more points,
#' the derivative at an interior point is obtained from the quadratic
#' interpolant through the previous, current, and next observations. This gives
#' the standard centered three-point derivative for equally spaced data and its
#' corresponding unequal-spacing form for irregular `x`.
#'
#' Endpoints use the corresponding one-sided three-point formula. With exactly
#' two points, the secant slope is the only estimable derivative and is assigned
#' to both positions.
#'
#' `x` values must be strictly unique after grouping. The function does not
#' average duplicate positions, discard rows, or replace a zero spacing by
#' `eps`.
#'
#' ## Elasticity
#'
#' Elasticity is dimensionless:
#'
#' \deqn{
#'   E(x) = f'(x)x/f(x).
#' }
#'
#' For the GAM method, both \eqn{f'(x)} and \eqn{f(x)} come from the same fitted
#' smooth. For the direct finite method, \eqn{f(x)} is the supplied `y`.
#'
#' Elasticity can be numerically unstable near zero. EpiExposure therefore
#' returns `NA` when the relevant curve value is within a scale-aware tolerance
#' of zero. No arbitrary constant is added to the denominator.
#'
#' ## Critical points
#'
#' `critical = TRUE` searches for changes in the sign of the estimated
#' derivative. Near-zero derivative values are treated as a bridge between the
#' nearest non-zero derivative signs. Only genuine positive-to-negative or
#' negative-to-positive transitions are marked.
#'
#' Therefore:
#'
#' - positive -> negative: local maximum;
#' - negative -> positive: local minimum.
#'
#' Flat regions without a sign reversal are not automatically labelled as
#' extrema. Boundary extrema are also not inferred from a one-sided derivative
#' alone.
#'
#' ## Relationship to the EpiExposure exact-history contract
#'
#' `df_sensitivity()` does not consume raw fitted exposure histories. It works
#' on an already evaluated curve, whose number of rows may legitimately differ
#' from the original history length. Consequently, it does **not** require
#'
#' \deqn{
#'   nrow(data) = max\_lag + 1.
#' }
#'
#' If current EpiExposure temporal metadata are present on `data`, however, the
#' function validates and propagates them. Specifically,
#'
#' \deqn{
#'   history\_length = max\_lag + 1
#' }
#'
#' and the stored history contract must be
#' `"all_fitted_exposures_same_exact_max_lag_plus_one"`.
#'
#' This validates the provenance of an EpiExposure-derived curve without
#' incorrectly treating curve-evaluation rows as original exposure-history
#' observations.
#'
#' ## Uncertainty
#'
#' `df_sensitivity()` does not generate coefficient/posterior draws. If the
#' supplied data contain draw-specific curves, include the draw identifier in
#' `scenario_var` so each draw is differentiated independently. The function
#' does not silently pool repeated `x` positions across draws.
#'
#' @export
df_sensitivity <- function(
    data,
    x = "value",
    y = "prediction",
    scenario_var = "scenario",
    k = 10,
    method = c("analytical", "finite"),
    smooth_basis = c("cs", "cr", "tp", "ts"),
    elasticity = TRUE,
    critical = TRUE,
    eps = 1e-6
) {

  # ==========================================================================
  # SMALL VALIDATORS
  # ==========================================================================

  valid_name <- function(z) {
    is.character(z) &&
      length(z) == 1L &&
      !is.na(z) &&
      nzchar(z)
  }

  valid_flag <- function(z) {
    is.logical(z) &&
      length(z) == 1L &&
      !is.na(z)
  }

  valid_integer_scalar <- function(z) {
    is.numeric(z) &&
      length(z) == 1L &&
      !is.na(z) &&
      is.finite(z) &&
      z == as.integer(z)
  }

  group_description <- function(d, vars) {
    if (is.null(vars)) {
      return("the supplied curve")
    }

    values <- vapply(
      vars,
      function(v) {
        paste0(v, "=", as.character(d[[v]][1L]))
      },
      character(1)
    )

    paste(values, collapse = ", ")
  }

  scale_tolerance <- function(z, relative_eps) {
    z <- as.numeric(z)
    finite_z <- z[is.finite(z)]

    if (!length(finite_z)) {
      return(NA_real_)
    }

    magnitude <- max(abs(finite_z))

    if (magnitude == 0) {
      return(relative_eps)
    }

    relative_eps * magnitude
  }

  # ==========================================================================
  # ARGUMENT MATCHING AND GENERAL VALIDATION
  # ==========================================================================

  method <- match.arg(method)
  smooth_basis <- match.arg(smooth_basis)

  if (!is.data.frame(data) ||
      !nrow(data)) {
    stop(
      "`data` must be a non-empty data.frame.",
      call. = FALSE
    )
  }

  if (!valid_name(x)) {
    stop(
      "`x` must be one non-empty character column name.",
      call. = FALSE
    )
  }

  if (!valid_name(y)) {
    stop(
      "`y` must be one non-empty character column name.",
      call. = FALSE
    )
  }

  if (identical(x, y)) {
    stop(
      "`x` and `y` must identify different columns.",
      call. = FALSE
    )
  }

  missing_xy <- setdiff(
    c(x, y),
    names(data)
  )

  if (length(missing_xy)) {
    stop(
      "`data` is missing required column(s): ",
      paste(missing_xy, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[x]]) ||
      anyNA(data[[x]]) ||
      any(!is.finite(data[[x]]))) {
    stop(
      "`", x, "` must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }

  if (!is.numeric(data[[y]]) ||
      anyNA(data[[y]]) ||
      any(!is.finite(data[[y]]))) {
    stop(
      "`", y, "` must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }

  if (!is.null(scenario_var)) {
    if (!is.character(scenario_var) ||
        !length(scenario_var) ||
        anyNA(scenario_var) ||
        any(!nzchar(scenario_var)) ||
        anyDuplicated(scenario_var)) {
      stop(
        "`scenario_var` must be NULL or contain unique, non-empty grouping ",
        "column names.",
        call. = FALSE
      )
    }

    missing_groups <- setdiff(
      scenario_var,
      names(data)
    )

    if (length(missing_groups)) {
      stop(
        "Grouping column(s) not found in `data`: ",
        paste(missing_groups, collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (any(scenario_var %in% c(x, y))) {
      stop(
        "`scenario_var` cannot include the `x` or `y` column.",
        call. = FALSE
      )
    }

    for (v in scenario_var) {
      if (anyNA(data[[v]])) {
        stop(
          "Grouping column '",
          v,
          "' cannot contain missing values.",
          call. = FALSE
        )
      }

      if (is.character(data[[v]]) &&
          any(!nzchar(data[[v]]))) {
        stop(
          "Grouping column '",
          v,
          "' cannot contain empty character labels.",
          call. = FALSE
        )
      }
    }
  }

  if (!valid_integer_scalar(k) ||
      k < 1L) {
    stop(
      "`k` must be one positive integer.",
      call. = FALSE
    )
  }
  k <- as.integer(k)

  if (identical(method, "analytical") &&
      k < 3L) {
    stop(
      "`method = 'analytical'` requires `k >= 3`.",
      call. = FALSE
    )
  }

  if (!valid_flag(elasticity)) {
    stop(
      "`elasticity` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!valid_flag(critical)) {
    stop(
      "`critical` must be TRUE or FALSE.",
      call. = FALSE
    )
  }

  if (!is.numeric(eps) ||
      length(eps) != 1L ||
      is.na(eps) ||
      !is.finite(eps) ||
      eps <= 0 ||
      eps >= 1) {
    stop(
      "`eps` must be one finite numeric value satisfying 0 < eps < 1.",
      call. = FALSE
    )
  }
  eps <- as.numeric(eps)

  if (identical(method, "analytical") &&
      !requireNamespace("mgcv", quietly = TRUE)) {
    stop(
      "Package 'mgcv' is required for `method = 'analytical'`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # OPTIONAL EPIEXPOSURE TEMPORAL-METADATA VALIDATION
  # ==========================================================================

  max_lag_metadata <- attr(
    data,
    "epiexposure_max_lag",
    exact = TRUE
  )

  history_length_metadata <- attr(
    data,
    "epiexposure_history_length",
    exact = TRUE
  )

  history_contract_metadata <- attr(
    data,
    "epiexposure_history_contract",
    exact = TRUE
  )

  temporal_present <- c(
    !is.null(max_lag_metadata),
    !is.null(history_length_metadata),
    !is.null(history_contract_metadata)
  )

  has_temporal_metadata <- all(temporal_present)

  if (any(temporal_present) &&
      !has_temporal_metadata) {
    stop(
      "`data` contains incomplete EpiExposure temporal metadata. When any of ",
      "`epiexposure_max_lag`, `epiexposure_history_length`, or ",
      "`epiexposure_history_contract` is present, all three must be present.",
      call. = FALSE
    )
  }

  if (has_temporal_metadata) {
    if (!valid_integer_scalar(max_lag_metadata) ||
        max_lag_metadata < 0L) {
      stop(
        "`data` contains invalid `epiexposure_max_lag` metadata.",
        call. = FALSE
      )
    }
    max_lag_metadata <- as.integer(max_lag_metadata)

    if (!valid_integer_scalar(history_length_metadata) ||
        history_length_metadata < 1L) {
      stop(
        "`data` contains invalid `epiexposure_history_length` metadata.",
        call. = FALSE
      )
    }
    history_length_metadata <- as.integer(
      history_length_metadata
    )

    if (!identical(
      history_length_metadata,
      max_lag_metadata + 1L
    )) {
      stop(
        "EpiExposure temporal metadata are inconsistent: ",
        "`history_length` must equal `max_lag + 1`. Expected ",
        max_lag_metadata + 1L,
        " but found ",
        history_length_metadata,
        ".",
        call. = FALSE
      )
    }

    if (!identical(
      history_contract_metadata,
      "all_fitted_exposures_same_exact_max_lag_plus_one"
    )) {
      stop(
        "`data` does not satisfy the current EpiExposure exact-history ",
        "contract.",
        call. = FALSE
      )
    }
  }

  # ==========================================================================
  # GROUP INDEX
  # ==========================================================================

  if (is.null(scenario_var)) {
    group_rows <- list(
      all = seq_len(nrow(data))
    )

  } else {
    interaction_args <- c(
      unname(
        data[
          ,
          scenario_var,
          drop = FALSE
        ]
      ),
      list(
        drop = TRUE,
        lex.order = TRUE,
        sep = "\r"
      )
    )

    group_factor <- do.call(
      interaction,
      interaction_args
    )

    group_rows <- split(
      seq_len(nrow(data)),
      group_factor,
      drop = TRUE
    )

    first_rows <- vapply(
      group_rows,
      function(idx) idx[1L],
      integer(1)
    )

    group_rows <- group_rows[
      order(first_rows)
    ]
  }

  if (!length(group_rows)) {
    stop(
      "No curve groups were available after applying `scenario_var`.",
      call. = FALSE
    )
  }

  # ==========================================================================
  # PER-GROUP CURVE VALIDATION
  # ==========================================================================

  validate_curve <- function(d) {
    group_label <- group_description(
      d,
      scenario_var
    )

    n_current <- nrow(d)

    if (n_current < 2L) {
      stop(
        "At least two curve points are required within ",
        group_label,
        "; received ",
        n_current,
        ".",
        call. = FALSE
      )
    }

    x_values <- as.numeric(
      d[[x]]
    )

    if (anyDuplicated(x_values)) {
      duplicate_values <- unique(
        x_values[
          duplicated(x_values) |
            duplicated(
              x_values,
              fromLast = TRUE
            )
        ]
      )

      stop(
        "Duplicate `",
        x,
        "` values were detected within ",
        group_label,
        ": ",
        paste(
          utils::head(
            duplicate_values,
            10L
          ),
          collapse = ", "
        ),
        if (length(duplicate_values) > 10L) {
          ", ..."
        } else {
          ""
        },
        ". `df_sensitivity()` does not aggregate duplicate curve positions ",
        "silently.",
        call. = FALSE
      )
    }

    d <- d[
      order(x_values),
      ,
      drop = FALSE
    ]
    rownames(d) <- NULL

    x_values <- as.numeric(
      d[[x]]
    )

    dx <- diff(
      x_values
    )

    if (any(dx <= 0) ||
        any(!is.finite(dx))) {
      stop(
        "`",
        x,
        "` must increase strictly after ordering within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    x_span <- max(x_values) -
      min(x_values)

    if (!is.finite(x_span) ||
        x_span <= 0) {
      stop(
        "`",
        x,
        "` must span at least two distinct finite values within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    machine_gap_tolerance <-
      10 *
      .Machine$double.eps *
      max(
        abs(x_values),
        x_span,
        .Machine$double.xmin
      )

    if (any(dx <= machine_gap_tolerance)) {
      stop(
        "At least two `",
        x,
        "` positions are numerically indistinguishable within ",
        group_label,
        ". Increase the spacing or rescale the predictor before computing ",
        "derivatives.",
        call. = FALSE
      )
    }

    if (identical(method, "analytical") &&
        n_current < 3L) {
      stop(
        "`method = 'analytical'` requires at least three distinct `",
        x,
        "` values within ",
        group_label,
        "; received ",
        n_current,
        ". Use `method = 'finite'` for a two-point curve.",
        call. = FALSE
      )
    }

    d
  }

  # ==========================================================================
  # GAM-BASED DERIVATIVE
  # ==========================================================================

  compute_derivative_gam <- function(d) {
    group_label <- group_description(
      d,
      scenario_var
    )

    x_values <- as.numeric(
      d[[x]]
    )

    y_values <- as.numeric(
      d[[y]]
    )

    n_unique_x <- length(
      unique(x_values)
    )

    k_use <- min(
      k,
      n_unique_x
    )

    if (k_use < 3L) {
      stop(
        "The effective GAM basis dimension is smaller than 3 within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    fit_data <- data.frame(
      .epix_x = x_values,
      .epix_y = y_values,
      check.names = FALSE
    )

    # A constant supplied curve has an exact zero derivative. Avoid asking the
    # smoother to identify curvature that is mathematically absent.
    if (length(unique(y_values)) == 1L) {
      d$sensitivity <- rep(
        0,
        nrow(d)
      )
      d$sensitivity_curve_value <- y_values

      attr(
        d,
        ".epix_effective_k"
      ) <- k_use

      return(d)
    }

    formula_environment <- new.env(
      parent = environment()
    )
    formula_environment$s <- mgcv::s
    formula_environment$k_use <- k_use
    formula_environment$smooth_basis <- smooth_basis

    gam_formula <- stats::as.formula(
      ".epix_y ~ s(.epix_x, k = k_use, bs = smooth_basis)",
      env = formula_environment
    )

    fit <- tryCatch(
      mgcv::gam(
        formula = gam_formula,
        data = fit_data,
        method = "REML"
      ),
      error = function(e) {
        stop(
          "Could not fit the GAM sensitivity curve with `smooth_basis = '",
          smooth_basis,
          "'` within ",
          group_label,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    coefficients <- stats::coef(
      fit
    )

    if (!length(coefficients) ||
        anyNA(coefficients) ||
        any(!is.finite(coefficients))) {
      stop(
        "The GAM sensitivity fit returned non-finite coefficients within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    fitted_curve <- tryCatch(
      as.numeric(
        stats::predict(
          fit,
          newdata = fit_data,
          type = "response"
        )
      ),
      error = function(e) {
        stop(
          "Could not evaluate the fitted GAM sensitivity curve within ",
          group_label,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (length(fitted_curve) != nrow(d) ||
        anyNA(fitted_curve) ||
        any(!is.finite(fitted_curve))) {
      stop(
        "The GAM sensitivity curve produced invalid fitted values within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    x_min <- min(
      x_values
    )
    x_max <- max(
      x_values
    )
    x_span <- x_max - x_min

    perturbation <-
      max(
        eps * x_span,
        10 *
          .Machine$double.eps *
          max(
            abs(x_values),
            x_span,
            .Machine$double.xmin
          )
      )

    x_lower <- pmax(
      x_min,
      x_values - perturbation
    )

    x_upper <- pmin(
      x_max,
      x_values + perturbation
    )

    denominator <- x_upper -
      x_lower

    if (any(denominator <= 0) ||
        any(!is.finite(denominator))) {
      stop(
        "Could not construct a stable within-range perturbation for the GAM ",
        "derivative within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    lower_data <- data.frame(
      .epix_x = x_lower,
      .epix_y = y_values,
      check.names = FALSE
    )

    upper_data <- data.frame(
      .epix_x = x_upper,
      .epix_y = y_values,
      check.names = FALSE
    )

    X_lower <- tryCatch(
      stats::predict(
        fit,
        newdata = lower_data,
        type = "lpmatrix"
      ),
      error = function(e) {
        stop(
          "Could not evaluate the lower GAM `lpmatrix` within ",
          group_label,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    X_upper <- tryCatch(
      stats::predict(
        fit,
        newdata = upper_data,
        type = "lpmatrix"
      ),
      error = function(e) {
        stop(
          "Could not evaluate the upper GAM `lpmatrix` within ",
          group_label,
          ": ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )

    if (!is.matrix(X_lower) ||
        !is.matrix(X_upper) ||
        !identical(
          dim(X_lower),
          dim(X_upper)
        ) ||
        nrow(X_lower) != nrow(d) ||
        ncol(X_lower) != length(coefficients)) {
      stop(
        "The GAM prediction matrices are inconsistent with the fitted ",
        "coefficient vector within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    derivative_matrix <-
      sweep(
        X_upper - X_lower,
        MARGIN = 1L,
        STATS = denominator,
        FUN = "/"
      )

    sensitivity <- as.numeric(
      derivative_matrix %*%
        coefficients
    )

    if (length(sensitivity) != nrow(d) ||
        anyNA(sensitivity) ||
        any(!is.finite(sensitivity))) {
      stop(
        "The GAM derivative calculation produced non-finite sensitivity values ",
        "within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    d$sensitivity <- sensitivity
    d$sensitivity_curve_value <- fitted_curve

    attr(
      d,
      ".epix_effective_k"
    ) <- k_use

    d
  }

  # ==========================================================================
  # DIRECT FINITE-DIFFERENCE DERIVATIVE
  # ==========================================================================

  compute_derivative_finite <- function(d) {
    group_label <- group_description(
      d,
      scenario_var
    )

    xv <- as.numeric(
      d[[x]]
    )

    yv <- as.numeric(
      d[[y]]
    )

    n_points <- length(
      xv
    )

    sensitivity <- numeric(
      n_points
    )

    if (n_points == 2L) {
      slope <- (
        yv[2L] - yv[1L]
      ) / (
        xv[2L] - xv[1L]
      )

      sensitivity[] <- slope

    } else {
      # Left boundary: derivative of the quadratic interpolant through points
      # 1, 2, and 3, evaluated at x1.
      h1 <- xv[2L] - xv[1L]
      h2 <- xv[3L] - xv[2L]

      sensitivity[1L] <-
        -(2 * h1 + h2) /
        (h1 * (h1 + h2)) *
        yv[1L] +
        (h1 + h2) /
        (h1 * h2) *
        yv[2L] -
        h1 /
        (h2 * (h1 + h2)) *
        yv[3L]

      # Interior points: derivative of the local quadratic interpolant on an
      # unequal grid.
      if (n_points > 2L) {
        for (i in 2L:(n_points - 1L)) {
          h_left <- xv[i] -
            xv[i - 1L]

          h_right <- xv[i + 1L] -
            xv[i]

          sensitivity[i] <-
            -h_right /
            (
              h_left *
                (h_left + h_right)
            ) *
            yv[i - 1L] +
            (
              h_right - h_left
            ) /
            (
              h_left *
                h_right
            ) *
            yv[i] +
            h_left /
            (
              h_right *
                (h_left + h_right)
            ) *
            yv[i + 1L]
        }
      }

      # Right boundary: mirrored one-sided quadratic derivative.
      h_left2 <- xv[n_points - 1L] -
        xv[n_points - 2L]

      h_left1 <- xv[n_points] -
        xv[n_points - 1L]

      sensitivity[n_points] <-
        h_left1 /
        (
          h_left2 *
            (h_left2 + h_left1)
        ) *
        yv[n_points - 2L] -
        (
          h_left2 + h_left1
        ) /
        (
          h_left2 *
            h_left1
        ) *
        yv[n_points - 1L] +
        (
          2 * h_left1 + h_left2
        ) /
        (
          h_left1 *
            (h_left2 + h_left1)
        ) *
        yv[n_points]
    }

    if (anyNA(sensitivity) ||
        any(!is.finite(sensitivity))) {
      stop(
        "The direct finite-difference calculation produced non-finite ",
        "sensitivity values within ",
        group_label,
        ".",
        call. = FALSE
      )
    }

    d$sensitivity <- sensitivity
    d$sensitivity_curve_value <- yv

    attr(
      d,
      ".epix_effective_k"
    ) <- NA_integer_

    d
  }

  # ==========================================================================
  # ELASTICITY
  # ==========================================================================

  add_elasticity <- function(d) {
    if (!elasticity) {
      return(d)
    }

    denominator <- as.numeric(
      d$sensitivity_curve_value
    )

    denominator_tolerance <- scale_tolerance(
      denominator,
      eps
    )

    defined <- abs(denominator) >
      denominator_tolerance

    elasticity_values <- rep(
      NA_real_,
      nrow(d)
    )

    elasticity_values[defined] <-
      d$sensitivity[defined] *
      (
        as.numeric(
          d[[x]][defined]
        ) /
          denominator[defined]
      )

    if (any(
      !is.na(elasticity_values) &
      !is.finite(elasticity_values)
    )) {
      stop(
        "Elasticity calculation produced a non-finite value.",
        call. = FALSE
      )
    }

    d$elasticity <- elasticity_values

    d
  }

  # ==========================================================================
  # CRITICAL-POINT DETECTION
  # ==========================================================================

  add_critical_points <- function(d) {
    if (!critical) {
      return(d)
    }

    xv <- as.numeric(
      d[[x]]
    )

    s <- as.numeric(
      d$sensitivity
    )

    derivative_tolerance <- scale_tolerance(
      s,
      eps
    )

    sign_class <- rep(
      0L,
      length(s)
    )

    sign_class[
      s > derivative_tolerance
    ] <- 1L

    sign_class[
      s < -derivative_tolerance
    ] <- -1L

    nonzero_index <- which(
      sign_class != 0L
    )

    critical_flag <- rep(
      FALSE,
      length(s)
    )

    critical_type <- rep(
      NA_character_,
      length(s)
    )

    critical_x <- rep(
      NA_real_,
      length(s)
    )

    if (length(nonzero_index) >= 2L) {
      for (j in seq_len(
        length(nonzero_index) - 1L
      )) {
        i_left <- nonzero_index[j]
        i_right <- nonzero_index[j + 1L]

        left_sign <- sign_class[
          i_left
        ]
        right_sign <- sign_class[
          i_right
        ]

        if (left_sign ==
            right_sign) {
          next
        }

        s_left <- s[
          i_left
        ]
        s_right <- s[
          i_right
        ]

        x_left <- xv[
          i_left
        ]
        x_right <- xv[
          i_right
        ]

        root <- if (
          is.finite(s_left) &&
          is.finite(s_right) &&
          !identical(
            s_left,
            s_right
          )
        ) {
          x_left -
            s_left *
            (
              x_right - x_left
            ) /
            (
              s_right - s_left
            )
        } else {
          mean(
            c(
              x_left,
              x_right
            )
          )
        }

        if (!is.finite(root) ||
            root < x_left ||
            root > x_right) {
          root <- mean(
            c(
              x_left,
              x_right
            )
          )
        }

        candidate_rows <- seq.int(
          i_left,
          i_right
        )

        nearest <- candidate_rows[
          which.min(
            abs(
              xv[candidate_rows] -
                root
            )
          )
        ]

        critical_flag[
          nearest
        ] <- TRUE

        critical_type[
          nearest
        ] <- if (
          left_sign > 0L &&
          right_sign < 0L
        ) {
          "local_maximum"
        } else {
          "local_minimum"
        }

        critical_x[
          nearest
        ] <- root
      }
    }

    d$critical <- critical_flag
    d$critical_type <- critical_type
    d$critical_x <- critical_x

    d
  }

  # ==========================================================================
  # PROCESS EACH CURVE
  # ==========================================================================

  result_list <- vector(
    "list",
    length(group_rows)
  )

  effective_k <- integer(
    length(group_rows)
  )

  for (i in seq_along(group_rows)) {
    d <- data[
      group_rows[[i]],
      ,
      drop = FALSE
    ]

    d <- validate_curve(
      d
    )

    if (identical(
      method,
      "analytical"
    )) {
      d <- compute_derivative_gam(
        d
      )

    } else {
      d <- compute_derivative_finite(
        d
      )
    }

    effective_k[i] <- attr(
      d,
      ".epix_effective_k",
      exact = TRUE
    )

    attr(
      d,
      ".epix_effective_k"
    ) <- NULL

    d <- add_elasticity(
      d
    )

    d <- add_critical_points(
      d
    )

    result_list[[i]] <- d
  }

  out <- do.call(
    rbind,
    result_list
  )

  rownames(out) <- NULL

  # ==========================================================================
  # FINAL ORDER
  # ==========================================================================

  if (is.null(scenario_var)) {
    out <- out[
      order(
        out[[x]]
      ),
      ,
      drop = FALSE
    ]

  } else {
    order_args <- c(
      lapply(
        scenario_var,
        function(v) out[[v]]
      ),
      list(
        out[[x]]
      )
    )

    out <- out[
      do.call(
        order,
        order_args
      ),
      ,
      drop = FALSE
    ]
  }

  rownames(out) <- NULL

  # ==========================================================================
  # OUTPUT METADATA
  # ==========================================================================

  attr(
    out,
    "epiexposure_sensitivity_method"
  ) <- method

  attr(
    out,
    "epiexposure_sensitivity_x"
  ) <- x

  attr(
    out,
    "epiexposure_sensitivity_y"
  ) <- y

  attr(
    out,
    "epiexposure_sensitivity_group_vars"
  ) <- scenario_var

  attr(
    out,
    "epiexposure_sensitivity_smooth_basis"
  ) <- if (
    identical(
      method,
      "analytical"
    )
  ) {
    smooth_basis
  } else {
    NULL
  }

  attr(
    out,
    "epiexposure_sensitivity_smooth_basis_requested"
  ) <- smooth_basis

  attr(
    out,
    "epiexposure_sensitivity_k_requested"
  ) <- if (
    identical(
      method,
      "analytical"
    )
  ) {
    k
  } else {
    NULL
  }

  attr(
    out,
    "epiexposure_sensitivity_k_effective"
  ) <- if (
    identical(
      method,
      "analytical"
    )
  ) {
    effective_k
  } else {
    NULL
  }

  attr(
    out,
    "epiexposure_sensitivity_eps"
  ) <- eps

  attr(
    out,
    "epiexposure_sensitivity_elasticity"
  ) <- elasticity

  attr(
    out,
    "epiexposure_sensitivity_critical"
  ) <- critical

  attr(
    out,
    "epiexposure_sensitivity_derivative_contract"
  ) <- if (
    identical(
      method,
      "analytical"
    )
  ) {
    paste0(
      "gam_reml_",
      smooth_basis,
      "_lpmatrix_within_range_finite_difference"
    )
  } else {
    "direct_three_point_finite_difference_supporting_irregular_x"
  }

  attr(
    out,
    "epiexposure_sensitivity_smooth_basis_contract"
  ) <- if (
    identical(
      method,
      "analytical"
    )
  ) {
    "auxiliary_postprocessing_gam_basis_not_dlnm_basis"
  } else {
    "unused_for_direct_finite_difference"
  }

  attr(
    out,
    "epiexposure_sensitivity_elasticity_contract"
  ) <- if (
    identical(
      method,
      "analytical"
    )
  ) {
    "gam_derivative_times_x_over_gam_fitted_curve"
  } else {
    "finite_derivative_times_x_over_supplied_curve"
  }

  attr(
    out,
    "epiexposure_sensitivity_critical_contract"
  ) <- if (critical) {
    "derivative_sign_reversal_with_interpolated_zero_crossing"
  } else {
    NULL
  }

  attr(
    out,
    "epiexposure_sensitivity_contract"
  ) <- "postprocessed_curve_derivative_not_raw_exposure_history"

  if (has_temporal_metadata) {
    attr(
      out,
      "epiexposure_max_lag"
    ) <- max_lag_metadata

    attr(
      out,
      "epiexposure_history_length"
    ) <- history_length_metadata

    attr(
      out,
      "epiexposure_history_contract"
    ) <- history_contract_metadata
  }

  out
}
