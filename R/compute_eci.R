#' Compute Exposure Cumulative Impact (ECI - the exposure profile with the fitted model coefficients)#' Compute Exposure Cumulative Impact (ECI_weighted`) for one or more exposure variables.
#'
#' It supports both deterministic estimation and uncertainty propagation.
#' When `uncertainty = TRUE`, `ECI_weighted` is recomputed across simulated or
#' posterior draws of the model coefficients.
#'
#' If `output = "summary"`, the central estimate is computed as the median of
#' the simulated ECI values, while interval limits are obtained from empirical
#' quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** when `uncertainty = TRUE` and `output = "summary"`,
#' `ECI_weighted` represents the *central estimate*, computed as the median
#' of the simulated distribution on the requested `scale`.
#'
#' @param profile Numeric vector or list of numeric vectors representing one or
#'   more lag profiles. If `fit = NULL`, `profile` must be a single numeric
#'   vector. If `fit` is provided and `var` contains multiple variables,
#'   `profile` can be a named list of profiles (one per variable), or a list
#'   that will be matched to `var` by position.
#'
#' @param epi_data Optional data frame containing observed exposure histories.
#'   Used only when `fit` is provided and the user wants the function to extract
#'   the profile(s) from a grouped dataset instead of supplying `profile`
#'   directly.
#'
#'   **IMPORTANT:** data in `epi_data` must already be ordered in the correct
#'   chronological sequence within each group (oldest → most recent), with no
#'   misplaced days or misordered observations. This function does **not**
#'   reorder the observed data internally.
#'
#' @param group Optional column name used to identify groups in `epi_data`
#'   (e.g., `"epi_id"`).
#'
#' @param group_level Optional value or vector of values inside `group` used to
#'   select the profile(s) from `epi_data`. This represents the sub-level(s) of
#'   the grouping variable (e.g., one or more epidemic IDs). It can be numeric,
#'   character, or factor.
#'
#'   If `group_level = NULL`, the function computes ECI estimates for **all**
#'   observed levels of `group`.
#'
#' @param fit Fitted model from `fit_epidlnm()`. Required to compute
#'   `ECI_weighted`. If `NULL`, only `ECI_raw` is returned.
#'
#' @param var Optional character scalar or vector indicating which fitted
#'   exposure variable(s) should be used. Required when `fit` contains more
#'   than one exposure variable. Names in `var` must exactly match the variable
#'   names stored in the fitted model.
#'
#' @param scale Character. Scale for `ECI_weighted`:
#'   - `"link"`: linear predictor scale (default)
#'   - `"response"`: inverse-link transformed scale
#'   - `"percent"`: relative change scale, computed as `(exp(eta) - 1) * 100`
#'
#' @param reverse Logical. If `TRUE`, reverses the input profile(s) before analysis.
#'   Use this when the supplied profile(s) are ordered from oldest to most recent
#'   (chronological order). If `FALSE`, the function assumes the first value
#'   corresponds to lag 0 (most recent) and the last value to lag max (oldest).
#'
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty.
#'
#' @param output Character. `"summary"` or `"samples"`.
#'
#' @param n_samples Integer. Number of samples used for uncertainty quantification.
#'
#' @return A data.frame.
#'
#' - If `fit = NULL`, returns:
#'   - `ECI_raw`
#'   - `ECI_weighted = NA`
#'
#' - If `uncertainty = FALSE`, returns:
#'   - `ECI_raw`
#'   - `ECI_weighted`
#'   - and, when multiple variables are requested, a `var` column
#'   - and, when `epi_data` is used, a grouping column named as `group`
#'
#' - If `uncertainty = TRUE` and `output = "summary"`, returns:
#'   - `ECI_raw`
#'   - `ECI_weighted`
#'   - `sd`
#'   - `lower`
#'   - `upper`
#'   - and, when multiple variables are requested, a `var` column
#'   - and, when `epi_data` is used, a grouping column named as `group`
#'
#' - If `uncertainty = TRUE` and `output = "samples"`, returns:
#'   - `sample`
#'   - `ECI_raw`
#'   - `ECI_weighted`
#'   - and, when multiple variables are requested, a `var` column
#'   - and, when `epi_data` is used, a grouping column named as `group`
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use simulation from the asymptotic coefficient distribution
#'
#' The use of the median as the central estimate improves robustness under
#' asymmetric or non-normal simulated ECI distributions.
#'
#' For `scale = "response"`, the inverse link is obtained from the fitted model
#' metadata when available.
#'
#' For `scale = "percent"`, the transformation is defined as:
#' \deqn{(exp(\eta) - 1) * 100}
#' which is most directly interpretable for log-linked or relative-effect contexts.
#'
#' @export
compute_eci <- function(
    profile = NULL,
    epi_data = NULL,
    group = NULL,
    group_level = NULL,
    fit = NULL,
    var = NULL,
    scale = c("link", "response", "percent"),
    reverse = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  scale <- match.arg(scale)
  output <- match.arg(output)

  # -----------------------
  # Validation
  # -----------------------
  if (!is.logical(reverse) || length(reverse) != 1L) {
    stop("`reverse` must be TRUE or FALSE.")
  }

  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }

  if (!is.numeric(n_samples) || length(n_samples) != 1L || !is.finite(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }

  n_samples <- as.integer(n_samples)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(0)
    stats::sd(x)
  }

  # ----------------------------------------------------------
  # helper: robust matching of brms draw names
  # ----------------------------------------------------------
  match_brms_draw_names <- function(cb_names_ref, draw_colnames) {

    bnames <- paste0("b_", cb_names_ref)

    if (all(bnames %in% draw_colnames)) {
      return(bnames)
    }

    raw_match <- cb_names_ref[cb_names_ref %in% draw_colnames]
    if (length(raw_match) == length(cb_names_ref)) {
      return(raw_match)
    }

    stop("Could not match brms posterior draw names to crossbasis columns.")
  }

  # ----------------------------------------------------------
  # Helper: order coefficient names robustly
  # ----------------------------------------------------------
  sort_cb_names <- function(x) {
    if (length(x) == 0) return(x)

    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    idx[is.na(idx)] <- seq_along(x)

    x[order(idx)]
  }

  # ----------------------------------------------------------
  # Helper: inverse link
  # ----------------------------------------------------------
  get_linkinv <- function(family_fit) {

    if (is.character(family_fit)) {
      if (family_fit %in% c("beta", "binomial")) return(plogis)
      if (family_fit %in% c("poisson", "gamma")) return(exp)
      if (family_fit %in% c("gaussian")) return(identity)
    }

    if (inherits(family_fit, "family") && !is.null(family_fit$linkinv)) {
      return(family_fit$linkinv)
    }

    identity
  }

  transform_scale <- function(x, scale, linkinv) {
    if (scale == "link") return(x)
    if (scale == "response") return(linkinv(x))
    if (scale == "percent") return((exp(x) - 1) * 100)
    stop("Unsupported scale.")
  }

  # ----------------------------------------------------------
  # Helper: deterministic coefficients + vcov
  # ----------------------------------------------------------
  extract_coef_vcov <- function(model) {

    if (inherits(model, "glmmTMB")) {
      return(list(
        beta = glmmTMB::fixef(model)$cond,
        vcov = as.matrix(stats::vcov(model)$cond)
      ))
    }

    if (inherits(model, "merMod")) {
      return(list(
        beta = lme4::fixef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "lme")) {
      return(list(
        beta = nlme::fixef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "gls")) {
      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "gam")) {
      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "HLfit")) {
      V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
      if (is.null(V)) stop("Could not extract vcov from spaMM model.")
      return(list(
        beta = spaMM::fixef(model),
        vcov = V
      ))
    }

    if (inherits(model, "brmsfit")) {
      b <- brms::fixef(model)
      beta <- b[, "Estimate"]
      names(beta) <- rownames(b)
      return(list(
        beta = beta,
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      V <- diag(model$summary.fixed$sd^2)
      rownames(V) <- names(beta)
      colnames(V) <- names(beta)
      return(list(beta = beta, vcov = V))
    }

    if (inherits(model, "bdlnm")) {
      beta <- model$coefficients.summary[, "mean"]
      V <- stats::cov(t(model$coefficients))
      return(list(beta = beta, vcov = V))
    }

    beta <- stats::coef(model)
    if (!is.numeric(beta)) {
      stop("Could not extract numeric coefficients from fitted model.")
    }

    return(list(
      beta = beta,
      vcov = as.matrix(stats::vcov(model))
    ))
  }

  # ----------------------------------------------------------
  # Helper: get cb names
  # ----------------------------------------------------------
  get_cb_names <- function(coef_names, var_sel, cb_cols_fit = NULL) {

    cb_names <- grep(paste0("^cb_", var_sel, "_"), coef_names, value = TRUE)

    if (length(cb_names) == 0 && !is.null(cb_cols_fit)) {
      cb_names <- intersect(
        grep(paste0("^cb_", var_sel, "_"), cb_cols_fit, value = TRUE),
        coef_names
      )
    }

    cb_names <- sort_cb_names(cb_names)

    cb_names
  }

  # ----------------------------------------------------------
  # Helper: posterior / simulated draws by engine
  # ----------------------------------------------------------
  extract_beta_draws <- function(model, cb_names_ref, n_samples, cv) {

    if (inherits(model, "brmsfit")) {

      draws <- as.matrix(brms::as_draws_matrix(model))
      draw_names <- match_brms_draw_names(cb_names_ref, colnames(draws))

      draws <- draws[, draw_names, drop = FALSE]

      if (nrow(draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(nrow(draws)), n_samples)
        draws <- draws[keep, , drop = FALSE]
      }

      colnames(draws) <- cb_names_ref
      return(draws)
    }

    if (inherits(model, "inla")) {

      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("Package 'INLA' is required for INLA uncertainty quantification.")
      }

      posterior <- tryCatch(
        INLA::inla.posterior.sample(n = n_samples, result = model),
        error = function(e) NULL
      )

      if (is.null(posterior)) {
        stop("INLA posterior samples could not be drawn. Ensure the model was fitted with control.compute = list(config = TRUE).")
      }

      beta_draws <- do.call(rbind, lapply(posterior, function(s) {
        latent <- s$latent
        names(latent) <- gsub(":1$", "", names(latent))
        vals <- latent[cb_names_ref]
        as.numeric(vals)
      }))

      colnames(beta_draws) <- cb_names_ref
      return(beta_draws)
    }

    if (inherits(model, "bdlnm")) {

      beta_draws <- model$coefficients

      if (is.null(dim(beta_draws))) {
        beta_draws <- matrix(beta_draws, ncol = 1)
      }

      if (!all(cb_names_ref %in% rownames(beta_draws))) {
        stop("Mismatch between bdlnm posterior draws and crossbasis structure.")
      }

      beta_draws <- beta_draws[cb_names_ref, , drop = FALSE]

      if (ncol(beta_draws) > n_samples) {
        set.seed(1)
        keep <- sample(seq_len(ncol(beta_draws)), n_samples)
        beta_draws <- beta_draws[, keep, drop = FALSE]
      }

      out <- t(beta_draws)
      colnames(out) <- cb_names_ref
      return(out)
    }

    beta_hat <- cv$beta[cb_names_ref]
    V_hat <- cv$vcov[cb_names_ref, cb_names_ref, drop = FALSE]

    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty approximation.")
    }

    beta_draws <- MASS::mvrnorm(
      n = n_samples,
      mu = beta_hat,
      Sigma = V_hat
    )

    if (is.null(dim(beta_draws))) {
      beta_draws <- matrix(beta_draws, nrow = 1)
    }

    colnames(beta_draws) <- cb_names_ref
    return(beta_draws)
  }

  # ----------------------------------------------------------
  # No model -> raw only
  # ----------------------------------------------------------
  if (is.null(fit)) {

    if (is.null(profile)) {
      stop("When `fit = NULL`, `profile` must be provided.")
    }

    if (!is.numeric(profile)) {
      stop("When `fit = NULL`, `profile` must be a numeric vector.")
    }

    if (any(!is.finite(profile))) {
      stop("`profile` must contain only finite values.")
    }

    if (!is.null(epi_data)) {
      warning("`epi_data`, `group`, and `group_level` are ignored when `fit = NULL`.")
    }

    if (!is.null(var)) {
      warning("`var` is ignored when `fit = NULL`.")
    }

    if (reverse) {
      profile <- rev(profile)
      message(
        "Profile reversed: interpreting input as oldest \u2192 most recent. ",
        "After reversal, first value = lag 0 (most recent); last value = lag max (oldest)."
      )
    } else {
      message(
        "Assuming profile is lag-ordered: first value = lag 0 (most recent); ",
        "last value = lag max (oldest)."
      )
    }

    return(data.frame(
      ECI_raw = sum(profile),
      ECI_weighted = NA_real_,
      scale = scale
    ))
  }

  # ----------------------------------------------------------
  # Extract model metadata
  # ----------------------------------------------------------
  spec <- attr(fit, "epiexposure_spec")
  vars_fit <- attr(fit, "epiexposure_vars")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  family_fit <- attr(fit, "epiexposure_family")

  if (is.null(spec)) stop("`fit` does not contain `epiexposure_spec`.")
  if (is.null(vars_fit) || length(vars_fit) == 0) stop("`fit` does not contain `epiexposure_vars`.")

  # ----------------------------------------------------------
  # Resolve `var`
  # ----------------------------------------------------------
  if (is.null(var)) {
    if (length(vars_fit) == 1) {
      var <- vars_fit[1]
    } else {
      stop(
        "This fitted model contains multiple exposure variables (",
        paste(vars_fit, collapse = ", "),
        "). Please provide `var` using exactly the same variable name(s) stored in the model."
      )
    }
  }

  if (!is.character(var) || length(var) < 1) {
    stop("`var` must be a character scalar or vector when `fit` is provided.")
  }

  bad_var <- setdiff(var, vars_fit)
  if (length(bad_var) > 0) {
    stop(
      "`var` contains name(s) not found in the fitted model: ",
      paste(bad_var, collapse = ", "),
      ". Use exactly the same variable name(s) stored in the model: ",
      paste(vars_fit, collapse = ", "),
      "."
    )
  }

  # ----------------------------------------------------------
  # Decide input strategy
  # ----------------------------------------------------------
  using_profile <- !is.null(profile)
  using_data <- !is.null(epi_data)

  if (using_profile && using_data) {
    stop("Provide either `profile` OR `epi_data`, not both.")
  }

  if (!using_profile && !using_data) {
    stop("When `fit` is provided, you must supply either `profile` OR `epi_data` + `group`.")
  }

  # ----------------------------------------------------------
  # Normalize profile(s) from direct input
  # ----------------------------------------------------------
  if (using_profile) {

    if (is.numeric(profile)) {

      if (length(var) != 1) {
        stop(
          "When `profile` is a numeric vector and `fit` is provided, `var` must have length 1. ",
          "If you want to compute ECI for multiple variables, provide `profile` as a list and `var = c(...)`."
        )
      }

      if (any(!is.finite(profile))) {
        stop("`profile` must contain only finite values.")
      }

      profiles <- setNames(list(profile), var)

    } else if (is.list(profile)) {

      if (length(profile) != length(var) && is.null(names(profile))) {
        stop(
          "When `profile` is an unnamed list, its length must match the length of `var`."
        )
      }

      for (i in seq_along(profile)) {
        if (!is.numeric(profile[[i]]) || any(!is.finite(profile[[i]]))) {
          stop("All elements in `profile` must be finite numeric vectors.")
        }
      }

      if (is.null(names(profile))) {
        profiles <- profile
        names(profiles) <- var
      } else {
        missing_profiles <- setdiff(var, names(profile))
        if (length(missing_profiles) > 0) {
          stop(
            "The following variable(s) in `var` are missing from the named `profile` list: ",
            paste(missing_profiles, collapse = ", "),
            "."
          )
        }
        profiles <- profile[var]
      }

    } else {
      stop("`profile` must be a numeric vector or a list of numeric vectors.")
    }

    if (reverse) {
      profiles <- lapply(profiles, rev)
      message(
        "Profile(s) reversed: interpreting input as oldest \u2192 most recent. ",
        "After reversal, first value = lag 0 (most recent); last value = lag max (oldest)."
      )
    } else {
      message(
        "Assuming supplied profile(s) are lag-ordered: first value = lag 0 (most recent); ",
        "last value = lag max (oldest)."
      )
    }
  }

  # ----------------------------------------------------------
  # Build group levels to use from epi_data
  # ----------------------------------------------------------
  if (using_data) {

    if (!is.data.frame(epi_data)) {
      stop("`epi_data` must be a data.frame.")
    }

    if (is.null(group) || !is.character(group) || length(group) != 1L) {
      stop("When using `epi_data`, `group` must be a single column name.")
    }

    if (!group %in% names(epi_data)) {
      stop("`group` ('", group, "') was not found in `epi_data`.")
    }

    missing_data_vars <- setdiff(var, names(epi_data))
    if (length(missing_data_vars) > 0) {
      stop(
        "The following variable(s) requested in `var` are not present in `epi_data`: ",
        paste(missing_data_vars, collapse = ", "),
        "."
      )
    }

    if (is.null(group_level)) {
      levels_to_use <- unique(epi_data[[group]])
    } else {
      missing_levels <- setdiff(group_level, unique(epi_data[[group]]))
      if (length(missing_levels) > 0) {
        stop(
          "The following `group_level` value(s) were not found inside `epi_data[['", group, "']]`: ",
          paste(missing_levels, collapse = ", "),
          "."
        )
      }
      levels_to_use <- group_level
    }

    if (length(levels_to_use) == 0) {
      stop("No valid group levels were found to compute ECI.")
    }
  }

  # ----------------------------------------------------------
  # Extract deterministic coefficients once
  # ----------------------------------------------------------
  cv <- extract_coef_vcov(fit)
  beta_full <- cv$beta
  linkinv <- get_linkinv(family_fit)

  # ----------------------------------------------------------
  # Worker by variable
  # ----------------------------------------------------------
  compute_one_var <- function(var_sel, profile_sel) {

    spec_v <- spec[[var_sel]]

    if (is.null(spec_v)) {
      stop("Missing spec for variable: ", var_sel)
    }

    lag_max_use <- as.integer(max(spec_v$lag_max))

    if (length(profile_sel) != lag_max_use + 1L) {
      stop(
        "`profile` length for variable '", var_sel,
        "' must be equal to lag_max + 1 (expected ", lag_max_use + 1L, ")."
      )
    }

    eci_raw <- sum(profile_sel)

    cb <- dlnm::crossbasis(
      profile_sel,
      lag = lag_max_use,
      argvar = spec_v$argvar,
      arglag = spec_v$arglag
    )

    cb_row <- as.numeric(cb[lag_max_use + 1L, ])

    cb_names_ref <- get_cb_names(names(beta_full), var_sel, cb_cols_fit)

    if (length(cb_names_ref) != length(cb_row)) {
      stop(
        "Mismatch between coefficient vector and crossbasis structure for variable '",
        var_sel, "'."
      )
    }

    # -----------------------
    # Deterministic
    # -----------------------
    if (!uncertainty) {

      beta_cb <- as.numeric(beta_full[cb_names_ref])

      eta_weighted <- sum(cb_row * beta_cb)

      eci_weighted <- transform_scale(
        x = eta_weighted,
        scale = scale,
        linkinv = linkinv
      )

      out <- data.frame(
        var = var_sel,
        ECI_raw = eci_raw,
        ECI_weighted = eci_weighted,
        scale = scale,
        stringsAsFactors = FALSE
      )

      return(out)
    }

    # -----------------------
    # Uncertainty
    # -----------------------
    beta_draws <- extract_beta_draws(
      model = fit,
      cb_names_ref = cb_names_ref,
      n_samples = n_samples,
      cv = cv
    )

    draw_sd <- apply(beta_draws, 2, stats::sd)
    if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
      warning(
        "Near-zero variability detected in coefficient draws for variable '", var_sel,
        "'. ECI intervals may collapse to a single value."
      )
    }

    eta_draws <- as.numeric(beta_draws %*% cb_row)

    eci_draws <- transform_scale(
      x = eta_draws,
      scale = scale,
      linkinv = linkinv
    )

    if (output == "samples") {
      return(data.frame(
        var = var_sel,
        sample = seq_along(eci_draws),
        ECI_raw = eci_raw,
        ECI_weighted = eci_draws,
        scale = scale,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      var = var_sel,
      ECI_raw = eci_raw,
      ECI_weighted = stats::median(eci_draws, na.rm = TRUE),
      sd = safe_sd(eci_draws),
      lower = safe_quantile(eci_draws)[1],
      upper = safe_quantile(eci_draws)[2],
      scale = scale,
      stringsAsFactors = FALSE
    )
  }

  # ----------------------------------------------------------
  # Execute by input mode
  # ----------------------------------------------------------
  if (using_profile) {

    out_list <- lapply(var, function(v) compute_one_var(v, profiles[[v]]))
    out <- do.call(rbind, out_list)
    rownames(out) <- NULL

    if (length(var) == 1 && is.numeric(profile)) {
      out$var <- NULL
    }

    return(out)
  }

  # ----------------------------------------------------------
  # Execute over all requested group levels
  # ----------------------------------------------------------
  out_by_group <- lapply(levels_to_use, function(gl) {

    data_sub <- epi_data[epi_data[[group]] == gl, , drop = FALSE]

    profiles <- lapply(var, function(v) data_sub[[v]])
    names(profiles) <- var

    if (reverse) {
      profiles <- lapply(profiles, rev)
    }

    tmp_list <- lapply(var, function(v) compute_one_var(v, profiles[[v]]))
    tmp <- do.call(rbind, tmp_list)
    rownames(tmp) <- NULL

    tmp[[group]] <- gl

    # put grouping column first
    tmp <- tmp[, c(group, setdiff(names(tmp), group)), drop = FALSE]

    tmp
  })

  out <- do.call(rbind, out_by_group)
  rownames(out) <- NULL

  return(out)
}
#'
#' This function computes both the raw Exposure Cumulative Impact (`ECI_raw`)
#' and, when a fitted model is provided, the model-weighted Exposure Cumulative
