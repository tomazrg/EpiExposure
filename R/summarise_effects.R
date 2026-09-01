#' Summarises distributed lag nonlinear model effects
#'
#' Returns lag-specific or period-specific effect summaries based on the
#' exposure specification stored by `fit_epidlnm()`. Exposure histories are
#' interpreted chronologically: the first observation is the oldest exposure
#' and the final observation is the most recent exposure, corresponding
#' internally to lag 0.
#'
#' For non-ordinal models, effects are calculated relative to a reference
#' (baseline) exposure condition. The `"linear"` scale returns the contrast
#' on the linear predictor scale (`eta`). The `"exponentiated"` scale returns
#' `exp(eta)`, and the `"percent"` scale returns `(exp(eta) - 1) * 100`.
#'
#' Effects therefore describe how strongly a given exposure condition is
#' associated with the outcome relative to the selected reference exposure.
#' Values greater than zero indicate a greater relative effect than the
#' reference condition, whereas values below zero indicate a smaller relative
#' effect.
#'
#' Under a log link, exponentiated and percent effects can be interpreted as
#' relative changes in the expected outcome. Under a logit link, they represent
#' relative exposure effects with respect to the reference condition and should
#' not be interpreted as direct changes in the response variable.
#'
#' For ordinal `brms` models, the function returns the predicted probability of
#' every response category. `effect_measure` is not applied to ordinal category
#' probabilities.
#'
#' The `delta` column is defined consistently for every non-ordinal family as
#' `predicted - baseline` on the response scale. Its interpretation follows the
#' response distribution: a difference in predicted means for Gaussian and
#' Gamma models, a difference in expected counts for Poisson and Negative
#' Binomial models, and a difference in predicted proportions or probabilities
#' for Beta and Binomial models. For Beta and Binomial models, multiplying
#' `delta` by 100 gives the difference in percentage points.
#'
#' @param fit Fitted model returned by `fit_epidlnm()`.
#' @param data Long-format exposure data containing `epi_id`, `time`, and the
#'   fitted exposure variables.
#' @param var Exposure-variable name(s), or `NULL` for all fitted exposures.
#' @param scale Character. `"lag"` for lag-specific effects or
#'   `"period"` for cumulative effects aggregated across lag intervals.
#' @param lag_periods Optional data.frame with `period`, `lag_start`,
#'   and `lag_end`, required for user-defined period summaries.
#' @param probs Unique probabilities used to generate the exposure grid for
#'   variables without a custom `at` specification. The default preserves the
#'   current quantile-based behavior.
#' @param at Optional exposure values at which effects are evaluated. Use
#'   `NULL` to generate quantile-based grids from `probs` for every variable.
#'   A numeric vector is allowed when exactly one exposure variable is
#'   requested. For multiple variables, use a named list, for example
#'   `list(tmax = seq(29, 34, length.out = 1000), rain = seq(0, 50, by = 0.5))`.
#'   The list may be partial: variables omitted from `at` continue to use the
#'   quantile-based grid defined by `probs`. To evaluate effects at the unique
#'   observed values, supply `sort(unique(data[[variable]]))` for that variable.
#' @param ref List with `method = "median"`, `"percentile"`, or `"fixed"`, and
#'   an optional `value`.
#' @param effect_measure Character. `"linear"`, `"exponentiated"`, or
#'   `"percent"`.
#' @param incremental Logical. For `scale = "period"`, calculate
#'   cumulative effects incrementally from lag 0 through each lag
#'   instead of using user-defined periods.
#' @param uncertainty Logical. Propagate coefficient uncertainty draw by draw.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Positive integer number of draws. At least two are required
#'   when `uncertainty = TRUE`.
#' @param diagnostics Logical. Available only for lag-specific non-ordinal
#'   effects. Returns integrated lag rankings and contribution metrics.
#'   Rankings depend on `effect_measure`.
#' @param seed Optional finite integer seed for reproducible sampling.
#'
#' @details
#' The cross-basis is reconstructed exclusively from `epiexposure_spec`; basis
#' arguments are not re-specified by the user. When `uncertainty = FALSE`, the
#' central coefficients supplied by the fitted engine are used. When
#' `uncertainty = TRUE`, every draw is transformed separately and summaries use
#' the median and empirical 2.5% and 97.5% quantiles. The baseline response is
#' intentionally held fixed, consistently with other EpiExposure functions.
#'
#' @return A data.frame of effects or category probabilities. For non-ordinal
#'   models, `predicted` is the response-scale prediction, `baseline` is the
#'   response-scale reference prediction, and `delta` is their absolute
#'   difference (`predicted - baseline`) in the natural units of the response.
#'   With `diagnostics = TRUE`, a list with `effects` and `diagnostics`.
#' @export
summarise_effects <- function(
    fit,
    data,
    var = NULL,
    scale = c("lag", "period"),
    lag_periods = NULL,
    probs = seq(0.05, 0.95, by = 0.01),
    at = NULL,
    ref = list(method = "median", value = NULL),
    effect_measure = c("linear", "exponentiated", "percent"),
    incremental = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    diagnostics = FALSE,
    seed = NULL
) {
  scale <- match.arg(scale)
  effect_measure <- match.arg(effect_measure)
  output <- match.arg(output)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!all(c("epi_id", "time") %in% names(data))) {
    stop("`data` must contain at least 'epi_id' and 'time'.")
  }
  if (!is.logical(incremental) || length(incremental) != 1L || is.na(incremental)) {
    stop("`incremental` must be TRUE or FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L || is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.logical(diagnostics) || length(diagnostics) != 1L || is.na(diagnostics)) {
    stop("`diagnostics` must be TRUE or FALSE.")
  }
  if (diagnostics && scale != "lag") {
    stop("`diagnostics = TRUE` is available only when `scale = 'lag'`.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L ||
      !is.finite(n_samples) || n_samples <= 0 ||
      n_samples != as.integer(n_samples)) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)
  if (uncertainty && n_samples < 2L) {
    stop("`n_samples` must be at least 2 when `uncertainty = TRUE`.")
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
        seed != as.integer(seed)) {
      stop("`seed` must be NULL or one finite integer.")
    }
    set.seed(as.integer(seed))
  }
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop("`probs` must contain finite probabilities between 0 and 1.")
  }
  if (anyDuplicated(probs)) stop("`probs` must contain unique probabilities.")
  probs <- sort(probs)
  if (!is.null(at) && !is.numeric(at) && !is.list(at)) {
    stop("`at` must be NULL, a numeric vector, or a named list.")
  }
  if (!is.list(ref) || is.null(ref$method) ||
      !is.character(ref$method) || length(ref$method) != 1L ||
      is.na(ref$method) || !nzchar(ref$method)) {
    stop("`ref` must contain one valid character `method`.")
  }
  ref_method <- match.arg(ref$method, c("median", "percentile", "fixed"))
  if (scale == "period" && !incremental && is.null(lag_periods)) {
    stop("`lag_periods` is required when scale = 'period' and incremental = FALSE.")
  }
  if (!is.null(lag_periods)) {
    required <- c("period", "lag_start", "lag_end")
    if (!is.data.frame(lag_periods) || !all(required %in% names(lag_periods))) {
      stop("`lag_periods` must be a data.frame containing period, lag_start, and lag_end.")
    }
    if (anyNA(lag_periods$period) || any(lag_periods$period == "") ||
        anyDuplicated(lag_periods$period)) {
      stop("`lag_periods$period` must contain unique non-empty labels.")
    }
    for (column in c("lag_start", "lag_end")) {
      x <- lag_periods[[column]]
      if (!is.numeric(x) || any(!is.finite(x)) || any(x != as.integer(x)) || any(x < 0)) {
        stop("`lag_periods$", column, "` must contain non-negative finite integers.")
      }
    }
    if (any(lag_periods$lag_start > lag_periods$lag_end)) {
      stop("Every lag period must satisfy `lag_start <= lag_end`.")
    }
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b
  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    x <- x[is.finite(x)]
    if (!length(x)) return(rep(NA_real_, length(probs)))
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }
  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(0)
    stats::sd(x)
  }
  sort_cb_names <- function(x) {
    if (!length(x)) return(x)
    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    missing_idx <- is.na(idx)
    idx[missing_idx] <- seq_along(x)[missing_idx]
    x[order(idx)]
  }
  match_brms_draw_names <- function(cb_names_ref, draw_names) {
    prefixed <- paste0("b_", cb_names_ref)
    if (all(prefixed %in% draw_names)) return(prefixed)
    if (all(cb_names_ref %in% draw_names)) return(cb_names_ref)
    stop("Could not match brms posterior draws to cross-basis coefficients.")
  }
  get_cb_names <- function(coefficient_names, variable, expected, stored = NULL) {
    matched <- character(0)
    if (!is.null(stored)) {
      ordered <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), stored, value = TRUE
      ))
      matched <- ordered[ordered %in% coefficient_names]
    }
    if (!length(matched)) {
      matched <- sort_cb_names(grep(
        paste0("^cb_", variable, "_"), coefficient_names, value = TRUE
      ))
    }
    if (length(matched) != expected) {
      stop(
        "Could not align coefficients with the cross-basis for variable '",
        variable, "'. Expected ", expected, " coefficients but found ",
        length(matched), "."
      )
    }
    matched
  }
  check_vcov_names <- function(V, coefficient_names) {
    if (is.null(rownames(V)) || is.null(colnames(V))) {
      stop("The coefficient covariance matrix has no row or column names.")
    }
    missing_names <- setdiff(
      coefficient_names,
      intersect(rownames(V), colnames(V))
    )
    if (length(missing_names)) {
      stop(
        "Cross-basis coefficients missing from the covariance matrix: ",
        paste(missing_names, collapse = ", "), "."
      )
    }
    invisible(TRUE)
  }
  check_lag_coverage <- function(dat, maximum_lag) {
    required_n <- maximum_lag + 1L
    temporal_summary <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(
        n_rows = dplyr::n(),
        n_time = dplyr::n_distinct(time),
        has_missing_time = any(is.na(time)),
        .groups = "drop"
      )
    missing_time <- temporal_summary |> dplyr::filter(has_missing_time)
    if (nrow(missing_time)) {
      stop("Missing `time` values detected. Example epi_id: ",
           paste(utils::head(missing_time$epi_id, 5), collapse = ", "), ".")
    }
    duplicate_time <- temporal_summary |> dplyr::filter(n_rows != n_time)
    if (nrow(duplicate_time)) {
      stop("Duplicated `time` values detected within epidemics. Example epi_id: ",
           paste(utils::head(duplicate_time$epi_id, 5), collapse = ", "), ".")
    }
    insufficient <- temporal_summary |> dplyr::filter(n_time < required_n)
    if (nrow(insufficient)) {
      stop(
        "Insufficient temporal coverage. Required observations per epi_id: ",
        required_n, ". Example epi_id: ",
        paste(utils::head(insufficient$epi_id, 5), collapse = ", "), "."
      )
    }
    invisible(TRUE)
  }
  build_pooled_series <- function(dat, variable, separator_n) {
    ids <- unique(dat$epi_id)
    pooled <- vector("list", length(ids))
    for (i in seq_along(ids)) {
      values <- dat |>
        dplyr::filter(epi_id == ids[i]) |>
        dplyr::arrange(time) |>
        dplyr::pull(.data[[variable]])
      pooled[[i]] <- c(values, rep(NA_real_, separator_n))
    }
    unlist(pooled)
  }

  extract_coef_vcov <- function(model) {
    if (inherits(model, "glmmTMB")) return(list(
      beta = glmmTMB::fixef(model)$cond,
      vcov = as.matrix(stats::vcov(model)$cond)
    ))
    if (inherits(model, "merMod")) return(list(
      beta = lme4::fixef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "lme")) return(list(
      beta = nlme::fixef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "gls")) return(list(
      beta = stats::coef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (is.list(model) && !is.null(model$gam) && inherits(model$gam, "gam")) {
      return(list(
        beta = stats::coef(model$gam),
        vcov = as.matrix(stats::vcov(model$gam))
      ))
    }
    if (inherits(model, "glm") || inherits(model, "gam")) return(list(
      beta = stats::coef(model), vcov = as.matrix(stats::vcov(model))
    ))
    if (inherits(model, "HLfit")) {
      V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
      if (is.null(V)) stop("Could not extract vcov from the spaMM model.")
      return(list(beta = spaMM::fixef(model), vcov = V))
    }
    if (inherits(model, "brmsfit")) {
      fixed <- brms::fixef(model)
      beta <- fixed[, "Estimate"]
      names(beta) <- rownames(fixed)
      return(list(beta = beta, vcov = as.matrix(stats::vcov(model))))
    }
    if (inherits(model, "inla")) {
      beta <- model$summary.fixed$mean
      if (is.null(names(beta))) names(beta) <- rownames(model$summary.fixed)
      V <- diag(model$summary.fixed$sd^2)
      dimnames(V) <- list(names(beta), names(beta))
      return(list(beta = beta, vcov = V))
    }
    if (inherits(model, "bdlnm")) {
      if (is.null(model$coefficients.summary) || is.null(model$coefficients)) {
        stop("The bdlnm model lacks coefficient summaries or posterior draws.")
      }
      beta <- model$coefficients.summary[, "mean"]
      if (is.null(names(beta))) names(beta) <- rownames(model$coefficients.summary)
      V <- stats::cov(t(model$coefficients))
      if (is.null(rownames(V)) && !is.null(names(beta))) dimnames(V) <- list(names(beta), names(beta))
      return(list(beta = beta, vcov = V))
    }
    beta <- stats::coef(model)
    if (!is.numeric(beta)) stop("Could not extract numeric model coefficients.")
    V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
    if (is.null(V)) stop("Could not extract the coefficient covariance matrix.")
    list(beta = beta, vcov = V)
  }

  extract_beta_draws <- function(model, names_ref, n, coefficient_info) {
    if (inherits(model, "brmsfit")) {
      if (!requireNamespace("posterior", quietly = TRUE)) {
        stop("Package 'posterior' is required for brms uncertainty.")
      }
      draws <- posterior::as_draws_matrix(model)
      draw_names <- match_brms_draw_names(names_ref, colnames(draws))
      draws <- as.matrix(draws[, draw_names, drop = FALSE])
      if (nrow(draws) > n) {
        draws <- draws[sample(seq_len(nrow(draws)), n), , drop = FALSE]
      }
      colnames(draws) <- names_ref
      return(draws)
    }
    if (inherits(model, "inla")) {
      if (!requireNamespace("INLA", quietly = TRUE)) {
        stop("Package 'INLA' is required for INLA uncertainty.")
      }
      posterior_samples <- tryCatch(
        INLA::inla.posterior.sample(n = n, result = model),
        error = function(e) NULL
      )
      if (is.null(posterior_samples)) {
        stop("Could not draw INLA posterior samples. Fit with `control.compute = list(config = TRUE)`.")
      }
      draws <- do.call(rbind, lapply(posterior_samples, function(sample_object) {
        latent <- sample_object$latent
        names(latent) <- gsub(":1$", "", names(latent))
        values <- latent[names_ref]
        if (anyNA(values)) {
          stop("Could not match INLA draws to cross-basis coefficients.")
        }
        as.numeric(values)
      }))
      colnames(draws) <- names_ref
      return(draws)
    }
    if (inherits(model, "bdlnm")) {
      draws <- model$coefficients
      if (is.null(dim(draws))) draws <- matrix(draws, ncol = 1L)
      if (is.null(rownames(draws)) || !all(names_ref %in% rownames(draws))) {
        stop("Could not match bdlnm draws to cross-basis coefficients.")
      }
      draws <- draws[names_ref, , drop = FALSE]
      if (ncol(draws) > n) {
        draws <- draws[, sample(seq_len(ncol(draws)), n), drop = FALSE]
      }
      draws <- t(draws)
      colnames(draws) <- names_ref
      return(draws)
    }
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package 'MASS' is required for frequentist uncertainty.")
    }
    check_vcov_names(coefficient_info$vcov, names_ref)
    draws <- MASS::mvrnorm(
      n = n,
      mu = coefficient_info$beta[names_ref],
      Sigma = coefficient_info$vcov[names_ref, names_ref, drop = FALSE]
    )
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1L)
    colnames(draws) <- names_ref
    draws
  }

  get_linkfun <- function(link_name) {
    switch(
      tolower(link_name),
      identity = identity,
      log = log,
      logit = stats::qlogis,
      probit = stats::qnorm,
      cloglog = function(mu) log(-log1p(-mu)),
      inverse = function(mu) 1 / mu,
      stop("Unsupported link function for link '", link_name, "'.")
    )
  }
  get_linkinv <- function(link_name) {
    switch(
      tolower(link_name),
      identity = identity,
      log = exp,
      logit = stats::plogis,
      probit = stats::pnorm,
      cloglog = function(eta) 1 - exp(-exp(eta)),
      inverse = function(eta) 1 / eta,
      stop("Unsupported inverse-link function for link '", link_name, "'.")
    )
  }

  # Integrated lag diagnostics -------------------------------------------
  # The input must be deterministic lag effects or sample-level lag
  # effects. A single diagnostics table combines lag ranking and contribution
  # metrics, so the returned object has only two top-level components.
  compute_lag_diagnostics <- function(lag_df) {
    if (!all(c("var", "lag", "effect") %in% names(lag_df))) {
      stop("Internal diagnostics require columns 'var', 'lag', and 'effect'.")
    }
    if (!is.numeric(lag_df$effect)) {
      stop("The internal lag `effect` column must be numeric.")
    }

    compute_one_diagnostic <- function(df) {
      lag_metrics <- df |>
        dplyr::group_by(lag) |>
        dplyr::summarise(
          mean_effect = mean(effect, na.rm = TRUE),
          max_effect = {
            valid <- effect[is.finite(effect)]
            if (!length(valid)) NA_real_ else valid[which.max(abs(valid))]
          },
          max_abs_effect = {
            valid <- effect[is.finite(effect)]
            if (!length(valid)) NA_real_ else max(abs(valid))
          },
          absmean_effect = mean(abs(effect), na.rm = TRUE),
          eta_lag = sum(effect, na.rm = TRUE),
          .groups = "drop"
        )

      lag_metrics <- lag_metrics |>
        dplyr::mutate(
          score_max = max_abs_effect,
          score_mean = abs(mean_effect),
          score_absmean = absmean_effect,
          rank_max = dplyr::min_rank(dplyr::desc(score_max)),
          rank_mean = dplyr::min_rank(dplyr::desc(score_mean)),
          rank_absmean = dplyr::min_rank(dplyr::desc(score_absmean))
        )

      absolute_denominator <- sum(abs(lag_metrics$eta_lag), na.rm = TRUE)
      signed_denominator <- sum(lag_metrics$eta_lag, na.rm = TRUE)
      stability_tolerance <- sqrt(.Machine$double.eps) *
        max(1, absolute_denominator)
      signed_stable <- is.finite(signed_denominator) &&
        abs(signed_denominator) > stability_tolerance

      lag_metrics |>
        dplyr::mutate(
          contribution_absolute = if (is.finite(absolute_denominator) &&
                                      absolute_denominator > 0) {
            abs(eta_lag) / absolute_denominator
          } else {
            NA_real_
          },
          contribution_absolute_percent = 100 * contribution_absolute,
          contribution_signed = if (signed_stable) {
            eta_lag / signed_denominator
          } else {
            NA_real_
          },
          contribution_signed_percent = 100 * contribution_signed,
          signed_contribution_available = signed_stable
        )
    }

    has_samples <- "sample" %in% names(lag_df)

    if (!has_samples) {
      diagnostics_out <- lag_df |>
        dplyr::group_by(var) |>
        dplyr::group_modify(~compute_one_diagnostic(.x)) |>
        dplyr::ungroup() |>
        dplyr::arrange(var, rank_absmean, lag)
      return(as.data.frame(diagnostics_out))
    }

    sample_diagnostics <- lag_df |>
      dplyr::group_by(sample, var) |>
      dplyr::group_modify(~compute_one_diagnostic(.x)) |>
      dplyr::ungroup()

    numeric_metrics <- c(
      "mean_effect", "max_effect", "max_abs_effect", "absmean_effect",
      "eta_lag", "score_max", "score_mean", "score_absmean",
      "rank_max", "rank_mean", "rank_absmean",
      "contribution_absolute", "contribution_absolute_percent",
      "contribution_signed", "contribution_signed_percent"
    )

    diagnostics_out <- sample_diagnostics |>
      dplyr::group_by(var, lag) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(numeric_metrics),
          list(
            estimate = ~stats::median(.x, na.rm = TRUE),
            sd = ~safe_sd(.x),
            lower = ~safe_quantile(.x)[1],
            upper = ~safe_quantile(.x)[2]
          ),
          .names = "{.col}_{.fn}"
        ),
        signed_contribution_available = all(
          signed_contribution_available,
          na.rm = TRUE
        ),
        .groups = "drop"
      )

    for (metric in numeric_metrics) {
      estimate_column <- paste0(metric, "_estimate")
      names(diagnostics_out)[names(diagnostics_out) == estimate_column] <- metric
    }

    diagnostics_out <- diagnostics_out |>
      dplyr::arrange(var, rank_absmean, lag)
    as.data.frame(diagnostics_out)
  }

  fit_spec <- attr(fit, "epiexposure_spec")
  fit_vars <- attr(fit, "epiexposure_vars")
  family_fit <- attr(fit, "epiexposure_family")
  family_name <- attr(fit, "epiexposure_family_name")
  link_name <- attr(fit, "epiexposure_link")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  data_template <- attr(fit, "epiexposure_data_template")
  basis_objects <- attr(fit, "epiexposure_basis_objects")

  if (is.null(family_name) || !is.character(family_name) || length(family_name) != 1L) {
    stop("`fit` does not contain valid `epiexposure_family_name` metadata.")
  }
  family_name <- tolower(family_name)
  if (
    family_name == "ordinal" &&
    effect_measure != "linear"
  ) {
    warning(
      "`effect_measure = '", effect_measure,
      "' is ignored for ordinal models. ",
      "Category probabilities are returned instead.",
      call. = FALSE
    )
  }
  if (is.null(link_name) || !is.character(link_name) || length(link_name) != 1L) {
    stop("`fit` does not contain valid `epiexposure_link` metadata.")
  }
  link_name <- tolower(link_name)
  if ((is.null(fit_vars) || !length(fit_vars)) && is.list(basis_objects) && length(basis_objects)) {
    fit_vars <- names(basis_objects)
  }
  if (family_name == "ordinal" && diagnostics) {
    stop("`diagnostics = TRUE` is not available for ordinal category probabilities.")
  }
  if (is.null(var)) {
    variables <- fit_vars
    if (is.null(variables) || !length(variables)) {
      stop("`var` is NULL and no exposure variables were found in model metadata.")
    }
  } else {
    if (!is.character(var) || !length(var) || anyNA(var) || any(var == "")) {
      stop("`var` must be NULL or a non-empty character vector.")
    }
    variables <- var
  }
  if (anyDuplicated(variables)) stop("`var` must contain unique variable names.")
  unknown_model_vars <- if (is.null(fit_vars)) character(0) else setdiff(variables, fit_vars)
  if (length(unknown_model_vars)) {
    stop("Variables not found in fitted-model metadata: ",
         paste(unknown_model_vars, collapse = ", "), ".")
  }
  missing_data_vars <- setdiff(variables, names(data))
  if (length(missing_data_vars)) {
    stop("Variables not present in `data`: ",
         paste(missing_data_vars, collapse = ", "), ".")
  }

  # Validate and standardise optional exposure-specific evaluation grids.
  # `at = NULL` preserves the original quantile-based behaviour. A partial
  # named list overrides only the listed variables; all others use `probs`.
  if (is.numeric(at)) {
    if (length(variables) != 1L) {
      stop(
        "A numeric `at` can be used only when exactly one exposure variable ",
        "is requested. For multiple variables, supply `at` as a named list."
      )
    }
    if (!length(at) || any(!is.finite(at))) {
      stop("A numeric `at` must contain at least one finite value.")
    }
    if (anyDuplicated(at)) {
      stop("A numeric `at` must contain unique values.")
    }
    at <- sort(as.numeric(at))
  } else if (is.list(at)) {
    if (!length(at) || is.null(names(at)) || anyNA(names(at)) ||
        any(names(at) == "") || anyDuplicated(names(at))) {
      stop(
        "When supplied as a list, `at` must be a non-empty named list ",
        "with unique, non-empty variable names."
      )
    }
    unknown_at_variables <- setdiff(names(at), variables)
    if (length(unknown_at_variables)) {
      stop(
        "`at` contains variables that were not requested: ",
        paste(unknown_at_variables, collapse = ", "), "."
      )
    }
    for (current_variable in names(at)) {
      current_values <- at[[current_variable]]
      if (!is.numeric(current_values) || !length(current_values) ||
          any(!is.finite(current_values))) {
        stop(
          "`at[['", current_variable,
          "']]` must contain at least one finite numeric value."
        )
      }
      if (anyDuplicated(current_values)) {
        stop(
          "`at[['", current_variable,
          "']]` must contain unique values."
        )
      }
      at[[current_variable]] <- sort(as.numeric(current_values))
    }
  }

  transform_effect <- function(eta) {
    if (effect_measure == "linear") return(eta)
    if (effect_measure == "exponentiated") return(exp(eta))
    (exp(eta) - 1) * 100
  }

  summarise_ordinal_variable <- function(variable, variable_spec, at_values, center_value) {
    if (!inherits(fit, "brmsfit")) {
      stop("Ordinal category probabilities currently require a `brmsfit` model.")
    }
    if (!requireNamespace("brms", quietly = TRUE)) stop("Package 'brms' is required.")

    make_reference_profiles <- function() {
      profiles <- lapply(fit_vars, function(v) {
        lv <- as.integer(max(fit_spec[[v]]$max_lag))
        rep(stats::median(data[[v]], na.rm = TRUE), lv + 1L)
      })
      names(profiles) <- fit_vars
      profiles
    }
    make_newdata <- function(profiles) {
      nd <- data_template[1, , drop = FALSE]
      for (v in fit_vars) {
        sp <- fit_spec[[v]]
        lv <- as.integer(max(sp$max_lag))
        cbv <- dlnm::crossbasis(profiles[[v]], lag = lv, argvar = sp$argvar, arglag = sp$arglag)
        row <- as.numeric(cbv[lv + 1L, , drop = TRUE])
        stored <- sort_cb_names(grep(paste0("^cb_", v, "_"), cb_cols_fit, value = TRUE))
        if (!length(stored)) stored <- sort_cb_names(grep(paste0("^cb_", v, "_"), names(data_template), value = TRUE))
        if (length(stored) != length(row)) stop("Could not align ordinal cross-basis columns for '", v, "'.")
        nd[stored] <- as.data.frame(as.list(stats::setNames(row, stored)))
      }
      nd
    }
    conditions <- list(); metadata <- list(); k <- 0L
    lag_values <- 0:as.integer(max(variable_spec$max_lag))
    refs <- make_reference_profiles()
    if (scale == "lag") {
      for (value in at_values) for (lag in lag_values) {
        k <- k + 1L; pr <- refs
        pos <- length(pr[[variable]]) - lag
        pr[[variable]][pos] <- value
        conditions[[k]] <- make_newdata(pr)
        metadata[[k]] <- data.frame(var=variable, lag=lag, scale="lag", value=value)
      }
    } else if (incremental) {
      for (value in at_values) for (lag in lag_values) {
        k <- k + 1L; pr <- refs
        positions <- length(pr[[variable]]) - (0:lag)
        pr[[variable]][positions] <- value
        conditions[[k]] <- make_newdata(pr)
        metadata[[k]] <- data.frame(var=variable, lag=lag, period=paste0("0-",lag), scale="period", value=value)
      }
    } else {
      for (value in at_values) for (i in seq_len(nrow(lag_periods))) {
        k <- k + 1L; pr <- refs
        lags <- lag_periods$lag_start[i]:lag_periods$lag_end[i]
        positions <- length(pr[[variable]]) - lags
        pr[[variable]][positions] <- value
        conditions[[k]] <- make_newdata(pr)
        metadata[[k]] <- data.frame(var=variable, period=lag_periods$period[i], scale="period", value=value)
      }
    }
    nd <- do.call(rbind, conditions)
    ep <- brms::posterior_epred(
      fit, newdata = nd, re_formula = NA,
      ndraws = if (uncertainty) n_samples else NULL
    )
    d <- dim(ep)
    if (length(d) != 3L) stop("Expected ordinal probabilities with dimensions draws x observations x categories.")
    category_names <- dimnames(ep)[[3]]
    if (is.null(category_names)) category_names <- as.character(seq_len(d[3]))
    meta <- do.call(rbind, metadata); rownames(meta) <- NULL
    if (!uncertainty) {
      probs_mean <- apply(ep, c(2,3), mean, na.rm=TRUE)
      out <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
        cbind(meta[i,,drop=FALSE], category=category_names, probability=as.numeric(probs_mean[i,]))
      }))
      return(list(effects=as.data.frame(out), diagnostics_source=NULL))
    }
    samples <- do.call(rbind, lapply(seq_len(d[1]), function(draw) {
      do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
        cbind(sample=draw, meta[i,,drop=FALSE], category=category_names, probability=as.numeric(ep[draw,i,]))
      }))
    }))
    if (output == "samples") return(list(effects=as.data.frame(samples), diagnostics_source=NULL))
    groups <- intersect(c("var","lag","period","scale","value","category"), names(samples))
    summary <- samples |>
      dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |>
      dplyr::summarise(
        probability=stats::median(probability,na.rm=TRUE),
        probability_sd=safe_sd(probability),
        probability_lower=safe_quantile(probability)[1],
        probability_upper=safe_quantile(probability)[2],
        .groups="drop"
      )
    list(effects=as.data.frame(summary), diagnostics_source=NULL)
  }

  summarise_one_variable <- function(variable) {
    if (is.null(fit_spec) || !is.list(fit_spec) || is.null(fit_spec[[variable]])) {
      stop("No valid `epiexposure_spec` metadata for variable '", variable, "'.")
    }
    variable_spec <- fit_spec[[variable]]
    if (is.null(variable_spec$max_lag) || is.null(variable_spec$argvar) ||
        is.null(variable_spec$arglag)) {
      stop("Incomplete exposure specification for variable '", variable, "'.")
    }
    max_lag_use <- as.integer(max(variable_spec$max_lag))
    argvar <- variable_spec$argvar
    arglag <- variable_spec$arglag

    if (!is.null(lag_periods) && any(lag_periods$lag_end > max_lag_use)) {
      stop("`lag_periods` exceeds the maximum lag for variable '", variable, "'.")
    }
    check_lag_coverage(data, max_lag_use)
    pooled <- build_pooled_series(data, variable, max_lag_use)
    cb <- dlnm::crossbasis(
      pooled,
      lag = max_lag_use,
      argvar = argvar,
      arglag = arglag
    )

    coefficient_info <- extract_coef_vcov(fit)
    coefficients <- coefficient_info$beta
    covariance <- coefficient_info$vcov
    if (is.null(names(coefficients))) {
      stop("The fitted coefficient vector has no names.")
    }

    stored_columns <- cb_cols_fit
    if (is.null(stored_columns) && !is.null(data_template)) {
      stored_columns <- grep("^cb_", names(data_template), value = TRUE)
    }
    cb_names_ref <- get_cb_names(
      coefficient_names = names(coefficients),
      variable = variable,
      expected = ncol(cb),
      stored = stored_columns
    )
    check_vcov_names(covariance, cb_names_ref)
    beta <- coefficients[cb_names_ref]
    vcov_cb <- covariance[cb_names_ref, cb_names_ref, drop = FALSE]

    x_all <- data[[variable]]
    if (!is.numeric(x_all) || all(!is.finite(x_all))) {
      stop("Exposure variable '", variable, "' must contain finite numeric values.")
    }
    variable_at <- NULL
    if (is.numeric(at)) {
      variable_at <- at
    } else if (is.list(at) && variable %in% names(at)) {
      variable_at <- at[[variable]]
    }

    if (is.null(variable_at)) {
      at_values <- sort(unique(as.numeric(stats::quantile(
        x_all,
        probs = probs,
        na.rm = TRUE,
        names = FALSE
      ))))
    } else {
      at_values <- sort(unique(as.numeric(variable_at)))
    }

    if (!length(at_values) || any(!is.finite(at_values))) {
      stop(
        "No valid exposure evaluation values were available for variable '",
        variable, "'."
      )
    }

    observed_range <- range(x_all, na.rm = TRUE)
    outside_observed_range <- at_values < observed_range[1] |
      at_values > observed_range[2]
    if (any(outside_observed_range)) {
      stop(
        "`at` contains values outside the observed range for variable '",
        variable, "'. Observed range: [",
        format(observed_range[1]), ", ", format(observed_range[2]), "]."
      )
    }

    center_value <- switch(
      ref_method,
      median = stats::median(x_all, na.rm = TRUE),
      percentile = {
        if (is.null(ref$value) || length(ref$value) != 1L ||
            !is.numeric(ref$value) || ref$value < 0 || ref$value > 1) {
          stop("For ref$method = 'percentile', `ref$value` must be one probability.")
        }
        stats::quantile(x_all, ref$value, na.rm = TRUE)
      },
      fixed = {
        if (is.null(ref$value) || length(ref$value) != 1L ||
            !is.numeric(ref$value) || !is.finite(ref$value)) {
          stop("For ref$method = 'fixed', `ref$value` must be one finite number.")
        }
        ref$value
      },
      stop("Invalid ref$method. Use 'median', 'percentile', or 'fixed'.")
    )
    center_value <- as.numeric(center_value)

    if (family_name == "ordinal") {
      return(summarise_ordinal_variable(variable, variable_spec, at_values, center_value))
    }

    point_prediction <- dlnm::crosspred(
      cb,
      coef = beta,
      vcov = vcov_cb,
      at = at_values,
      cen = center_value,
      bylag = 1
    )
    point_eta <- point_prediction$matfit
    lag_index <- suppressWarnings(as.integer(gsub("lag", "", colnames(point_eta))))
    if (anyNA(lag_index)) lag_index <- 0:(ncol(point_eta) - 1L)
    lag_order <- order(lag_index)
    lag_index <- lag_index[lag_order]
    point_eta <- point_eta[, lag_order, drop = FALSE]

    linkfun <- get_linkfun(link_name)
    linkinv <- get_linkinv(link_name)
    baseline_response <- NA_real_
    if (!is.null(fit_vars) && !is.null(fit_spec) &&
        all(fit_vars %in% names(data))) {
      reference_profiles <- lapply(fit_vars, function(current_variable) {
        current_spec <- fit_spec[[current_variable]]
        current_lag <- as.integer(max(current_spec$max_lag))
        reference_value <- stats::median(data[[current_variable]], na.rm = TRUE)
        rep(reference_value, current_lag + 1L)
      })
      names(reference_profiles) <- fit_vars
      baseline_response <- tryCatch(
        as.numeric(predict_outcome(
          fit = fit,
          profiles = reference_profiles,
          re = "population",
          type = "response",
          uncertainty = FALSE
        )$prediction[1]),
        error = function(e) {
          warning("Baseline prediction failed: ", conditionMessage(e), call. = FALSE)
          NA_real_
        }
      )
    }

    build_output <- function(eta_matrix) {
      if (scale == "lag") {
        transformed_eta <- eta_matrix
        grid <- expand.grid(
          value = at_values,
          lag = lag_index,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        result <- data.frame(
          var = variable,
          lag = grid$lag,
          scale = "lag",
          value = grid$value,
          eta = as.vector(transformed_eta),
          stringsAsFactors = FALSE
        )
      } else if (incremental) {
        transformed_eta <- t(apply(eta_matrix, 1, cumsum))
        grid <- expand.grid(
          value = at_values,
          lag = lag_index,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        result <- data.frame(
          var = variable,
          lag = grid$lag,
          period = paste0("0-", grid$lag),
          scale = "period",
          value = grid$value,
          eta = as.vector(transformed_eta),
          stringsAsFactors = FALSE
        )
      } else {
        required_period_columns <- c("period", "lag_start", "lag_end")
        if (!all(required_period_columns %in% names(lag_periods))) {
          stop("`lag_periods` must contain period, lag_start, and lag_end.")
        }
        return(purrr::map_dfr(seq_len(nrow(lag_periods)), function(i) {
          selected_lags <- which(
            lag_index >= lag_periods$lag_start[i] &
              lag_index <= lag_periods$lag_end[i]
          )
          if (!length(selected_lags)) {
            stop(
              "No lags found for period '", lag_periods$period[i],
              "' (", lag_periods$lag_start[i], "-", lag_periods$lag_end[i], ")."
            )
          }
          period_eta <- rowSums(eta_matrix[, selected_lags, drop = FALSE])
          period_result <- data.frame(
            var = variable,
            period = lag_periods$period[i],
            scale = "period",
            value = at_values,
            eta = period_eta,
            stringsAsFactors = FALSE
          )
          period_result$effect <- transform_effect(period_result$eta)
          period_result$baseline <- baseline_response
          if (is.na(baseline_response)) {
            period_result$predicted <- NA_real_
            period_result$delta <- NA_real_
          } else {
            period_result$predicted <- linkinv(
              linkfun(baseline_response) + period_result$eta
            )
            period_result$delta <- period_result$predicted - baseline_response
          }
          period_result
        }))
      }

      result$effect <- transform_effect(result$eta)
      result$baseline <- baseline_response
      if (is.na(baseline_response)) {
        result$predicted <- NA_real_
        result$delta <- NA_real_
      } else {
        result$predicted <- linkinv(linkfun(baseline_response) + result$eta)
        result$delta <- result$predicted - baseline_response
      }
      result
    }

    point_output <- build_output(point_eta)

    if (!uncertainty) {
      return(list(
        effects = point_output,
        diagnostics_source = if (diagnostics) point_output else NULL
      ))
    }

    beta_draws <- extract_beta_draws(
      model = fit,
      names_ref = cb_names_ref,
      n = n_samples,
      coefficient_info = coefficient_info
    )
    draw_sd <- apply(beta_draws, 2, stats::sd)
    if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
      warning("Near-zero coefficient-draw variability for variable '",
              variable, "'. Intervals may collapse.", call. = FALSE)
    }

    eta_list <- lapply(seq_len(nrow(beta_draws)), function(i) {
      draw_prediction <- dlnm::crosspred(
        cb,
        coef = beta_draws[i, ],
        vcov = diag(0, ncol(beta_draws)),
        at = at_values,
        cen = center_value,
        bylag = 1
      )
      draw_prediction$matfit[, lag_order, drop = FALSE]
    })
    sample_output <- dplyr::bind_rows(lapply(seq_along(eta_list), function(i) {
      current <- build_output(eta_list[[i]])
      current$sample <- i
      current[, c("sample", setdiff(names(current), "sample")), drop = FALSE]
    }))

    if (output == "samples") {
      return(list(
        effects = sample_output,
        diagnostics_source = if (diagnostics) sample_output else NULL
      ))
    }

    grouping_columns <- intersect(
      c("var", "lag", "period", "scale", "value"),
      names(sample_output)
    )
    metric_columns <- intersect(
      c("eta", "effect", "delta", "predicted"),
      names(sample_output)
    )
    summary_output <- sample_output |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) |>
      dplyr::summarise(
        dplyr::across(
          dplyr::all_of(metric_columns),
          list(
            estimate = ~stats::median(.x, na.rm = TRUE),
            sd = ~safe_sd(.x),
            lower = ~safe_quantile(.x)[1],
            upper = ~safe_quantile(.x)[2]
          ),
          .names = "{.col}_{.fn}"
        ),
        baseline = dplyr::first(.data$baseline),
        .groups = "drop"
      )
    for (metric in metric_columns) {
      estimate_column <- paste0(metric, "_estimate")
      names(summary_output)[names(summary_output) == estimate_column] <- metric
    }

    list(
      effects = as.data.frame(summary_output),
      diagnostics_source = if (diagnostics) sample_output else NULL
    )
  }

  variable_results <- lapply(variables, summarise_one_variable)
  output_data <- dplyr::bind_rows(lapply(variable_results, `[[`, "effects"))
  rownames(output_data) <- NULL

  # Preserve the original return type and structure unless diagnostics are
  # explicitly requested.
  if (!diagnostics) return(output_data)

  diagnostics_source <- dplyr::bind_rows(
    lapply(variable_results, `[[`, "diagnostics_source")
  )
  diagnostics_data <- compute_lag_diagnostics(diagnostics_source)
  rownames(diagnostics_data) <- NULL

  structure(
    list(
      effects = output_data,
      diagnostics = diagnostics_data
    ),
    class = c("epiexposure_effects", "list")
  )
}
