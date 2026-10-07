#' Predict outcomes from fitted EpiExposure DLNM models
#'
#' Predicts the population-level expected outcome associated with either
#' user-defined chronological exposure profiles or a new longitudinal exposure
#' data set. The function uses only the standardized model metadata and internal
#' prediction helpers created by `fit_epidlnm()`; family, link, cross-basis
#' definitions, coefficient mappings, and engine behavior are never inferred
#' again inside this public function.
#'
#' @param fit Fitted model returned by the current `fit_epidlnm()`.
#' @param profiles Optional exposure profiles. Supply either `profiles` or
#'   `newdata`, but not both.
#'
#'   For a model containing one exposure, a single profile may be supplied as a
#'   numeric vector. For one or more exposures, use a named list whose names
#'   match the fitted exposures exactly. Each exposure element may be:
#'
#'   \itemize{
#'     \item one numeric vector, representing one chronological history;
#'     \item a list of numeric vectors, representing multiple histories;
#'     \item a numeric matrix or data frame, with one profile per row.
#'   }
#'
#'   Histories must be ordered from the oldest observation to the most recent
#'   observation. The most recent value corresponds to lag 0. Each history
#'   must contain exactly `max_lag+1` observations. Longer histories are not
#'   truncated and shorter histories are not padded.
#'   When several fitted exposures are supplied, profiles are
#'   matched by position rather than combined as a Cartesian product. An
#'   exposure with one profile is recycled across exposures that contain
#'   multiple profiles.
#' @param newdata Optional long-format data frame containing new observed
#'   exposure histories. Supply either `newdata` or `profiles`, but not both.
#'   The response variable is not required and, if present, is not used for
#'   prediction. All fitted exposure variables must be present. spaMM spatial
#'   coordinates and grouping factors are not required for fixed-component
#'   prediction; the same applies to `profiles`.
#' @param group Character scalar naming the column that identifies independent
#'   exposure histories in `newdata`, or `NULL` to treat all rows as one
#'   history. This argument is used only to split the longitudinal prediction
#'   data. It does not request group-specific random-effect or conditional
#'   spaMM spatial predictions. Coordinates and spatial-group variables are
#'   not required for population-level prediction.
#' @param time Character scalar naming the chronological numeric time column in
#'   `newdata`. Within each prediction group, times must be unique, complete,
#'   and equally spaced because one DLNM lag must represent the same time
#'   interval throughout the prediction data.
#' @param type Prediction scale: `"response"` for the expected outcome on its
#'   natural response scale, or `"link"` for the corresponding population-level
#'   linear predictor.
#' @param uncertainty Logical. If `FALSE`, return the deterministic prediction
#'   obtained from the central fixed/population parameter estimate. If `TRUE`,
#'   propagate fixed/population parameter uncertainty draw by draw.
#' @param output Character. When `uncertainty = TRUE`, `"summary"` returns the
#'   median, standard deviation, and empirical interval across expected-response
#'   draws, whereas `"samples"` returns every draw. When
#'   `uncertainty = FALSE`, only `"summary"` is valid and the returned data
#'   contain the deterministic `prediction`.
#' @param n_samples Positive integer number of parameter draws used when
#'   `uncertainty = TRUE`. At least two are required.
#' @param probs Numeric vector of length two defining the lower and upper
#'   empirical uncertainty quantiles. The default `c(0.025, 0.975)` gives a
#'   95 percent interval. Used only when `uncertainty = TRUE` and
#'   `output = "summary"`.
#' @param seed `NULL` or one strictly positive finite integer used to make
#'   parameter-draw sampling reproducible. For INLA-backed predictions, the
#'   value controls both the R random-number generator and the native
#'   `INLA::inla.posterior.sample()` seed. The caller's global random-number
#'   state is restored when the function exits.
#' @param extrapolation Behavior when prediction exposures extend beyond the
#'   exposure range stored with the fitted cross-basis: `"warn"` (default),
#'   `"error"`, or `"allow"`. This controls notification only; the fitted basis
#'   is never re-estimated from prediction data.
#'
#' @return A data frame.
#'
#'   For multiple explicit profiles, a `profile` column identifies matched
#'   profile combinations. For `newdata` with `group` supplied, the grouping
#'   column identifies each predicted history.
#'
#'   With `uncertainty = FALSE`, the result contains `prediction`.
#'
#'   With `uncertainty = TRUE` and `output = "summary"`, the result contains
#'   `prediction`, `prediction_sd`, `prediction_lower`, and
#'   `prediction_upper`.
#'
#'   With `uncertainty = TRUE` and `output = "samples"`, the result contains
#'   `sample` and `prediction`, with one row for each parameter draw and
#'   prediction history.
#'
#' @details
#' ## Prediction estimand
#'
#' EpiExposure v1 defines prediction as the population-level expected response.
#' Models may be fitted with conventional random intercepts and, for spaMM,
#' spatially autocorrelated Matérn random effects. All non-fixed effects are
#' excluded from predictions. This includes `(1 | epi_id)`, `(1 | block_id)`,
#' `Matern(1 | x_coord + y_coord)`, and independent field realizations such as
#' `Matern(1 | x_coord + y_coord %in% year)`; they are fitted but contribute
#' zero to the default prediction. In mixed-model notation, the prediction
#' target is based on
#'
#' \deqn{\eta = X\beta}
#'
#' rather than
#'
#' \deqn{\eta = X\beta + b_i.}
#'
#' Here, "population-level" means that conventional and spatial random
#' effects are set to zero or excluded from the prediction component. It does
#' not mean integration over their distributions. For nonlinear links, the
#' resulting expected response need not equal a marginally integrated mean.
#'
#' With `type = "response"`, the deterministic prediction is
#'
#' \deqn{\hat{\mu} = g^{-1}(X\hat{\beta}),}
#'
#' where \eqn{\hat{\beta}} is the harmonized central fixed/population parameter
#' estimate. Frequentist engines use their fitted fixed-effect estimates;
#' Bayesian engines use posterior means of the fixed/population coefficients.
#' In particular, deterministic Bayesian prediction is not defined as the
#' median of posterior expected-response draws.
#'
#' ## Uncertainty
#'
#' With `uncertainty = TRUE`, parameter uncertainty is propagated using
#'
#' \deqn{\mu^{(s)} = g^{-1}(X\beta^{(s)})}
#'
#' for each draw \eqn{s}. The inverse link is applied separately to every draw
#' before any summary is computed. The same coefficient draw is applied jointly
#' to all prediction rows, preserving covariance among predictions sharing the
#' same fitted parameters.
#'
#' For frequentist engines, fixed-effect draws are obtained from the estimated
#' joint coefficient covariance matrix. For Bayesian engines, posterior or
#' approximate-posterior fixed-effect draws are used. Group-specific
#' random-effect or spatial-field uncertainty is not propagated because
#' predictions depend exclusively on fixed/population parameters.
#'
#' The uncertainty distribution therefore describes uncertainty in the expected
#' response. It does not simulate a new observed outcome and does not add
#' residual, observation, process, dispersion, or posterior-predictive noise.
#'
#' ## New observed exposure data
#'
#' `newdata` is transformed using the cross-basis definition stored in `fit`.
#' Knots, boundary knots, spline definitions, and lag-basis parameters are not
#' re-estimated from the new observations. Each prediction group is sorted by
#' `time`, transformed with the fitted basis, and represented by the final
#' cross-basis row.
#'
#' This is the same train-to-test principle used for out-of-fold prediction in
#' `find_bestfit()`: the fitted/training basis defines the transformation and
#' the new or held-out exposure history is only projected through that stored
#' basis. spaMM models with Matérn terms do not require coordinates, spatial
#' groups, or levels of conventional random effects in prediction `newdata` or
#' `profiles`: these variables do not enter the fixed-component prediction.
#' For spatial spaMM models the internal point-prediction helper computes
#' `X %*% beta` directly; uncertainty draws likewise use only fixed-effect
#' coefficients, with no simulated spatial fields or new observations.
#'
#' ## Multiple explicit profiles
#'
#' Multiple profiles are paired by profile index. For example, if `tmean` and
#' `rain` each contain 100 profiles, profile 1 of `tmean` is combined with
#' profile 1 of `rain`, and so forth. The function does not create all
#' pairwise combinations. If one exposure contains a single profile and another
#' contains 100, the single profile is recycled 100 times.
#'
#' @export
predict_outcomes <- function(
    fit,
    profiles = NULL,
    newdata = NULL,
    group = "epi_id",
    time = "time",
    type = c("response", "link"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000L,
    probs = c(0.025, 0.975),
    seed = NULL,
    extrapolation = c("warn", "error", "allow")
) {

  # ARGUMENT VALIDATION AND STRICT MODEL CONTRACT

  type <- match.arg(type)
  output <- match.arg(output)
  extrapolation <- match.arg(extrapolation)

  if (is.null(fit)) {
    stop("`fit` cannot be NULL.", call. = FALSE)
  }

  metadata <- .get_epiexposure_metadata(fit)

  has_profiles <- !is.null(profiles)
  has_newdata <- !is.null(newdata)

  if (identical(has_profiles, has_newdata)) {
    stop(
      "Supply exactly one of `profiles` or `newdata`.",
      call. = FALSE
    )
  }

  if (!is.logical(uncertainty) || length(uncertainty) != 1L ||
      is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!uncertainty && identical(output, "samples")) {
    stop(
      "`output = 'samples'` requires `uncertainty = TRUE`. ",
      "Deterministic prediction contains one central expected-response estimate ",
      "and therefore has no parameter-draw sample index.",
      call. = FALSE
    )
  }

  if (uncertainty) {
    n_samples <- .epix_validate_n_samples(n_samples)

    if (n_samples < 2L) {
      stop(
        "`n_samples` must be at least 2 when `uncertainty = TRUE`.",
        call. = FALSE
      )
    }

    if (identical(output, "summary")) {
      probs <- .epix_validate_probs(probs)
    }
  }

  seed <- .epix_validate_seed(
    seed
  )

  # The current fitted-model metadata must encode exactly the prediction
  # contract agreed for EpiExposure v1. .get_epiexposure_metadata() already
  # stops if these attributes are absent or inconsistent; these checks are kept
  # here as explicit public-function invariants.
  if (!identical(metadata$prediction_level, "population")) {
    stop(
      "`predict_outcomes()` supports population-level prediction only.",
      call. = FALSE
    )
  }

  if (!identical(metadata$prediction_estimand, "expected_response")) {
    stop(
      "`predict_outcomes()` requires the fitted prediction estimand to be ",
      "'expected_response'.",
      call. = FALSE
    )
  }

  # LOCAL RNG HANDLING

  # A user-supplied seed makes parameter sampling reproducible without leaving
  # a changed .Random.seed in the caller's global environment.
  if (!is.null(seed)) {
    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    if (had_random_seed) {
      old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }

    on.exit(
      {
        if (had_random_seed) {
          assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    set.seed(seed)
  }

  # OUTPUT HELPERS

  attach_contract_attributes <- function(
    result,
    source,
    prediction_method,
    n_predictions
  ) {
    attr(result, "epiexposure_prediction_source") <- source
    attr(result, "epiexposure_prediction_level") <- "population"
    attr(result, "epiexposure_prediction_estimand") <- "expected_response"
    attr(result, "epiexposure_prediction_type") <- type
    attr(result, "epiexposure_uncertainty") <- uncertainty
    attr(result, "epiexposure_prediction_method") <- prediction_method
    attr(result, "epiexposure_n_predictions") <- as.integer(n_predictions)

    if (uncertainty) {
      attr(result, "epiexposure_n_samples") <- as.integer(n_samples)
      if (identical(output, "summary")) {
        attr(result, "epiexposure_probs") <- probs
        attr(result, "epiexposure_summary_center") <- "median"
      }
    }

    result
  }

  format_deterministic <- function(prediction, keys, source) {
    prediction_method <- attr(
      prediction,
      "epiexposure_prediction_method",
      exact = TRUE
    )

    if (is.null(prediction_method)) {
      prediction_method <- "central_expected_response"
    }

    prediction <- as.numeric(prediction)

    if (length(prediction) != nrow(keys) && nrow(keys) > 0L) {
      stop(
        "Prediction/key length mismatch while formatting deterministic output.",
        call. = FALSE
      )
    }

    values <- data.frame(
      prediction = prediction,
      stringsAsFactors = FALSE
    )

    result <- if (ncol(keys)) {
      cbind(keys, values)
    } else {
      values
    }

    rownames(result) <- NULL

    attach_contract_attributes(
      result = result,
      source = source,
      prediction_method = prediction_method,
      n_predictions = length(prediction)
    )
  }

  format_summary <- function(draws, keys, source) {
    summary_values <- .summarise_prediction_draws(
      draws = draws,
      probs = probs,
      value_name = "prediction"
    )

    if (nrow(summary_values) != nrow(keys) && nrow(keys) > 0L) {
      stop(
        "Prediction/key length mismatch while formatting uncertainty summary.",
        call. = FALSE
      )
    }

    summary_values$row_id <- NULL

    result <- if (ncol(keys)) {
      cbind(keys, summary_values)
    } else {
      summary_values
    }

    rownames(result) <- NULL

    attach_contract_attributes(
      result = result,
      source = source,
      prediction_method = "parameter_draws_xbeta",
      n_predictions = nrow(summary_values)
    )
  }

  format_samples <- function(draws, keys, source) {
    draws <- as.matrix(draws)

    if (ncol(draws) != nrow(keys) && nrow(keys) > 0L) {
      stop(
        "Prediction/key length mismatch while formatting parameter draws.",
        call. = FALSE
      )
    }

    blocks <- vector("list", ncol(draws))

    for (j in seq_len(ncol(draws))) {
      key_block <- if (ncol(keys)) {
        keys[rep(j, nrow(draws)), , drop = FALSE]
      } else {
        data.frame(row.names = seq_len(nrow(draws)))
      }

      block <- cbind(
        key_block,
        data.frame(
          sample = seq_len(nrow(draws)),
          prediction = as.numeric(draws[, j]),
          stringsAsFactors = FALSE
        )
      )

      blocks[[j]] <- block
    }

    result <- do.call(rbind, blocks)
    rownames(result) <- NULL

    attach_contract_attributes(
      result = result,
      source = source,
      prediction_method = "parameter_draws_xbeta",
      n_predictions = ncol(draws)
    )
  }

  # PROFILE INPUT NORMALIZATION


  as_profile_list <- function(x, variable) {

    # One numeric vector = one chronological profile.
    if (is.numeric(x) && is.null(dim(x))) {
      return(list(as.numeric(x)))
    }

    # Matrix/data.frame = one profile per row.
    if (is.matrix(x) || is.data.frame(x)) {
      x_matrix <- as.matrix(x)

      if (!is.numeric(x_matrix)) {
        stop(
          "Matrix/data.frame profiles for variable '", variable,
          "' must be numeric.",
          call. = FALSE
        )
      }

      return(
        lapply(
          seq_len(nrow(x_matrix)),
          function(i) as.numeric(x_matrix[i, , drop = TRUE])
        )
      )
    }

    # List = multiple explicitly supplied profiles.
    if (is.list(x)) {
      if (!length(x)) {
        stop(
          "No profiles were supplied for variable '", variable, "'.",
          call. = FALSE
        )
      }

      return(
        lapply(
          seq_along(x),
          function(i) {
            current <- x[[i]]

            if (is.data.frame(current) || is.matrix(current)) {
              current_matrix <- as.matrix(current)

              if (nrow(current_matrix) != 1L || !is.numeric(current_matrix)) {
                stop(
                  "Nested profile ", i, " for variable '", variable,
                  "' must contain exactly one numeric row.",
                  call. = FALSE
                )
              }

              current <- as.numeric(current_matrix[1, , drop = TRUE])
            }

            if (!is.numeric(current) || is.list(current) ||
                !is.null(dim(current))) {
              stop(
                "Profile ", i, " for variable '", variable,
                "' must be a numeric vector.",
                call. = FALSE
              )
            }

            as.numeric(current)
          }
        )
      )
    }

    stop(
      "Profiles for variable '", variable,
      "' must be a numeric vector, list, matrix, or data.frame.",
      call. = FALSE
    )
  }

  normalize_profiles <- function(input_profiles) {

    # Preserve the convenient single-exposure shorthand from the previous API.
    if (is.numeric(input_profiles) && is.null(dim(input_profiles))) {
      if (length(metadata$vars) != 1L) {
        stop(
          "A numeric `profiles` vector can only be used when the fitted model ",
          "contains one exposure. For multiple exposures, supply a named list.",
          call. = FALSE
        )
      }

      input_profiles <- stats::setNames(
        list(as.numeric(input_profiles)),
        metadata$vars
      )
    }

    if (!is.list(input_profiles) || is.null(names(input_profiles)) ||
        !length(input_profiles) || anyNA(names(input_profiles)) ||
        any(!nzchar(names(input_profiles))) ||
        anyDuplicated(names(input_profiles))) {
      stop(
        "`profiles` must be a numeric vector for a single-exposure model or ",
        "a named list with unique exposure names.",
        call. = FALSE
      )
    }

    missing_vars <- setdiff(metadata$vars, names(input_profiles))
    extra_vars <- setdiff(names(input_profiles), metadata$vars)

    if (length(missing_vars) || length(extra_vars)) {
      details <- c(
        if (length(missing_vars)) {
          paste0("missing: ", paste(missing_vars, collapse = ", "))
        },
        if (length(extra_vars)) {
          paste0("unexpected: ", paste(extra_vars, collapse = ", "))
        }
      )

      stop(
        "`profiles` names must match the fitted exposures exactly (",
        paste(details, collapse = "; "), ").",
        call. = FALSE
      )
    }

    # Reorder to the official fitted exposure order before any recycling.
    input_profiles <- input_profiles[metadata$vars]

    profile_sets <- lapply(
      metadata$vars,
      function(variable) as_profile_list(input_profiles[[variable]], variable)
    )
    names(profile_sets) <- metadata$vars

    profile_counts <- vapply(profile_sets, length, integer(1))
    n_profiles <- max(profile_counts)

    incompatible <- profile_counts != 1L & profile_counts != n_profiles

    if (any(incompatible)) {
      stop(
        "Every fitted exposure must provide either one profile or the common ",
        "maximum number of profiles. Received: ",
        paste(
          paste0(names(profile_counts), "=", profile_counts),
          collapse = ", "
        ),
        ".",
        call. = FALSE
      )
    }

    for (variable in metadata$vars) {
      if (length(profile_sets[[variable]]) == 1L && n_profiles > 1L) {
        profile_sets[[variable]] <- rep(
          profile_sets[[variable]],
          n_profiles
        )
      }
    }

    # Validate finite values before basis construction so malformed input fails
    # with a profile-specific message rather than inside dlnm.
    for (variable in metadata$vars) {
      expected_length <- metadata$spec[[variable]]$max_lag + 1L

      for (i in seq_len(n_profiles)) {
        current <- profile_sets[[variable]][[i]]

        if (!length(current) ||
            anyNA(current) ||
            any(!is.finite(current))) {
          stop(
            "Profile ",
            i,
            " for variable '",
            variable,
            "' must contain only finite numeric values.",
            call. = FALSE
          )
        }

        if (length(current) != expected_length) {
          stop(
            "Profile ",
            i,
            " for variable '",
            variable,
            "' contains ",
            length(current),
            " value(s), but exactly ",
            expected_length,
            " are required for max_lag = ",
            metadata$spec[[variable]]$max_lag,
            ". Histories are not truncated or padded.",
            call. = FALSE
          )
        }
      }
    }

    list(
      sets = profile_sets,
      n_profiles = n_profiles
    )
  }

  # BUILD CANONICAL PREDICTION DESIGN

  if (has_profiles) {

    normalized <- normalize_profiles(profiles)
    n_predictions <- normalized$n_profiles

    design_blocks <- vector("list", n_predictions)

    for (i in seq_len(n_predictions)) {
      current_profiles <- lapply(
        metadata$vars,
        function(variable) normalized$sets[[variable]][[i]]
      )
      names(current_profiles) <- metadata$vars

      design_blocks[[i]] <- .build_newdata_basis(
        fit = fit,
        profiles = current_profiles,
        extrapolation = extrapolation
      )
    }

    prediction_design <- do.call(rbind, design_blocks)
    rownames(prediction_design) <- NULL

    keys <- if (n_predictions > 1L) {
      data.frame(
        profile = seq_len(n_predictions),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(row.names = 1L)
    }

    source <- "profiles"

  } else {

    # .build_newdata_basis() performs the strict checks for fitted exposures,
    # temporal coverage, duplicated/irregular time points, stored basis reuse,
    # basis dimensionality, and extrapolation.
    prediction_design <- .build_newdata_basis(
      fit = fit,
      newdata = newdata,
      group = group,
      time = time,
      extrapolation = extrapolation
    )

    n_predictions <- nrow(prediction_design)

    if (is.null(group)) {
      keys <- data.frame(row.names = seq_len(n_predictions))
    } else {
      # Preserve the original type and first-occurrence order of the user's
      # grouping variable. The internal basis helper uses the same order.
      group_character <- as.character(newdata[[group]])
      first_rows <- which(!duplicated(group_character))
      keys <- newdata[first_rows, group, drop = FALSE]

      if (nrow(keys) != n_predictions) {
        stop(
          "Internal grouping-key mismatch while preparing prediction output.",
          call. = FALSE
        )
      }

      rownames(keys) <- NULL
    }

    source <- "newdata"
  }

  # DETERMINISTIC CENTRAL EXPECTED-RESPONSE PREDICTION

  if (!uncertainty) {
    point_prediction <- .predict_point_population(
      fit = fit,
      newdata = prediction_design,
      type = type,
      # No silent fallback: if a supported frequentist native predictor cannot
      # be used, the user is informed before the equivalent X-beta route is used.
      warn_fallback = TRUE
    )

    return(
      format_deterministic(
        prediction = point_prediction,
        keys = keys,
        source = source
      )
    )
  }

  # DRAW-BY-DRAW PARAMETER UNCERTAINTY

  prediction_draws <- .predict_draws_population(
    fit = fit,
    newdata = prediction_design,
    n_samples = n_samples,
    type = type,
    seed = seed
  )

  if (identical(output, "samples")) {
    return(
      format_samples(
        draws = prediction_draws,
        keys = keys,
        source = source
      )
    )
  }

  format_summary(
    draws = prediction_draws,
    keys = keys,
    source = source
  )
}
