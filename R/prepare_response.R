#' Prepare and validate the epidemic-level response for modeling
#'
#' Validates an epidemic-level response against the distribution family selected
#' by the user and creates the standardized `y_model` column required by
#' `fit_epidlnm()`.
#'
#' In the standard EpiExposure workflow:
#'
#' ```
#' templates <- define_exposures(...)
#' design    <- build_design(
#'   data = data,
#'   cb_templates = templates,
#'   include_response = TRUE
#' )
#' prepared  <- prepare_response(
#'   data = design,
#'   response = "y",
#'   family = "beta"
#' )
#' fit       <- fit_epidlnm(
#'   data = prepared,
#'   ...
#' )
#' ```
#'
#' `build_design(include_response = TRUE)` and `prepare_response()` have
#' different responsibilities. `build_design()` only carries the original
#' epidemic-level outcome into the design matrix. `prepare_response()` checks
#' whether that outcome is compatible with the requested statistical family and
#' performs only the explicitly documented scale conversion needed for Beta
#' percentages.
#'
#' `prepare_response()` does **not** choose a probability distribution from the
#' observed data and does not perform goodness-of-fit model selection. The
#' `family` is selected by the analyst and is then validated here.
#'
#' @param data Non-empty data frame containing the response variable. In the
#'   standard workflow this is the output of `build_design()`.
#' @param response Character scalar naming the original response variable in
#'   `data`. In the standard workflow this is `"y"`.
#'
#'   For compatibility with EpiExposure model-selection code, `response` may be
#'   omitted when `y_var` is supplied. If both are supplied they must identify
#'   the same column.
#' @param family Distribution family supplied as a supported character name or
#'   as a family-like object containing a valid `family` field.
#'
#'   Supported EpiExposure v1 canonical families are:
#'
#'   - `"beta"`: continuous proportions strictly inside `(0, 1)`;
#'   - `"binomial"`: binary outcomes coded exactly `0/1`;
#'   - `"poisson"`: non-negative integer counts;
#'   - `"negative_binomial"`: NB2 non-negative integer counts;
#'   - `"gaussian"`: finite continuous numeric values;
#'   - `"gamma"`: strictly positive continuous values.
#'
#'   Aliases such as `"bernoulli"`, `"normal"`, `"negbin"`, `"nbinom"`, and
#'   `"nbinom2"` are normalized to the corresponding canonical family.
#'
#'   `nbinom1` / NB1 are rejected explicitly. EpiExposure v1 standardizes
#'   `"negative_binomial"` to the NB2 variance parameterization. Ordinal
#'   outcomes are not supported in EpiExposure v1.
#' @param beta_scale Character scalar controlling interpretation of a Beta
#'   response. One of:
#'
#'   - `"auto"` (default): values already within `[0, 1]` are interpreted as
#'     proportions; otherwise, if all values lie within `[0, 100]` and at least
#'     one value exceeds 1, they are interpreted as percentages and divided by
#'     100;
#'   - `"proportion"`: require the supplied response to be on the `[0, 1]`
#'     scale;
#'   - `"percent"`: require values on the `[0, 100]` scale and divide by 100.
#'
#'   After scale handling, Beta regression still requires the **open** interval
#'   `(0, 1)`. Exact boundary values 0 or 1 are not altered silently.
#' @param y_var Optional compatibility alias for `response`. This is used by
#'   existing EpiExposure internal workflows. New user-facing code should prefer
#'   `response`.
#'
#' @return A data frame containing all original columns plus `y_model`.
#'
#'   Existing attributes on `data` (including the design/template metadata
#'   created by `build_design()`) are preserved. The following response metadata
#'   are added:
#'
#'   - `"response_family_name"`: canonical EpiExposure family;
#'   - `"response_name"`: original response-column name;
#'   - `"response_original_scale"`: response scale before preparation;
#'   - `"response_model_scale"`: scale represented by `y_model`;
#'   - `"response_transform"`: transformation applied, if any;
#'   - `"y_scale_mult"`: multiplier that maps the Beta model scale back to the
#'     original numerical scale (`100` after percent-to-proportion conversion,
#'     otherwise `1`);
#'   - `"response_nb_parameterization"`: `"NB2"` for negative-binomial
#'     responses, otherwise `NULL`;
#'   - `"response_contract"`: `"family_validated_response_v1"`.
#'
#' @details
#' ## Why this function remains necessary after `build_design()`
#'
#' With `include_response = TRUE`, `build_design()` ensures that one original
#' epidemic-level `y` value is carried into each design row. It deliberately
#' does not interpret that outcome statistically.
#'
#' `prepare_response()` is the family-aware layer between the design matrix and
#' `fit_epidlnm()`. `fit_epidlnm()` expects a standardized numeric `y_model`
#' column and checks that the response-family metadata agree with the family
#' used for fitting.
#'
#' ## Beta responses
#'
#' Beta regression requires
#'
#' \deqn{0 < y < 1.}
#'
#' When `beta_scale = "auto"`, the function retains the previous EpiExposure
#' convenience of converting a clear 0--100 percentage scale to proportions.
#' The conversion is recorded in output metadata.
#'
#' Earlier EpiExposure code silently replaced exact 0 and 1 by
#' `1e-5` and `1 - 1e-5`. That behavior is intentionally removed. Boundary
#' modification changes observed outcomes and the appropriate treatment depends
#' on the scientific data-generating process. Exact 0/1 values therefore
#' produce an explicit error rather than an undocumented numerical adjustment.
#'
#' ## Count responses
#'
#' Poisson and negative-binomial responses must be non-negative integers. Values
#' are validated but never truncated or rounded to make them valid.
#'
#' EpiExposure v1 uses NB2 whenever the canonical family is
#' `"negative_binomial"`. Explicit NB1 requests generate an error before model
#' fitting so different engines cannot silently use different
#' mean--variance relationships.
#'
#' ## Binomial, Gaussian, and Gamma responses
#'
#' Binomial responses must be numeric 0/1. Gaussian responses may be any finite
#' numeric values. Gamma responses must be strictly positive.
#'
#' ## No distribution inference
#'
#' This function validates a **chosen** family; it does not infer the best family
#' from histograms, moments, normality tests, overdispersion tests, or other
#' automatic rules. Distribution choice remains a modeling decision based on the
#' outcome definition and study design.
#'
#' @export
prepare_response <- function(
    data,
    response = NULL,
    family,
    beta_scale = c("auto", "proportion", "percent"),
    y_var = NULL
) {

  # BASIC VALIDATION

  if (!is.data.frame(data) || !nrow(data)) {
    stop(
      "`data` must be a non-empty data.frame.",
      call. = FALSE
    )
  }

  valid_name <- function(x) {
    is.character(x) &&
      length(x) == 1L &&
      !is.na(x) &&
      nzchar(x)
  }

  response_supplied <- !is.null(response)
  y_var_supplied <- !is.null(y_var)

  if (!response_supplied && !y_var_supplied) {
    stop(
      "Supply the response column through `response` (preferred) or `y_var`.",
      call. = FALSE
    )
  }

  if (response_supplied && !valid_name(response)) {
    stop(
      "`response` must be NULL or one non-empty column name.",
      call. = FALSE
    )
  }

  if (y_var_supplied && !valid_name(y_var)) {
    stop(
      "`y_var` must be NULL or one non-empty column name.",
      call. = FALSE
    )
  }

  if (response_supplied && y_var_supplied &&
      !identical(response, y_var)) {
    stop(
      "`response` and `y_var` were both supplied but identify different ",
      "columns. Supply only one name or make them identical.",
      call. = FALSE
    )
  }

  response_name <- if (response_supplied) {
    response
  } else {
    y_var
  }

  if (!response_name %in% names(data)) {
    stop(
      "Response column '",
      response_name,
      "' was not found in `data`.",
      call. = FALSE
    )
  }

  if (identical(response_name, "y_model")) {
    stop(
      "`response`/`y_var` must identify the original outcome, not the ",
      "prepared `y_model` column.",
      call. = FALSE
    )
  }

  if ("y_model" %in% names(data)) {
    stop(
      "`data` already contains `y_model`. `prepare_response()` should be ",
      "applied once to an unprepared response to avoid accidental ",
      "re-preparation or scale conversion.",
      call. = FALSE
    )
  }

  beta_scale <- match.arg(beta_scale)


  # STRICT EPIEXPOSURE v1 FAMILY RESOLUTION


  extract_family_label <- function(family_input) {
    if (is.character(family_input) &&
        length(family_input) == 1L &&
        !is.na(family_input) &&
        nzchar(family_input)) {
      return(family_input)
    }

    if (is.list(family_input) &&
        !is.null(family_input$family) &&
        length(family_input$family) >= 1L) {
      label <- as.character(family_input$family[[1L]])
      if (length(label) == 1L && !is.na(label) && nzchar(label)) {
        return(label)
      }
    }

    if (inherits(family_input, "family") &&
        !is.null(family_input$family)) {
      label <- as.character(family_input$family[[1L]])
      if (length(label) == 1L && !is.na(label) && nzchar(label)) {
        return(label)
      }
    }

    stop(
      "Unsupported `family` specification. Provide a supported EpiExposure v1 ",
      "family name or a family-like object containing a valid `family` field.",
      call. = FALSE
    )
  }

  normalize_family_label <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("[[:space:]-]+", "_", x)
    x <- gsub("[^a-z0-9_]", "", x)
    x
  }

  resolve_family_name_v1 <- function(family_input) {
    family_raw <- extract_family_label(family_input)
    normalized <- normalize_family_label(family_raw)

    # Explicit exclusions first so they cannot be swallowed by broader aliases.
    if (normalized %in% c(
      "ordinal",
      "cumulative",
      "ordered",
      "ordinal_logit"
    )) {
      stop(
        "Ordinal outcomes are not supported in EpiExposure v1.",
        call. = FALSE
      )
    }

    if (normalized %in% c(
      "nbinom1",
      "nb1",
      "negative_binomial_1",
      "negative_binomial1",
      "negbin1"
    )) {
      stop(
        "Negative-binomial NB1 is not supported in EpiExposure v1. ",
        "`negative_binomial` is standardized to NB2.",
        call. = FALSE
      )
    }

    if (normalized %in% c(
      "beta",
      "beta_family",
      "beta_proportion"
    )) {
      return("beta")
    }

    if (normalized %in% c(
      "binomial",
      "bernoulli"
    )) {
      return("binomial")
    }

    if (identical(normalized, "poisson")) {
      return("poisson")
    }

    if (identical(normalized, "gamma")) {
      return("gamma")
    }

    if (normalized %in% c(
      "gaussian",
      "normal"
    )) {
      return("gaussian")
    }

    if (normalized %in% c(
      "negative_binomial",
      "negativebinomial",
      "negative_binomial_2",
      "negative_binomial2",
      "negbin",
      "negbin2",
      "nbinom",
      "nbinom2"
    )) {
      return("negative_binomial")
    }

    # Family objects from some packages embed a dispersion value or other
    # suffix in the printed family label, e.g. "Negative Binomial(theta)".
    # Such generic Negative Binomial labels correspond to the package-wide
    # EpiExposure NB2 contract unless they explicitly identified NB1 above.
    if (grepl("^negative_binomial", normalized) ||
        grepl("^negativebinomial", normalized)) {
      return("negative_binomial")
    }

    stop(
      "Unsupported family: '",
      family_raw,
      "'. Supported EpiExposure v1 canonical families are: beta, binomial, ",
      "poisson, gamma, gaussian, and negative_binomial (NB2).",
      call. = FALSE
    )
  }

  family_name <- resolve_family_name_v1(family)

  supported_families <- c(
    "beta",
    "binomial",
    "poisson",
    "gamma",
    "gaussian",
    "negative_binomial"
  )

  if (!family_name %in% supported_families) {
    stop(
      "Unsupported EpiExposure v1 family: '",
      family_name,
      "'.",
      call. = FALSE
    )
  }

  # RESPONSE EXTRACTION

  original_y <- data[[response_name]]

  if (is.factor(original_y)) {
    stop(
      "The response variable cannot be a factor in EpiExposure v1. Supply the ",
      "numeric outcome required by the selected family.",
      call. = FALSE
    )
  }

  if (!is.numeric(original_y) &&
      !is.integer(original_y)) {
    stop(
      "The response variable must be numeric or integer.",
      call. = FALSE
    )
  }

  y <- as.numeric(original_y)

  if (!length(y)) {
    stop(
      "The response variable cannot be empty.",
      call. = FALSE
    )
  }

  if (anyNA(y) ||
      any(!is.finite(y))) {
    stop(
      "The response variable must contain only finite, non-missing values.",
      call. = FALSE
    )
  }

  original_range <- range(y)
  y_scale_mult <- 1
  original_scale <- "native"
  model_scale <- "native"
  response_transform <- "none"
  nb_parameterization <- NULL

  # BETA

  if (identical(family_name, "beta")) {
    if (any(y < 0)) {
      stop(
        "For `family = 'beta'`, response values cannot be negative.",
        call. = FALSE
      )
    }

    if (identical(beta_scale, "auto")) {
      if (all(y >= 0 & y <= 1)) {
        resolved_beta_scale <- "proportion"
      } else if (
        all(y >= 0 & y <= 100) &&
        any(y > 1)
      ) {
        resolved_beta_scale <- "percent"
      } else {
        stop(
          "`beta_scale = 'auto'` could not interpret the response. Beta data ",
          "must be proportions in [0, 1] or percentages in [0, 100]. ",
          "Specify `beta_scale` explicitly if needed.",
          call. = FALSE
        )
      }
    } else {
      resolved_beta_scale <- beta_scale
    }

    if (identical(resolved_beta_scale, "proportion")) {
      if (any(y < 0 | y > 1)) {
        stop(
          "`beta_scale = 'proportion'` requires all supplied response values ",
          "to lie in [0, 1].",
          call. = FALSE
        )
      }

      original_scale <- "proportion"
      model_scale <- "proportion"

    } else if (identical(resolved_beta_scale, "percent")) {
      if (any(y < 0 | y > 100)) {
        stop(
          "`beta_scale = 'percent'` requires all supplied response values to ",
          "lie in [0, 100].",
          call. = FALSE
        )
      }

      y <- y / 100
      y_scale_mult <- 100
      original_scale <- "percent"
      model_scale <- "proportion"
      response_transform <- "divide_by_100"
    }

    if (any(y <= 0 | y >= 1)) {
      boundary_values <- sort(unique(y[y <= 0 | y >= 1]))

      stop(
        "`family = 'beta'` requires every modeled response to lie strictly ",
        "inside (0, 1) after scale conversion. Boundary value(s) detected: ",
        paste(utils::head(boundary_values, 10L), collapse = ", "),
        if (length(boundary_values) > 10L) ", ..." else "",
        ". EpiExposure no longer moves 0/1 values inward with an arbitrary ",
        "epsilon. Handle boundary observations explicitly according to the ",
        "scientific outcome/model before fitting.",
        call. = FALSE
      )
    }
  }

  # BINOMIAL

  if (identical(family_name, "binomial")) {
    if (!all(y %in% c(0, 1))) {
      invalid_values <- sort(unique(y[!y %in% c(0, 1)]))

      stop(
        "For `family = 'binomial'`, every response must be coded exactly 0 or ",
        "1. Invalid value(s): ",
        paste(utils::head(invalid_values, 10L), collapse = ", "),
        if (length(invalid_values) > 10L) ", ..." else "",
        ".",
        call. = FALSE
      )
    }

    y <- as.integer(y)
    original_scale <- "binary_0_1"
    model_scale <- "binary_0_1"
  }

  # POISSON AND NEGATIVE BINOMIAL (NB2)

  if (family_name %in% c(
    "poisson",
    "negative_binomial"
  )) {
    if (any(y < 0)) {
      stop(
        "For `family = '",
        family_name,
        "'`, every response must be non-negative.",
        call. = FALSE
      )
    }

    tolerance <- sqrt(.Machine$double.eps)
    non_integer <- abs(y - round(y)) > tolerance

    if (any(non_integer)) {
      invalid_values <- sort(unique(y[non_integer]))

      stop(
        "For `family = '",
        family_name,
        "'`, every response must be an integer count. Invalid value(s): ",
        paste(utils::head(invalid_values, 10L), collapse = ", "),
        if (length(invalid_values) > 10L) ", ..." else "",
        ". Values are not rounded automatically.",
        call. = FALSE
      )
    }

    # Keep numeric storage rather than coercing to R integer, avoiding needless
    # integer-overflow risk for large valid counts.
    y <- round(y)
    original_scale <- "count"
    model_scale <- "count"

    if (identical(
      family_name,
      "negative_binomial"
    )) {
      nb_parameterization <- "NB2"
    }
  }

  # GAUSSIAN

  if (identical(family_name, "gaussian")) {
    y <- as.numeric(y)
    original_scale <- "continuous"
    model_scale <- "continuous"
  }

  # GAMMA

  if (identical(family_name, "gamma")) {
    if (any(y <= 0)) {
      invalid_values <- sort(unique(y[y <= 0]))

      stop(
        "For `family = 'gamma'`, every response must be strictly greater than ",
        "0. Invalid value(s): ",
        paste(utils::head(invalid_values, 10L), collapse = ", "),
        if (length(invalid_values) > 10L) ", ..." else "",
        ".",
        call. = FALSE
      )
    }

    y <- as.numeric(y)
    original_scale <- "positive_continuous"
    model_scale <- "positive_continuous"
  }

  # FINAL VALIDATION AND OUTPUT

  if (length(y) != nrow(data) ||
      anyNA(y) ||
      any(!is.finite(y))) {
    stop(
      "Internal response preparation failed to produce one finite `y_model` ",
      "value per input row.",
      call. = FALSE
    )
  }

  data$y_model <- y

  # Preserve all build_design() and other existing data-frame attributes; the
  # assignments below add only response-specific metadata.
  attr(data, "y_scale_mult") <- y_scale_mult
  attr(data, "response_family_name") <- family_name
  attr(data, "response_name") <- response_name
  attr(data, "response_original_scale") <- original_scale
  attr(data, "response_model_scale") <- model_scale
  attr(data, "response_transform") <- response_transform
  attr(data, "response_original_range") <- original_range
  attr(data, "response_model_range") <- range(y)
  attr(data, "response_nb_parameterization") <- nb_parameterization
  attr(data, "response_contract") <- "family_validated_response_v1"

  data
}
