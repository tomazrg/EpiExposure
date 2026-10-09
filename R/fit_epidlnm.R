#' Fit a harmonized DLNM inferential model
#'
#' Fits a distributed lag nonlinear model (DLNM) from an epidemic-level design
#' matrix and stores the standardized EpiExposure metadata required by
#' downstream prediction, effect summarization, scenario simulation, and model
#' comparison functions.
#'
#' @param data A data.frame containing `y_model` and the fitted DLNM design
#'   columns. For all engines except `bdlnm`, cross-basis columns must follow the
#'   package convention `cb_<variable>_<index>`. For `bdlnm`, epidemic-level,
#'   matrix-form cross-basis objects are supplied through `basis_objects`
#'   instead. Data prepared by `prepare_response()` may carry the
#'   `response_family_name` attribute; when present, it must agree with `family`.
#'
#' @param model_engine Character scalar identifying the modeling engine. One of
#'   `"glm"`, `"glmmTMB"`, `"gam"`, `"gamm"`, `"gls"`, `"spamm"`,
#'   `"brms"`, `"inla"`, or `"bdlnm"`.
#' @param family Response distribution supplied as a supported canonical name
#'   or as a family object from which a supported canonical family and link can
#'   be identified. The supported EpiExposure v1 families are `"beta"`,
#'   `"binomial"`, `"poisson"`, `"gamma"`, `"gaussian"`, and
#'   `"negative_binomial"`. The canonical negative-binomial model is the
#'   quadratic-variance NB2 parameterization. NB1 and ordinal models are not
#'   supported in this version and produce explicit errors.
#' @param random_effect Optional character scalar naming one grouping column.
#'   In EpiExposure v1, this argument represents a conventional random
#'   intercept only. Random slopes, nested or crossed random-effect
#'   specifications, and arbitrary engine-specific random-effect expressions
#'   are not supported.
#'
#'   For the Bayesian `brms` engine, the random intercept is fitted as
#'   `(1 | group)`. For the INLA-backed `inla` and `bdlnm` engines, it is fitted
#'   as `INLA::f(group, model = "iid")`. EpiExposure v1 does not expose other
#'   INLA latent structures through this argument, including `"rw1"`, `"rw2"`,
#'   `"ar1"`, `"besag"`, `"bym"`, `"bym2"`, or SPDE-based effects.
#'
#'   INLA-backed fits internally re-index arbitrary grouping labels to
#'   consecutive integers while retaining the original column name in
#'   EpiExposure metadata. `glm` and `gls` do not support this argument. `gamm`
#'   requires it in the current EpiExposure interface because its generated
#'   fixed formula contains no other smooth or random term.
#'
#' @param random_effect_prior Optional named list defining the hyperprior for
#'   the precision of the random intercept when `model_engine = "inla"` or
#'   `"bdlnm"`. The list is passed unchanged to the `hyper` argument of the
#'   internally generated `INLA::f(..., model = "iid")` term.
#'
#'   For example,
#'   `list(prec = list(prior = "pc.prec", param = c(1, 0.01)))`.
#'   `NULL` retains the engine's default prior. A non-`NULL` value requires
#'   `random_effect` and is currently supported only for the INLA-backed
#'   engines.
#'
#'   This argument customizes the hyperprior of the supported IID random
#'   intercept; it does not select or configure other INLA latent models.
#'   Structures such as `"rw1"`, `"rw2"`, `"ar1"`, `"besag"`, `"bym"`,
#'   `"bym2"`, and SPDE effects are not supported by the EpiExposure v1
#'   random-effect interface.
#'
#' @param spatial_effect NULL (default) or two distinct names of numeric,
#'   finite coordinate columns in the epidemic-level design. Only
#'   `model_engine = "spamm"` supports spatial autocorrelation. Coordinates may
#'   repeat across different epidemics; at least two distinct pairs are needed.
#' @param spatial_structure Spatial correlation structure; currently only
#'   `"matern"` (case-insensitive) when `spatial_effect` is supplied. This
#'   argument does not change non-spatial fits.
#' @param spatial_group NULL for one shared spatial field, or the name of one
#'   grouping column (e.g., `"year"`) for independent Matérn realizations with
#'   shared correlation parameters. Requires `spatial_effect`; the source column
#'   may be factor, character, or integer and is converted in a local fit copy.
#'   `random_effect` remains a separate, conventional random intercept and both
#'   effects can be fitted together. Coordinates and spatial groups are not used
#'   as epidemic-profile identifiers.
#' @param epiexposure_spec Named list describing the cross-basis construction
#'   for every fitted exposure. This metadata is required so downstream
#'   functions reconstruct exactly the fitted exposure-lag basis rather than
#'   re-inferring `df`, knots, or basis functions. Each element must contain
#'   `max_lag`, `argvar`, and `arglag`.
#'
#'   Under the EpiExposure exact-profile contract, each `max_lag` must be one
#'   non-negative integer and all fitted exposures must use the same
#'   `max_lag`. Consequently, every original exposure profile used to build the
#'   epidemic-level design must contain exactly `max_lag + 1` observations for
#'   every fitted exposure. `fit_epidlnm()` receives the already collapsed
#'   epidemic-level design and therefore validates the common fitted `max_lag`;
#'   exact original profile length is validated upstream by the exposure/design
#'   construction functions.
#' @param basis_objects Optional named list of `dlnm::crossbasis()` objects.
#'   For `model_engine = "bdlnm"`, these must be the epidemic-level,
#'   matrix-form cross-basis objects containing exactly one row per model row.
#'   They are normally recovered from the
#'   `epiexposure_bdlnm_basis_objects` attribute produced by `build_design()`,
#'   but may also be supplied explicitly. For other engines, when supplied,
#'   their names must match the fitted exposure variables.
#'
#' @param ... Named additional arguments passed to the selected engine. Core
#'   arguments managed by EpiExposure (`formula`/`model`, `data`, `family`, and
#'   the engine-specific random-effect argument) cannot be supplied again in
#'   `...`. For INLA-backed fits (`inla` and `bdlnm`),
#'   `control.compute$config = TRUE` is required by the EpiExposure uncertainty
#'   contract and is added automatically when absent. If a conflicting value
#'   is supplied, the function stops.
#'
#' @return A fitted model object with a strict metadata contract. The following
#'   attributes are attached for downstream functions:
#'   `epiexposure_family_name`, `epiexposure_link`, `epiexposure_engine`,
#'   `epiexposure_cb_cols`, `epiexposure_vars`,
#'   `epiexposure_data_template`, `epiexposure_id_col`, `epiexposure_spec`, and
#'   `epiexposure_basis_objects`. Additional attributes record the standardized
#'   family parameterization, random-intercept structure, the INLA latent
#'   random-effect model in `epiexposure_random_effect_model`, the optional
#'   `epiexposure_random_effect_prior`, and the spatial specification
#'   (`epiexposure_spatial_effect`, `epiexposure_spatial_structure`,
#'   `epiexposure_spatial_group`, and `epiexposure_spatial_term`), common fitted
#'   `max_lag`, expected profile length (`max_lag + 1`), the exact-profile
#'   contract, the validated temporal step in
#'   `epiexposure_time_step`, and the EpiExposure v1 prediction contract.
#'
#' @details
#' ## Engine harmonization
#'
#' `fit_epidlnm()` standardizes the *statistical target* across engines without
#' requiring every engine to use the same computational implementation. The
#' model is fitted using the selected engine, while downstream EpiExposure v1
#' predictions are defined as population-level expected responses. Random
#' effects may therefore affect estimation during model fitting, but downstream
#' predictions exclude group-specific random-effect contributions.
#'
#' The random-intercept translation used at fitting is:
#'
#' * `glmmTMB`, `brms`, and `spaMM`: `(1 | group)`;
#' * `gam`: `s(group, bs = "re")`;
#' * `gamm`: `random = list(group = ~1)`;
#' * `INLA` and `bdlnm`: `f(group, model = "iid")`, optionally extended to
#'   `f(group, model = "iid", hyper = random_effect_prior)` when a custom
#'   random-effect hyperprior is supplied.
#'
#' ## Random-effect scope in EpiExposure v1
#'
#' The EpiExposure v1 random-effect interface is intentionally restricted to
#' one conventional random intercept. For `brms`, this corresponds to
#' `(1 | group)`; for `inla` and `bdlnm`, it corresponds to
#' `f(group, model = "iid")`.
#'
#' The interface does not currently construct structured Bayesian latent
#' effects such as first- or second-order random walks (`"rw1"` or `"rw2"`),
#' autoregressive effects (`"ar1"`), areal spatial effects (`"besag"`, `"bym"`,
#' or `"bym2"`), or continuous spatial effects based on SPDE models. These
#' structures may be available in the underlying INLA framework but are outside
#' the harmonized fitting, prediction, validation, and uncertainty contract of
#' EpiExposure v1.
#'
#' The `spatial_effect` interface described below is separate from these
#' Bayesian latent structures and currently supports Matérn spatial covariance
#' through `spaMM` only.
#'
#' `glm` and `gls` are fixed-effect engines in this interface and reject a
#' non-`NULL` `random_effect`. The `mgcv::gamm()` implementation uses the
#' engine's native `random` argument rather than lme4-style syntax. Because
#' `gamm()` uses PQL for non-Gaussian responses and is specifically known to be
#' problematic for binary data, a warning is issued for the binomial family.
#'
#' ## Spatial covariance (spaMM only)
#'
#' Spatial terms are explicitly distinct from `random_effect`: the former
#' model correlated Matérn fields, whereas the latter is an ordinary intercept.
#' `spatial_effect = c("x_coord", "y_coord")` adds
#' `Matern(1 | x_coord + y_coord)`; adding `spatial_group = "year"` instead
#' uses `Matern(1 | x_coord + y_coord %in% year)`. `spatial_effect = NULL`
#' leaves every pre-existing engine and random-intercept fit unchanged.
#' Coordinates and spatial grouping affect estimation, but not EpiExposure's
#' default fixed-component prediction, which sets all non-fixed contributions
#' to zero rather than integrating them over their distributions.
#'
#' ## Families and links
#'
#' Canonical families are translated to engine-native representations. The
#' canonical `negative_binomial` family is standardized to an NB2-type model
#' (`Var(Y) = mu + mu^2 / shape`) where supported. An explicit NB1 request is
#' rejected instead of being silently converted.
#'
#' Automatic family availability is engine-specific:
#'
#' * `glm`: Gaussian, Binomial, Poisson, Gamma;
#' * `glmmTMB`: all six EpiExposure v1 families;
#' * `gam`: all six EpiExposure v1 families;
#' * `gamm`: Gaussian, Binomial, Poisson, Gamma;
#' * `gls`: Gaussian with identity link only;
#' * `spaMM`: all six EpiExposure v1 families;
#' * `brms`: all six EpiExposure v1 families;
#' * `INLA`: all six EpiExposure v1 families, subject to the installed INLA
#'   likelihood implementation;
#' * `bdlnm`: all six EpiExposure v1 families through its INLA backend.
#'
#' When a family object supplies a non-default link, EpiExposure records that
#' exact link in `epiexposure_link` and reconstructs the corresponding
#' engine-native family. For INLA-based engines the same link is forwarded via
#' `control.family$control.link$model`. Unsupported family-link combinations
#' fail explicitly in EpiExposure or in the selected engine; no Gaussian or
#' identity-link fallback is used.
#'
#' To specify a non-default link, supply a supported family object,
#' such as family = stats::Gamma(link = "inverse");
#' character family names use the EpiExposure default link.
#'
#' For `gam`, Beta and negative-binomial models use `mgcv` extended families.
#' When needed, `method = "REML"` is supplied by default to satisfy the
#' supported fitting route; an explicitly supplied incompatible method is
#' rejected. A GAM random intercept also defaults to REML unless the user
#' explicitly selects another method.
#'
#' For `model_engine = "gamm"`, the returned object contains the native
#' `gam` and `lme` components and inherits from `"epiexposure_gamm"`.
#' Calling `summary()` returns the fixed/population GAM summary by default;
#' use `component = "lme"` or `"both"` to inspect the mixed-model component.
#'
#' ## Response-scale validation
#'
#' To avoid engine-specific silent coercion, `y_model` is checked against the
#' canonical family before fitting: Beta responses must lie strictly inside
#' `(0, 1)`; Binomial responses must be coded `0/1`; Poisson and negative-
#' binomial responses must be non-negative integer counts; Gamma responses must
#' be strictly positive; Gaussian responses need only be finite numeric values.
#'
#' ## Common lag and exact-profile contract
#'
#' EpiExposure models use one common retrospective lag window for every fitted
#' exposure. If the common maximum lag is `L`, all original exposure profiles
#' must contain exactly
#'
#' \deqn{
#'   L + 1
#' }
#'
#' equally spaced observations per epidemiological unit and per exposure.
#' Exposures with different fitted `max_lag` values are rejected here.
#'
#' Because `fit_epidlnm()` is called after `build_design()` has collapsed each
#' complete profile to one epidemic-level model row, the original long-format
#' row count is no longer available at this stage. Exact `L + 1` temporal
#' coverage must therefore be enforced by `define_exposures()` and
#' `build_design()`. This function records the common lag/profile contract in
#' the fitted model metadata so downstream functions can enforce it.
#'
#' ## Metadata and downstream prediction contract
#'
#' The fitted exposure specification is part of the model contract, not an
#' optional hint. Downstream functions must use the stored family, link, engine,
#' cross-basis mapping, exposure specification, and data template. They should
#' not infer a link from the family, rebuild a cross-basis from user-supplied
#' degrees of freedom, align coefficients by position, or substitute a
#' Gaussian/identity model when metadata are missing.
#'
#' EpiExposure v1 prediction functions are intended to return the expected
#' outcome under the fitted model, not a newly simulated observation. They use
#' population-level predictions: random effects are included in model fitting
#' when requested here, but group-specific random effects are not included in
#' downstream predictions. This function itself only fits the model and does
#' not create outcome predictions or uncertainty summaries.
#'
#' For `glm` and `glmmTMB`, models are fitted with the fully evaluated
#' EpiExposure design, but their stored calls are compacted after fitting.
#' The original `data` expression is retained for display instead of the
#' complete evaluated data frame and its EpiExposure attributes. For `glm`,
#' the stored function and family expressions are also compacted to prevent
#' `summary()` from printing the evaluated function definition and complete
#' family object. These changes affect only the stored calls used for display
#' and re-evaluation; they do not alter fitted coefficients, likelihoods,
#' covariance matrices, predictions, or EpiExposure metadata.
#'
#' Bayesian engines may generate posterior coefficient draws as part of their
#' native fitting procedure (notably `bdlnm`). Those stored coefficient draws
#' are model output and are available for later draw-by-draw uncertainty
#' propagation; they are not used here as a single deterministic prediction.
#'
#' Ordinal outcomes are deliberately excluded from EpiExposure v1 so fitting,
#' prediction, uncertainty, performance metrics, and ensemble behavior remain
#' harmonized across the supported response families.
#'
#' @examples
#' \dontrun{
#' # These alternatives assume an epidemic-level design `dat` and
#' # a matching cross-basis specification `spec` already exist:
#' fit_epidlnm(dat, "spamm", "poisson", random_effect = "epi_id",
#'             spatial_effect = NULL, epiexposure_spec = spec)
#' fit_epidlnm(dat, "spamm", "poisson", spatial_effect = c("x_coord", "y_coord"),
#'             epiexposure_spec = spec)
#' fit_epidlnm(dat, "spamm", "poisson", random_effect = "block_id",
#'             spatial_effect = c("x_coord", "y_coord"), epiexposure_spec = spec)
#' fit_epidlnm(dat, "spamm", "poisson", random_effect = "block_id",
#'             spatial_effect = c("x_coord", "y_coord"),
#'             spatial_group = "year", epiexposure_spec = spec)
#' }
#' @export
fit_epidlnm <- function(
    data,
    model_engine,
    family,
    random_effect = NULL,
    random_effect_prior = NULL,
    epiexposure_spec = NULL,
    basis_objects = NULL,
    spatial_effect = NULL,
    spatial_structure = "matern",
    spatial_group = NULL,
    ...
) {

  # SMALL INTERNAL HELPERS

  `%||%` <- function(a, b) if (!is.null(a)) a else b


  stopf <- function(...) stop(..., call. = FALSE)


  is_scalar_string <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }


  quote_name <- function(x) {
    if (!is_scalar_string(x)) stopf("Internal error: invalid name to quote.")
    if (grepl("`", x, fixed = TRUE)) {
      stopf("Column and basis names containing backticks are not supported.")
    }
    paste0("`", x, "`")
  }


  normalize_token <- function(x) {
    x <- tolower(trimws(as.character(x)[1]))
    x <- gsub("[[:space:]-]+", "_", x)
    gsub("[^a-z0-9_]", "", x)
  }


  extract_family_raw <- function(family_input) {
    if (is_scalar_string(family_input)) return(family_input)


    if (is.list(family_input) && !is.null(family_input$family) &&
        length(family_input$family) >= 1L) {
      return(as.character(family_input$family[[1]]))
    }


    stopf(
      "Unsupported `family` specification. Supply one supported canonical ",
      "family name or a family object containing a `family` field."
    )
  }


  resolve_family_info <- function(family_input) {
    raw <- extract_family_raw(family_input)
    z <- normalize_token(raw)


    if (z %in% c("ordinal", "cumulative")) {
      return(list(name = "ordinal", variant = "ordinal", raw = raw))
    }


    if (z %in% c(
      "nb1", "nbinom1", "negative_binomial_1", "negativebinomial1",
      "negative_binomial_type_1", "negativebinomialtype1"
    )) {
      return(list(name = "negative_binomial", variant = "NB1", raw = raw))
    }


    if (z %in% c(
      "nb2", "negbin", "nbinom", "nbinom2", "negative_binomial",
      "negative_binomial_2", "negativebinomial", "negativebinomial2",
      "negative_binomial_type_2", "negativebinomialtype2"
    )) {
      return(list(name = "negative_binomial", variant = "NB2", raw = raw))
    }


    if (grepl("negative_?binomial", z)) {
      stopf(
        "Unknown negative-binomial parameterization: '", raw, "'. ",
        "EpiExposure v1 accepts only the explicitly recognized NB2 variants."
      )
    }


    if (z %in% c("beta", "beta_family", "beta_proportion", "beta_regression", "betar", "beta_resp")) {
      return(list(name = "beta", variant = "mean_precision", raw = raw))
    }


    if (z %in% c("binomial", "bernoulli")) {
      return(list(name = "binomial", variant = "bernoulli", raw = raw))
    }


    if (z == "poisson") {
      return(list(name = "poisson", variant = "poisson", raw = raw))
    }


    if (z == "gamma") {
      return(list(name = "gamma", variant = "gamma", raw = raw))
    }


    if (z %in% c("gaussian", "normal")) {
      return(list(name = "gaussian", variant = "gaussian", raw = raw))
    }


    stopf(
      "Unsupported family: '", raw, "'. EpiExposure v1 supports: beta, ",
      "binomial, poisson, gamma, gaussian, and negative_binomial (NB2)."
    )
  }


  extract_input_link <- function(family_input) {
    if (is.list(family_input) && !is.null(family_input$link) &&
        length(family_input$link) >= 1L) {
      out <- tolower(trimws(as.character(family_input$link[[1]])))
      if (!is.na(out) && nzchar(out)) return(out)
    }
    NULL
  }


  default_family_link <- function(family_name) {
    switch(
      family_name,
      beta = "logit",
      binomial = "logit",
      poisson = "log",
      gamma = "log",
      gaussian = "identity",
      negative_binomial = "log",
      stopf("Could not determine the default link for family '", family_name, "'.")
    )
  }


  validate_link_name <- function(link_name) {
    if (!is_scalar_string(link_name)) stopf("Could not determine a valid model link.")
    link_name <- tolower(link_name)
    link_name
  }


  validate_family_link <- function(
    family_name,
    link_name,
    model_engine
  ) {
    engine_links <- list(
      glm = list(
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity")
      ),
      glmmTMB = list(
        beta = c("logit", "probit", "cloglog", "identity", "inverse", "sqrt"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      gam = list(
        beta = c("logit", "probit", "cloglog", "cauchit"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      gamm = list(
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity")
      ),
      gls = list(
        gaussian = "identity"
      ),
      spamm = list(
        beta = c("logit", "probit", "cloglog", "cauchit"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      brms = list(
        beta = c("logit", "probit", "cloglog", "cauchit"),
        gaussian = c("identity", "log", "inverse"),
        binomial = c("logit", "probit", "cauchit", "log", "cloglog"),
        poisson = c("log", "identity", "sqrt"),
        gamma = c("log", "inverse", "identity"),
        negative_binomial = c("log", "identity", "sqrt")
      ),
      inla = list(
        beta = "logit",
        gaussian = "identity",
        binomial = "logit",
        poisson = "log",
        gamma = "log",
        negative_binomial = "log"
      ),
      bdlnm = list(
        beta = "logit",
        gaussian = "identity",
        binomial = "logit",
        poisson = "log",
        gamma = "log",
        negative_binomial = "log"
      )
    )


    registered_families <- engine_links[[model_engine]]


    if (is.null(registered_families)) {
      stopf("Unsupported `model_engine`: '", model_engine, "'.")
    }


    allowed <- registered_families[[family_name]]


    if (is.null(allowed)) {
      stopf(
        "Family '", family_name, "' is not supported for ",
        "model_engine = '", model_engine, "'. ",
        "Supported families: ",
        paste(names(registered_families), collapse = ", "), "."
      )
    }


    # GLS is intentionally restricted by the EpiExposure contract.
    if (!link_name %in% allowed) {
      stopf(
        "Link '", link_name,
        "' is not supported for family '", family_name,
        "' with `model_engine = '", model_engine, "'. ",
        "Supported link(s): ",
        paste(allowed, collapse = ", "),
        "."
      )
    }


    invisible(TRUE)
  }


  check_dot_conflicts <- function(
    dots,
    reserved,
    engine
  ) {
    if (!length(dots)) {
      return(
        invisible(TRUE)
      )
    }

    dot_names <- names(
      dots
    )

    if (is.null(dot_names) ||
        anyNA(dot_names) ||
        any(!nzchar(dot_names))) {
      stopf(
        "All arguments supplied through `...` must be explicitly named."
      )
    }

    if (anyDuplicated(dot_names)) {
      duplicated_dot_names <- unique(
        dot_names[
          duplicated(dot_names)
        ]
      )

      stopf(
        "Arguments supplied through `...` must have unique names. ",
        "Duplicated argument(s): ",
        paste(
          duplicated_dot_names,
          collapse = ", "
        ),
        "."
      )
    }

    conflict <- intersect(
      dot_names,
      reserved
    )

    if (length(conflict)) {
      stopf(
        "Argument(s) managed internally by EpiExposure cannot be supplied ",
        "in `...` for engine '",
        engine,
        "': ",
        paste(
          conflict,
          collapse = ", "
        ),
        "."
      )
    }

    invisible(TRUE)
  }


  fit_with_context <- function(expr, engine, family_name, link_name) {
    tryCatch(
      expr,
      error = function(e) {
        stopf(
          "Model fitting failed for model_engine = '", engine,
          "', family = '", family_name, "' and link = '", link_name, "'.\n",
          "Original engine error: ", conditionMessage(e)
        )
      }
    )
  }

  # BASIC VALIDATION

  # Preserve the caller's original data expression for compact model calls.
  # This prevents `do.call()` from leaving the fully evaluated data.frame
  # embedded in engine call objects (notably `glmmTMB`), which would otherwise
  # make `summary()` print the complete design and all EpiExposure attributes.
  data_call <- substitute(data)
  family_call <- substitute(family)


  if (!is.data.frame(data)) stopf("`data` must be a data.frame.")
  if (!nrow(data)) stopf("`data` must contain at least one row.")


  if (!"y_model" %in% names(data)) {
    stopf("`data` must contain a column named 'y_model'.")
  }
  if (!is.numeric(data$y_model)) stopf("`y_model` must be numeric.")
  if (any(!is.finite(data$y_model))) {
    stopf("`y_model` must contain only finite values before model fitting.")
  }


  model_engine <- match.arg(
    model_engine,
    choices = c(
      "glm", "glmmTMB", "gam", "gamm", "gls",
      "spamm", "brms", "inla", "bdlnm"
    )
  )


  dots <- list(...)


  # Resolve spaMM spatial specification without changing the legacy intercept.
  # The shared helper is also used for the early find_bestfit() preflight.
  spatial_spec <- .epix_validate_spatial_spec(
    data = data,
    model_engine = model_engine,
    spatial_effect = spatial_effect,
    spatial_structure = spatial_structure,
    spatial_group = spatial_group,
    forbidden = c("y_model", "time")
  )
  spatial_effect <- spatial_spec$effect
  spatial_structure <- spatial_spec$structure
  spatial_group <- spatial_spec$group
  spatial_term <- spatial_spec$term


  family_info <- resolve_family_info(family)
  family_name <- family_info$name


  if (identical(family_name, "ordinal")) {
    stopf(
      "Ordinal outcomes are not supported in EpiExposure v1. This deliberate ",
      "restriction keeps fitting, prediction, uncertainty, performance metrics, ",
      "and ensemble behavior harmonized across supported families."
    )
  }


  if (identical(family_info$variant, "NB1")) {
    stopf(
      "Negative-binomial NB1 was requested, but EpiExposure v1 standardizes ",
      "`negative_binomial` to the quadratic-variance NB2 parameterization. ",
      "Use an NB2-compatible family specification."
    )
  }


  input_link <- extract_input_link(family)


  link_name <- validate_link_name(
    input_link %||% default_family_link(family_name)
  )


  validate_family_link(
    family_name = family_name,
    link_name = link_name,
    model_engine = model_engine
  )


  link_source <- if (is.null(input_link)) {
    "epiexposure_default"
  } else {
    "family_input"
  }

  # RESPONSE VALIDATION BY CANONICAL FAMILY

  y <- data$y_model


  if (family_name == "beta" && any(y <= 0 | y >= 1)) {
    stopf("`family = 'beta'` requires every `y_model` value to lie strictly in (0, 1).")
  }


  if (family_name == "binomial" && !all(y %in% c(0, 1))) {
    stopf(
      "`family = 'binomial'` in EpiExposure v1 represents a Bernoulli 0/1 ",
      "outcome; every `y_model` value must be coded 0 or 1."
    )
  }


  if (family_name %in% c("poisson", "negative_binomial") &&
      (any(y < 0) || any(y != floor(y)))) {
    stopf(
      "`family = '", family_name,
      "'` requires non-negative integer counts in `y_model`."
    )
  }


  if (family_name == "gamma" && any(y <= 0)) {
    stopf("`family = 'gamma'` requires strictly positive `y_model` values.")
  }

  # CHECK CONSISTENCY WITH prepare_response()

  prepared_family_name <- attr(data, "response_family_name")
  if (!is.null(prepared_family_name)) {
    prepared_info <- resolve_family_info(prepared_family_name)
    if (identical(prepared_info$name, "ordinal")) {
      stopf("The response was prepared as ordinal, which is not supported in EpiExposure v1.")
    }
    if (!identical(prepared_info$name, family_name)) {
      stopf(
        "The response was prepared for family '", prepared_info$name,
        "', but the model is being fitted with family '", family_name, "'."
      )
    }
  }

  # RANDOM-INTERCEPT CONTRACT

  if (!is.null(random_effect)) {
    if (!is_scalar_string(random_effect)) {
      stopf("`random_effect` must be NULL or one non-empty column name.")
    }
    if (!random_effect %in% names(data)) {
      stopf("`random_effect` ('", random_effect, "') was not found in `data`.")
    }
    if (anyNA(data[[random_effect]])) {
      stopf("`random_effect` cannot contain missing values.")
    }
    if (length(unique(data[[random_effect]])) < 2L) {
      stopf("`random_effect` must contain at least two observed grouping levels.")
    }
  }


  if (!is.null(random_effect) && model_engine %in% c("glm", "gls")) {
    stopf(
      "`model_engine = '", model_engine,
      "'` does not implement the EpiExposure random-intercept interface. ",
      "Use `random_effect = NULL` or select a supported mixed-model engine."
    )
  }


  if (identical(model_engine, "gamm") && is.null(random_effect)) {
    stopf(
      "`model_engine = 'gamm'` requires `random_effect` in the current ",
      "EpiExposure interface because the generated model contains no other ",
      "smooth/random term."
    )
  }


  # `mgcv` random-effect smooths and `nlme` random intercepts are most robust
  # when the grouping variable is represented as a factor. Preserve the user's
  # original data object and convert only the engine-specific fitting copy.
  fit_data <- data
  fit_random_effect <- random_effect


  # Convert the field-replication index in the fit copy only. It is not
  # interchangeable with, or a replacement for, the conventional intercept.
  if (!is.null(spatial_group)) {
    fit_data[[spatial_group]] <- factor(fit_data[[spatial_group]])
  }


  if (!is.null(random_effect) && model_engine %in% c("gam", "gamm")) {
    fit_data[[random_effect]] <- factor(fit_data[[random_effect]])
  }


  # INLA latent iid effects are indexed internally. Re-index arbitrary user
  # labels to consecutive integers for fitting while preserving the original
  # grouping-column name in EpiExposure metadata. Downstream v1 predictions
  # are population-level and therefore never require a group-specific index.
  if (!is.null(random_effect) && model_engine %in% c("inla", "bdlnm")) {
    internal_name <- ".epiexposure_re_index"
    while (internal_name %in% names(fit_data)) {
      internal_name <- paste0(internal_name, "_")
    }
    level_order <- unique(as.character(fit_data[[random_effect]]))
    fit_data[[internal_name]] <- match(as.character(fit_data[[random_effect]]), level_order)
    fit_random_effect <- internal_name
  }

  # RANDOM-EFFECT PRIOR CONTRACT

  if (!is.null(random_effect_prior)) {


    if (is.null(random_effect)) {
      stopf(
        "`random_effect_prior` requires a non-NULL `random_effect`."
      )
    }


    if (!model_engine %in% c(
      "inla",
      "bdlnm"
    )) {
      stopf(
        "`random_effect_prior` is currently supported only for ",
        "`model_engine = 'inla'` or `model_engine = 'bdlnm'`."
      )
    }


    if (!is.list(random_effect_prior) ||
        !length(random_effect_prior) ||
        is.null(names(random_effect_prior)) ||
        anyNA(names(random_effect_prior)) ||
        any(!nzchar(names(random_effect_prior))) ||
        anyDuplicated(names(random_effect_prior))) {

      stopf(
        "`random_effect_prior` must be NULL or a non-empty named list ",
        "accepted by `INLA::f(..., hyper = ...)`."
      )
    }
  }

  # NORMALIZE AND VALIDATE EXPOSURE SPECIFICATION

  normalize_spec <- function(spec) {
    if (is.null(spec)) {
      stopf(
        "`epiexposure_spec` is required. EpiExposure downstream functions must ",
        "reuse the fitted cross-basis specification rather than reconstructing it ",
        "from user-supplied degrees of freedom or defaults."
      )
    }


    if (!is.list(spec) || is.null(names(spec)) || !length(spec) ||
        anyNA(names(spec)) || any(names(spec) == "") || anyDuplicated(names(spec))) {
      stopf("`epiexposure_spec` must be a non-empty named list with unique names.")
    }


    for (nm in names(spec)) {
      current <- spec[[nm]]


      if (!is.list(current)) {
        stopf("Specification for variable '", nm, "' must be a list.")
      }


      max_lag_value <- current$max_lag


      if (is.null(max_lag_value) ||
          !is.numeric(max_lag_value) ||
          length(max_lag_value) != 1L ||
          is.na(max_lag_value) ||
          !is.finite(max_lag_value) ||
          max_lag_value < 0 ||
          max_lag_value != as.integer(max_lag_value)) {
        stopf(
          "`max_lag` for variable '", nm,
          "' must be one non-negative finite integer."
        )
      }


      if (is.null(current$argvar) || !is.list(current$argvar)) {
        stopf("Missing or invalid `argvar` for variable '", nm, "'.")
      }


      if (is.null(current$arglag) || !is.list(current$arglag)) {
        stopf("Missing or invalid `arglag` for variable '", nm, "'.")
      }


      current$max_lag <- as.integer(max_lag_value)
      spec[[nm]] <- current
    }


    fitted_max_lags <- vapply(
      spec,
      function(current) current$max_lag,
      integer(1)
    )


    if (length(unique(fitted_max_lags)) != 1L) {
      stopf(
        "All fitted exposure variables must use the same `max_lag` under the ",
        "EpiExposure exact-profile contract. Received: ",
        paste(
          paste0(names(fitted_max_lags), "=", fitted_max_lags),
          collapse = ", "
        ),
        "."
      )
    }


    spec
  }


  epiexposure_spec <- normalize_spec(epiexposure_spec)


  common_max_lag <- epiexposure_spec[[1L]]$max_lag
  expected_history_length <- common_max_lag + 1L

  fitted_time_step <- attr(
    data,
    "epiexposure_time_step",
    exact = TRUE
  )

  if (expected_history_length > 1L) {
    if (
      !is.numeric(fitted_time_step) ||
      length(fitted_time_step) != 1L ||
      is.na(fitted_time_step) ||
      !is.finite(fitted_time_step) ||
      fitted_time_step <= 0
    ) {
      stopf(
        "`data` is missing valid `epiexposure_time_step` metadata. ",
        "Rebuild the epidemic-level design with the current ",
        "`build_design()` implementation before fitting."
      )
    }

    fitted_time_step <- as.numeric(
      fitted_time_step
    )
  } else {
    if (
      !is.null(fitted_time_step) &&
      (
        !is.numeric(fitted_time_step) ||
        length(fitted_time_step) != 1L ||
        !is.na(fitted_time_step)
      )
    ) {
      stopf(
        "For `max_lag = 0`, `epiexposure_time_step` must be `NA` because ",
        "a one-observation profile does not define a temporal interval."
      )
    }

    fitted_time_step <- NA_real_
  }

  # CROSS-BASIS COLUMN / VARIABLE CONTRACT

  parse_cb_variable <- function(columns) {
    out <- sub("^cb_", "", columns)
    sub("_[0-9]+$", "", out)
  }


  sort_cb_cols <- function(columns) {
    if (!length(columns)) return(columns)
    variables <- parse_cb_variable(columns)
    variable_order <- unique(variables)
    index <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", columns)))
    if (anyNA(index)) {
      stopf(
        "Every EpiExposure cross-basis column must end in a numeric index: ",
        paste(columns[is.na(index)], collapse = ", "), "."
      )
    }
    columns[order(match(variables, variable_order), index)]
  }


  raw_cb_cols <- grep("^cb_", names(fit_data), value = TRUE)
  cb_cols <- sort_cb_cols(raw_cb_cols)


  if (!length(cb_cols) && model_engine != "bdlnm") {
    stopf("No `cb_*` columns were found in `data`.")
  }


  if (length(cb_cols)) {
    invalid_cb <- vapply(
      fit_data[cb_cols],
      function(x) !is.numeric(x) || any(!is.finite(x)),
      logical(1)
    )
    if (any(invalid_cb)) {
      stopf(
        "All fitted `cb_*` columns must contain finite numeric values. Invalid: ",
        paste(names(invalid_cb)[invalid_cb], collapse = ", "), "."
      )
    }
  }


  validate_basis_objects <- function(x, required = FALSE) {
    if (is.null(x)) {
      if (required) {
        stopf(
          "For `model_engine = 'bdlnm'`, `basis_objects` must be a non-empty ",
          "named list of epidemic-level, matrix-form `dlnm::crossbasis()` ",
          "objects containing one row per model row."
        )
      }


      return(NULL)
    }


    if (is.data.frame(x)) {
      stopf(
        "`basis_objects` must be a named list of `dlnm::crossbasis()` objects, ",
        "not a data.frame or an epidemic-level design matrix."
      )
    }


    if (!is.list(x) ||
        is.null(names(x)) ||
        !length(x) ||
        anyNA(names(x)) ||
        any(names(x) == "") ||
        anyDuplicated(names(x))) {
      stopf(
        "`basis_objects` must be NULL or a non-empty named list with ",
        "unique exposure-variable names."
      )
    }


    valid_crossbasis <- vapply(
      x,
      inherits,
      logical(1),
      what = "crossbasis"
    )


    if (any(!valid_crossbasis)) {
      stopf(
        "Every element of `basis_objects` must inherit from 'crossbasis'. ",
        "Invalid object(s): ",
        paste(names(x)[!valid_crossbasis], collapse = ", "),
        "."
      )
    }


    x
  }


  # Recover basis objects stored by build_design() when they are not
  # supplied explicitly by the user.
  basis_objects_source <- if (is.null(basis_objects)) {
    if (identical(model_engine, "bdlnm")) {
      "data_bdlnm_attribute"
    } else {
      "data_attribute"
    }
  } else {
    "user_argument"
  }


  if (is.null(basis_objects)) {


    if (identical(model_engine, "bdlnm")) {


      basis_objects <- attr(
        data,
        "epiexposure_bdlnm_basis_objects",
        exact = TRUE
      )


      if (is.null(basis_objects)) {
        stopf(
          "`model_engine = 'bdlnm'` requires epidemic-level cross-basis ",
          "objects produced by the current `build_design()` pipeline. ",
          "Rebuild the design with `build_design()` before fitting."
        )
      }


    } else {


      basis_objects <- attr(
        data,
        "epiexposure_basis_objects",
        exact = TRUE
      )
    }
  }

  if (is.null(basis_objects)) {
    basis_objects_source <- "none"
  }

  # Validate the final object regardless of whether it was supplied
  # explicitly or recovered from the design metadata.
  basis_objects <- validate_basis_objects(
    basis_objects,
    required = identical(model_engine, "bdlnm")
  )
  vars_inferred <- if (length(cb_cols)) {
    unique(parse_cb_variable(cb_cols))
  } else {
    names(basis_objects)
  }


  if (!length(vars_inferred)) {
    stopf("Could not determine the exposure variables used by the fitted model.")
  }


  if (!setequal(names(epiexposure_spec), vars_inferred)) {
    missing_spec <- setdiff(vars_inferred, names(epiexposure_spec))
    extra_spec <- setdiff(names(epiexposure_spec), vars_inferred)
    parts <- character(0)
    if (length(missing_spec)) {
      parts <- c(parts, paste0("missing: ", paste(missing_spec, collapse = ", ")))
    }
    if (length(extra_spec)) {
      parts <- c(parts, paste0("extra: ", paste(extra_spec, collapse = ", ")))
    }
    stopf(
      "`epiexposure_spec` names must match the fitted exposure variables exactly (",
      paste(parts, collapse = "; "), ")."
    )
  }
  epiexposure_spec <- epiexposure_spec[vars_inferred]


  if (!is.null(basis_objects)) {
    if (!setequal(names(basis_objects), vars_inferred)) {
      stopf(
        "`basis_objects` names must match the fitted exposure variables exactly. ",
        "Expected: ", paste(vars_inferred, collapse = ", "), "."
      )
    }


    basis_objects <- basis_objects[vars_inferred]


    for (nm in names(basis_objects)) {
      basis_object <- basis_objects[[nm]]


      if (inherits(basis_object, "crossbasis")) {
        basis_lag <- attr(basis_object, "lag")


        if (is.null(basis_lag) ||
            !is.numeric(basis_lag) ||
            !length(basis_lag) ||
            anyNA(basis_lag) ||
            any(!is.finite(basis_lag))) {
          stopf(
            "Cross-basis object for exposure '", nm,
            "' has invalid `lag` metadata."
          )
        }


        basis_max_lag <- as.integer(max(basis_lag))


        if (basis_max_lag != common_max_lag) {
          stopf(
            "Cross-basis lag metadata for exposure '", nm,
            "' does not match the common fitted `max_lag`. Expected ",
            common_max_lag, " but found ", basis_max_lag, "."
          )
        }
      }
    }
  }


  if (identical(model_engine, "bdlnm")) {


    invalid_nrow <- vapply(
      basis_objects,
      function(x) nrow(x) != nrow(fit_data),
      logical(1)
    )


    if (any(invalid_nrow)) {
      stopf(
        "For `model_engine = 'bdlnm'`, every cross-basis object must ",
        "contain one row per epidemic-level model row. `data` has ",
        nrow(fit_data),
        " rows, but incompatible basis object(s) were found: ",
        paste(
          names(basis_objects)[invalid_nrow],
          collapse = ", "
        ),
        "."
      )
    }

    # Ensure that every bdlnm cross-basis has valid and globally unique
    # internal column names. Native dlnm names such as v1.l1 are repeated
    # across exposure-specific cross-bases and may collide inside INLA.

    for (nm in names(basis_objects)) {

      current_basis <- basis_objects[[nm]]
      current_names <- colnames(
        current_basis
      )

      if (is.null(current_names)) {
        current_names <- paste0(
          "basis_",
          seq_len(
            ncol(current_basis)
          )
        )
      }

      if (length(current_names) != ncol(current_basis) ||
          anyNA(current_names) ||
          any(!nzchar(current_names)) ||
          anyDuplicated(current_names)) {

        stopf(
          "Cross-basis object for exposure '",
          nm,
          "' contains invalid internal column names."
        )
      }

      expected_prefix <- paste0(
        nm,
        "_"
      )

      has_expected_prefix <- startsWith(
        current_names,
        expected_prefix
      )

      if (any(has_expected_prefix) &&
          !all(has_expected_prefix)) {
        stopf(
          "Cross-basis object for exposure '",
          nm,
          "' contains partially prefixed internal column names. ",
          "Column names must either all include the exposure prefix or ",
          "all use the native cross-basis names."
        )
      }

      if (!any(has_expected_prefix)) {
        current_names <- paste0(
          expected_prefix,
          current_names
        )
      }

      colnames(
        current_basis
      ) <- current_names

      basis_objects[[nm]] <- current_basis
    }

    all_bdlnm_basis_names <- unlist(
      lapply(
        basis_objects,
        colnames
      ),
      use.names = FALSE
    )

    if (!length(all_bdlnm_basis_names) ||
        anyNA(all_bdlnm_basis_names) ||
        any(!nzchar(all_bdlnm_basis_names)) ||
        anyDuplicated(all_bdlnm_basis_names)) {

      duplicated_names <- unique(
        all_bdlnm_basis_names[
          duplicated(all_bdlnm_basis_names)
        ]
      )

      if (length(duplicated_names)) {
        stopf(
          "The bdlnm cross-basis objects contain duplicated internal ",
          "column names: ",
          paste(
            duplicated_names,
            collapse = ", "
          ),
          "."
        )
      }

      stopf(
        "The bdlnm cross-basis objects contain missing or empty ",
        "internal column names."
      )
    }
  }

  # ENGINE-SPECIFIC FAMILY CONSTRUCTION

  resolve_engine_family <- function(engine, family_name, link_name) {


    if (engine == "glm") {
      return(switch(
        family_name,
        gaussian = stats::gaussian(link = link_name),
        poisson = stats::poisson(link = link_name),
        gamma = stats::Gamma(link = link_name),
        binomial = stats::binomial(link = link_name),
        stopf(
          "Family '", family_name,
          "' is not supported by `model_engine = 'glm'` in EpiExposure v1."
        )
      ))
    }


    if (engine == "glmmTMB") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stopf("Package 'glmmTMB' is required for `model_engine = 'glmmTMB'`.")
      }
      return(switch(
        family_name,
        beta = glmmTMB::beta_family(link = link_name),
        gaussian = stats::gaussian(link = link_name),
        poisson = stats::poisson(link = link_name),
        gamma = stats::Gamma(link = link_name),
        binomial = stats::binomial(link = link_name),
        negative_binomial = glmmTMB::nbinom2(link = link_name),
        stopf("Unsupported glmmTMB family.")
      ))
    }


    if (engine == "gam") {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stopf("Package 'mgcv' is required for `model_engine = 'gam'`.")
      }
      return(switch(
        family_name,
        beta = do.call(mgcv::betar,list(link = link_name)),
        gaussian = stats::gaussian(link = link_name),
        poisson = stats::poisson(link = link_name),
        gamma = stats::Gamma(link = link_name),
        binomial = stats::binomial(link = link_name),
        negative_binomial = do.call(mgcv::nb,list(link = link_name)),
        stopf("Unsupported GAM family.")
      ))
    }


    if (engine == "gamm") {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stopf("Package 'mgcv' is required for `model_engine = 'gamm'`.")
      }
      return(switch(
        family_name,
        gaussian = stats::gaussian(link = link_name),
        poisson = stats::poisson(link = link_name),
        gamma = stats::Gamma(link = link_name),
        binomial = stats::binomial(link = link_name),
        stopf(
          "Family '", family_name, "' is not harmonized for `gamm`. ",
          "In EpiExposure v1, `gamm` supports Gaussian, Binomial, Poisson, and Gamma only."
        )
      ))
    }


    if (engine == "gls") {
      if (family_name != "gaussian") {
        stopf("`model_engine = 'gls'` supports only `family = 'gaussian'`.")
      }
      if (link_name != "identity") {
        stopf("`model_engine = 'gls'` requires the Gaussian identity link.")
      }
      return(stats::gaussian(link = "identity"))
    }


    if (engine == "spamm") {
      if (!requireNamespace("spaMM", quietly = TRUE)) {
        stopf("Package 'spaMM' is required for `model_engine = 'spamm'`.")
      }
      return(switch(
        family_name,
        beta = spaMM::beta_resp(link = link_name),
        gaussian = stats::gaussian(link = link_name),
        poisson = stats::poisson(link = link_name),
        gamma = stats::Gamma(link = link_name),
        binomial = stats::binomial(link = link_name),
        negative_binomial = spaMM::negbin2(link = link_name),
        stopf("Unsupported spaMM family.")
      ))
    }


    if (engine == "brms") {
      if (!requireNamespace("brms", quietly = TRUE)) {
        stopf("Package 'brms' is required for `model_engine = 'brms'`.")
      }
      return(switch(
        family_name,
        beta = brms::Beta(link = link_name),
        gaussian = stats::gaussian(link = link_name),
        poisson = stats::poisson(link = link_name),
        gamma = stats::Gamma(link = link_name),
        binomial = brms::bernoulli(link = link_name),
        negative_binomial = brms::negbinomial(link = link_name),
        stopf("Unsupported brms family.")
      ))
    }


    if (engine %in% c("inla", "bdlnm")) {
      return(switch(
        family_name,
        beta = "beta",
        gaussian = "gaussian",
        poisson = "poisson",
        gamma = "gamma",
        binomial = "binomial",
        negative_binomial = "nbinomial",
        stopf("Unsupported INLA-backed family.")
      ))
    }


    stopf("Unsupported `model_engine`: ", engine, ".")
  }


  engine_family <- resolve_engine_family(model_engine, family_name, link_name)


  if (model_engine %in% c("inla", "bdlnm")) {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stopf(
        "Package 'INLA' is required for ",
        "`model_engine = '", model_engine, "'."
      )
    }

    if (identical(model_engine, "bdlnm")) {
      if (!requireNamespace("bdlnm", quietly = TRUE)) {
        stopf(
          "Package 'bdlnm' is required for ",
          "`model_engine = 'bdlnm'`."
        )
      }

      if (!requireNamespace("sn", quietly = TRUE)) {
        stopf(
          "Package 'sn' is required for posterior sampling with ",
          "`model_engine = 'bdlnm'`. Install it with ",
          "`install.packages(\"sn\")`."
        )
      }
    }

    available_likelihoods <- names(
      INLA::inla.models()$likelihood
    )

    if (!engine_family %in% available_likelihoods) {
      stopf(
        "Likelihood '", engine_family,
        "' for EpiExposure family '", family_name,
        "' and model_engine = '", model_engine,
        "' is unavailable. ",
        "The current INLA installation does not provide this likelihood."
      )
    }
  }

  # ENGINE-SPECIFIC DOTS / LINK CONTROL

  if (model_engine == "gam") {
    if (!is.null(dots$method) && !is_scalar_string(dots$method)) {
      stopf("`method` supplied to `gam` must be one non-empty character value.")
    }


    # mgcv extended families have restricted smoothing-parameter estimation
    # routes. REML is also a stable default for random-effect smooths.
    if (is.null(dots$method) &&
        (family_name %in% c("beta", "negative_binomial") || !is.null(random_effect))) {
      dots$method <- "REML"
    }


    if (family_name == "beta" && !is.null(dots$method) &&
        !toupper(dots$method) %in% c("REML", "ML", "NCV")) {
      stopf(
        "`gam` with Beta family requires a supported extended-family method: ",
        "'REML', 'ML', or 'NCV'."
      )
    }


    if (family_name == "negative_binomial" && !is.null(dots$method) &&
        !toupper(dots$method) %in% c("REML", "NCV")) {
      stopf(
        "`gam` with `mgcv::nb()` requires `method = 'REML'` or `method = 'NCV'`."
      )
    }
  }


  prepare_inla_family_control <- function(input_dots) {
    control_family <- input_dots$control.family %||% list()
    if (!is.list(control_family)) {
      stopf("INLA `control.family` supplied through `...` must be a list.")
    }


    control_link <- control_family$control.link %||% list()
    if (!is.list(control_link)) {
      stopf("INLA `control.family$control.link` must be a list.")
    }


    if (!is.null(control_link$model)) {
      supplied_link <- tolower(as.character(control_link$model)[1])
      if (!identical(supplied_link, link_name)) {
        stopf(
          "Conflicting INLA link specifications: `family` implies '", link_name,
          "' but `control.family$control.link$model` is '", supplied_link, "'."
        )
      }
    }


    control_link$model <- link_name
    control_family$control.link <- control_link
    input_dots$control.family <- control_family
    input_dots
  }


  if (model_engine %in% c("inla", "bdlnm")) {
    dots <- prepare_inla_family_control(dots)
  }


  if (model_engine %in% c("inla", "bdlnm")) {
    control_compute <- dots$control.compute %||% list()
    if (!is.list(control_compute)) {
      stopf("INLA-backed engines require `control.compute` supplied through `...` to be a list.")
    }
    if (!is.null(control_compute$config) && !isTRUE(control_compute$config)) {
      stopf(
        "EpiExposure requires `control.compute$config = TRUE` for INLA-backed engines so ",
        "posterior coefficient draws can be generated for downstream uncertainty."
      )
    }
    control_compute$config <- TRUE
    dots$control.compute <- control_compute
  }

  # INLA RANDOM-EFFECT PRIOR ENVIRONMENT

  inla_formula_env <- new.env(
    parent = parent.frame()
  )


  if (!is.null(random_effect_prior)) {


    assign(
      ".epiexposure_random_effect_prior",
      random_effect_prior,
      envir = inla_formula_env
    )
  }

  # FORMULA CONSTRUCTION

  cb_rhs <- paste(vapply(cb_cols, quote_name, character(1)), collapse = " + ")
  fixed_formula_text <- if (nzchar(cb_rhs)) {
    paste0("y_model ~ 1 + ", cb_rhs)
  } else {
    "y_model ~ 1"
  }


  fixed_formula <- stats::as.formula(fixed_formula_text)


  mixed_formula <- fixed_formula
  gam_formula <- fixed_formula
  inla_formula <- fixed_formula


  if (!is.null(random_effect)) {
    qre <- quote_name(fit_random_effect)


    if (model_engine %in% c("glmmTMB", "spamm", "brms")) {
      mixed_formula <- stats::as.formula(
        paste0(fixed_formula_text, " + (1 | ", qre, ")")
      )
    }


    if (model_engine == "gam") {
      gam_formula <- stats::as.formula(
        paste0(fixed_formula_text, " + s(", qre, ", bs = 're')")
      )
    }


    if (model_engine == "inla") {


      if (is.null(random_effect_prior)) {


        random_term <- paste0(
          "f(",
          qre,
          ", model = 'iid')"
        )


      } else {


        random_term <- paste0(
          "f(",
          qre,
          ", model = 'iid', ",
          "hyper = .epiexposure_random_effect_prior)"
        )
      }


      inla_formula <- stats::as.formula(
        paste0(
          fixed_formula_text,
          " + ",
          random_term
        ),
        env = inla_formula_env
      )
    }
  }


  # The spaMM formula receives the spatial special term only when requested.
  # The spaMM namespace is used as its lexical parent, making its formula
  # handlers accessible without requiring users to attach library(spaMM).
  if (!is.null(spatial_term)) {
    if (!requireNamespace("spaMM", quietly = TRUE)) {
      stopf("Package 'spaMM' is required for spatial-effect fitting.")
    }
    spaMM_formula_text <- fixed_formula_text
    if (!is.null(random_effect)) {
      spaMM_formula_text <- paste0(
        spaMM_formula_text, " + (1 | ", quote_name(fit_random_effect), ")"
      )
    }
    mixed_formula <- stats::as.formula(
      paste0(spaMM_formula_text, " + ", spatial_term),
      env = new.env(parent = asNamespace("spaMM"))
    )
  }

  # STANDARDIZED METADATA

  family_parameterization <- switch(
    family_name,
    beta = "mean_precision",
    binomial = "Bernoulli_0_1",
    poisson = "Poisson_mean",
    gamma = "Gamma_mean_link",
    gaussian = "Gaussian_mean",
    negative_binomial = "NB2",
    NA_character_
  )


  random_structure <- if (is.null(random_effect)) {
    "none"
  } else {
    "random_intercept"
  }


  attach_epiexposure_meta <- function(model_obj) {
    basis_meta <- basis_objects
    if (identical(model_engine, "bdlnm") &&
        is.list(model_obj) && !is.null(model_obj$basis)) {
      basis_meta <- model_obj$basis
    }


    attr(model_obj, "epiexposure_engine") <- model_engine
    attr(model_obj, "epiexposure_family_input") <- family
    attr(model_obj, "epiexposure_family") <- engine_family
    attr(model_obj, "epiexposure_family_name") <- family_name
    attr(model_obj, "epiexposure_family_parameterization") <- family_parameterization
    attr(model_obj, "epiexposure_link") <- link_name
    attr(model_obj, "epiexposure_link_source") <- link_source
    attr(model_obj, "epiexposure_cb_cols") <- cb_cols
    attr(model_obj, "epiexposure_vars") <- vars_inferred
    attr(model_obj, "epiexposure_data_template") <- fit_data[1, , drop = FALSE]
    attr(model_obj, "epiexposure_id_col") <- random_effect
    attr(model_obj, "epiexposure_random_effect_fit_col") <- fit_random_effect
    attr(model_obj, "epiexposure_random_structure") <- random_structure
    attr(model_obj, "epiexposure_spatial_effect") <- spatial_effect
    attr(model_obj, "epiexposure_spatial_structure") <- spatial_structure
    attr(model_obj, "epiexposure_spatial_group") <- spatial_group
    attr(model_obj, "epiexposure_spatial_term") <- spatial_term
    attr(model_obj, "epiexposure_has_spatial_effect") <- !is.null(spatial_effect)
    attr(model_obj, "spatial_effect") <- spatial_effect
    attr(model_obj, "spatial_structure") <- spatial_structure
    attr(model_obj, "spatial_group") <- spatial_group
    attr(model_obj, "spatial_term") <- spatial_term
    attr(model_obj, "epiexposure_spec") <- epiexposure_spec
    attr(model_obj, "epiexposure_basis_objects") <- basis_meta
    attr(model_obj, "epiexposure_basis_objects_source") <- basis_objects_source
    attr(model_obj, "epiexposure_max_lag") <- common_max_lag
    attr(model_obj, "epiexposure_history_length") <- expected_history_length
    attr(model_obj,"epiexposure_history_contract") <- "all_fitted_exposures_same_exact_max_lag_plus_one"
    attr(model_obj,"epiexposure_time_step") <- fitted_time_step
    attr(model_obj, "epiexposure_prediction_level") <- "population"
    attr(model_obj, "epiexposure_prediction_estimand") <- "expected_response"
    attr(model_obj, "epiexposure_point_prediction_contract") <- "central_expected_response"
    attr(model_obj, "epiexposure_uncertainty_contract") <- "draw_by_draw_median_quantiles"
    attr(model_obj, "epiexposure_random_effect_prior") <- random_effect_prior
    attr(model_obj, "epiexposure_random_effect_model") <- if (
      !is.null(random_effect) && model_engine %in% c("inla","bdlnm")
    ) {
      "iid"
    } else {
      NULL
    }

    model_obj
  }

  # FIT: glm

  if (model_engine == "glm") {
    check_dot_conflicts(
      dots,
      c(
        "formula",
        "data",
        "family"
      ),
      model_engine
    )


    args <- c(
      list(
        formula = fixed_formula,
        data = fit_data,
        family = engine_family
      ),
      dots
    )


    model_obj <- fit_with_context(
      do.call(
        stats::glm,
        args
      ),
      model_engine,
      family_name,
      link_name
    )


    # `do.call()` receives evaluated arguments and may store the complete
    # function and family object inside the fitted call. Rebuild only the
    # stored call used for printing and re-evaluation, preserving all other
    # arguments passed through `...`.
    if (!is.null(model_obj$call) &&
        is.call(model_obj$call)) {


      call_parts <- as.list(
        model_obj$call
      )


      # Replace the embedded function definition with a compact function call.
      call_parts[[1L]] <- quote(
        stats::glm
      )


      # Preserve the fitted formula while displaying the original family and
      # data expressions supplied by the user.
      call_parts[["formula"]] <- fixed_formula
      call_parts[["family"]] <- family_call
      call_parts[["data"]] <- data_call


      model_obj$call <- as.call(
        call_parts
      )
    }


    return(
      attach_epiexposure_meta(
        model_obj
      )
    )
  }

  # FIT: glmmTMB

  if (model_engine == "glmmTMB") {
    check_dot_conflicts(dots, c("formula", "data", "family"), model_engine)
    args <- c(
      list(formula = mixed_formula, data = fit_data, family = engine_family),
      dots
    )


    model_obj <- fit_with_context(
      do.call(glmmTMB::glmmTMB, args),
      model_engine,
      family_name,
      link_name
    )


    # `do.call()` evaluates `fit_data` before calling glmmTMB. Consequently,
    # glmmTMB may store the entire evaluated data.frame inside `model_obj$call`.
    # That makes `summary(model_obj)` print hundreds/thousands of lines of
    # design values and attributes before the ordinary model summary.
    #
    # Replace only the stored `data` component of the call with the expression
    # originally supplied by the user (for example, `dat`). The fitted object
    # itself, its TMB structures, estimates, likelihood, vcov, and EpiExposure
    # metadata are untouched.
    if (!is.null(model_obj$call) && is.call(model_obj$call)) {
      call_parts <- as.list(model_obj$call)


      if ("data" %in% names(call_parts)) {
        call_parts[["data"]] <- data_call
        model_obj$call <- as.call(call_parts)
      }
    }


    return(attach_epiexposure_meta(model_obj))
  }

  # FIT: GAM

  if (model_engine == "gam") {
    check_dot_conflicts(dots, c("formula", "data", "family"), model_engine)
    args <- c(
      list(formula = gam_formula, data = fit_data, family = engine_family),
      dots
    )
    return(attach_epiexposure_meta(fit_with_context(
      do.call(mgcv::gam, args),
      model_engine,
      family_name,
      link_name
    )))
  }

  # FIT: GAMM

  if (model_engine == "gamm") {
    check_dot_conflicts(
      dots,
      c(
        "formula",
        "data",
        "family",
        "random"
      ),
      model_engine
    )


    if (family_name == "binomial") {
      warning(
        "`mgcv::gamm()` fits non-Gaussian models by PQL and mgcv specifically ",
        "warns that binary responses may perform poorly. Consider `gam` with a ",
        "random-effect smooth or another GLMM engine when appropriate.",
        call. = FALSE
      )
    }


    random_list <- stats::setNames(
      list(
        stats::as.formula(
          "~1"
        )
      ),
      fit_random_effect
    )


    args <- c(
      list(
        formula = fixed_formula,
        random = random_list,
        data = fit_data,
        family = engine_family
      ),
      dots
    )


    model_obj <- fit_with_context(
      do.call(
        mgcv::gamm,
        args
      ),
      model_engine,
      family_name,
      link_name
    )


    class(model_obj) <- c(
      "epiexposure_gamm",
      class(model_obj)
    )


    return(
      attach_epiexposure_meta(
        model_obj
      )
    )
  }

  # FIT: GLS

  if (model_engine == "gls") {


    if (!requireNamespace("nlme", quietly = TRUE)) {
      stopf("Package 'nlme' is required for `model_engine = 'gls'`.")
    }


    check_dot_conflicts(
      dots,
      c(
        "model",
        "data"
      ),
      model_engine
    )


    args <- c(
      list(
        model = fixed_formula,
        data  = fit_data
      ),
      dots
    )


    model_obj <- fit_with_context(
      do.call(
        nlme::gls,
        args
      ),
      model_engine,
      family_name,
      link_name
    )


    # `do.call()` receives the evaluated `fit_data` object. Consequently,
    # `nlme::gls()` may store the complete data.frame inside `model_obj$call`.
    # This causes `summary(model_obj)` to print the entire evaluated design
    # instead of the compact data expression originally supplied by the user.
    #
    # Replace only the stored `data` component used for display and
    # re-evaluation. This does not change fitted coefficients, residuals,
    # likelihood, covariance structures, predictions, or EpiExposure metadata.


    if (!is.null(model_obj$call) &&
        is.call(model_obj$call)) {


      call_parts <- as.list(
        model_obj$call
      )


      if ("data" %in% names(call_parts)) {


        call_parts[["data"]] <- data_call


        model_obj$call <- as.call(
          call_parts
        )


      }
    }


    return(
      attach_epiexposure_meta(
        model_obj
      )
    )
  }

  # FIT: spaMM

  if (model_engine == "spamm") {
    # `spatial_*` are interface arguments and must not enter `fitme()` dots.
    check_dot_conflicts(dots, c("formula", "data", "family"), model_engine)
    args <- c(
      list(formula = mixed_formula, data = fit_data, family = engine_family),
      dots
    )
    return(attach_epiexposure_meta(fit_with_context(
      do.call(spaMM::fitme, args),
      model_engine,
      family_name,
      link_name
    )))
  }

  # FIT: brms

  if (model_engine == "brms") {
    check_dot_conflicts(dots, c("formula", "data", "family"), model_engine)
    args <- c(
      list(formula = mixed_formula, data = fit_data, family = engine_family),
      dots
    )
    return(attach_epiexposure_meta(fit_with_context(
      do.call(brms::brm, args),
      model_engine,
      family_name,
      link_name
    )))
  }

  # FIT: INLA

  if (model_engine == "inla") {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stopf("Package 'INLA' is required for `model_engine = 'inla'`.")
    }
    check_dot_conflicts(dots, c("formula", "data", "family"), model_engine)
    args <- c(
      list(formula = inla_formula, data = fit_data, family = engine_family),
      dots
    )
    return(attach_epiexposure_meta(fit_with_context(
      do.call(INLA::inla, args),
      model_engine,
      family_name,
      link_name
    )))
  }

  # FIT: Bayesian DLNM (bdlnm)

  if (model_engine == "bdlnm") {

    check_dot_conflicts(dots, c("formula", "data", "family"), model_engine)


    basis_names <- names(basis_objects)
    basis_rhs <- paste(vapply(basis_names, quote_name, character(1)), collapse = " + ")
    bdlnm_formula_text <- paste0("y_model ~ 1 + ", basis_rhs)


    if (!is.null(random_effect)) {


      if (is.null(random_effect_prior)) {


        bdlnm_formula_text <- paste0(
          bdlnm_formula_text,
          " + f(",
          quote_name(fit_random_effect),
          ", model = 'iid')"
        )


      } else {


        bdlnm_formula_text <- paste0(
          bdlnm_formula_text,
          " + f(",
          quote_name(fit_random_effect),
          ", model = 'iid', ",
          "hyper = .epiexposure_random_effect_prior)"
        )
      }
    }


    bdlnm_formula <- stats::as.formula(bdlnm_formula_text)


    # bdlnm requires the original basis objects to be visible from the formula
    # environment. Preserve the caller as parent so user-supplied terms in `...`
    # keep their normal lookup behavior.
    eval_env <- new.env(parent = parent.frame())
    for (nm in basis_names) {
      assign(nm, basis_objects[[nm]], envir = eval_env)
    }


    if (!is.null(random_effect_prior)) {


      assign(
        ".epiexposure_random_effect_prior",
        random_effect_prior,
        envir = eval_env
      )
    }


    environment(bdlnm_formula) <- eval_env


    args <- c(
      list(formula = bdlnm_formula, data = fit_data, family = engine_family),
      dots
    )
    return(attach_epiexposure_meta(fit_with_context(
      #do.call(bdlnm::bdlnm, args),
      do.call(bdlnm::bdlnm, args, envir = eval_env),
      model_engine,
      family_name,
      link_name
    )))
  }


  stopf("Unsupported `model_engine`.")
}
