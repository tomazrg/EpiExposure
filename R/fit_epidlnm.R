#' Fit DLNM inferential model
#'
#' Fits a DLNM inferential model from an epidemic-level design matrix and
#' stores standardized EpiExposure metadata for downstream prediction,
#' effect summarization, and exposure-impact calculations.
#'
#' @param data Data frame containing `y_model` and the fitted `cb_*` terms.
#'   Data prepared with `prepare_response()` may contain the
#'   `response_family_name` attribute, which is checked against `family`.
#' @param model_engine Character scalar identifying the modeling engine:
#'   `"glm"`, `"glmmTMB"`, `"gam"`, `"gamm"`, `"gls"`, `"spamm"`,
#'   `"brms"`, `"inla"`, or `"bdlnm"`.
#' @param family Distribution family supplied as a supported character name or
#'   an engine-compatible family object. Canonical family names are `"beta"`,
#'   `"binomial"`, `"poisson"`, `"gamma"`, `"gaussian"`,
#'   `"negative_binomial"`, and `"ordinal"`. Aliases such as `"negbin"`,
#'   `"nbinom1"`, and `"nbinom2"` are normalized internally.
#' @param random_effect Optional character scalar naming a random-effect column.
#' @param epiexposure_spec Optional named list describing how each cross-basis
#'   was constructed. The recommended structure is:
#'   `list(tmean = list(lag_max = 85, argvar = list(...), arglag = list(...)))`.
#'   If `NULL`, the model can be fitted, but profile-based downstream functions
#'   may not be available.
#' @param basis_objects Optional named list of original cross-basis or one-basis
#'   objects. Required when `model_engine = "bdlnm"`.
#' @param ... Additional arguments passed to the selected modeling engine.
#'
#' @return A fitted model object with standardized EpiExposure attributes,
#'   including the engine, canonical family name, engine-specific family,
#'   link, ordered cross-basis columns, exposure variables, data template,
#'   grouping column, exposure specification, and basis objects.
#'
#' @details
#' This function fits models but does not generate posterior predictions or
#' uncertainty summaries. Consequently, it does not call
#' `posterior_linpred()`, `posterior_epred()`, or `posterior_predict()`, and no
#' single posterior draw is used as a deterministic estimate.
#'
#' Family specifications are normalized to a canonical family name and then
#' converted to the representation expected by the selected engine. Unknown or
#' unsupported engine-family combinations generate explicit errors. No silent
#' fallback to the Gaussian family is used.
#'
#' If `prepare_response()` stored a `response_family_name` attribute in `data`,
#' that family must match the canonical family used for model fitting.
#'
#' @export
fit_epidlnm <- function(
    data,
    model_engine,
    family,
    random_effect = NULL,
    epiexposure_spec = NULL,
    basis_objects = NULL,
    ...
) {

  # =========================================================
  # BASIC VALIDATION
  # =========================================================

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (!"y_model" %in% names(data)) {
    stop("`data` must contain a column named 'y_model'.")
  }

  if (!is.numeric(data$y_model)) {
    stop("`y_model` must be numeric.")
  }

  if (!length(data$y_model)) {
    stop("`y_model` cannot be empty.")
  }

  if (any(!is.finite(data$y_model))) {
    stop("`y_model` must contain only finite values before model fitting.")
  }

  model_engine <- match.arg(
    model_engine,
    choices = c(
      "glm", "glmmTMB", "gam", "gamm", "gls",
      "spamm", "brms", "inla", "bdlnm"
    )
  )

  if (!is.null(random_effect)) {
    if (!is.character(random_effect) || length(random_effect) != 1L ||
        is.na(random_effect) || !nzchar(random_effect)) {
      stop("`random_effect` must be NULL or one non-empty column name.")
    }

    if (!random_effect %in% names(data)) {
      stop(
        "`random_effect` ('", random_effect,
        "') was not found in `data`."
      )
    }

    if (anyNA(data[[random_effect]])) {
      stop("`random_effect` cannot contain missing values.")
    }
  }

  # =========================================================
  # FAMILY HELPERS
  # =========================================================

  resolve_family_name <- function(family_input) {
    family_raw <- NULL

    if (is.character(family_input) && length(family_input) == 1L &&
        !is.na(family_input) && nzchar(family_input)) {
      family_raw <- family_input
    } else if (is.list(family_input) && !is.null(family_input$family) &&
               length(family_input$family) >= 1L) {
      family_raw <- as.character(family_input$family[[1]])
    } else if (inherits(family_input, "family") &&
               !is.null(family_input$family)) {
      family_raw <- as.character(family_input$family[[1]])
    } else {
      stop(
        "Unsupported `family` specification. Provide a supported family name ",
        "or an engine-compatible family object containing a `family` field."
      )
    }

    normalized <- tolower(trimws(family_raw))
    normalized <- gsub("[[:space:]-]+", "_", normalized)
    normalized <- gsub("[^a-z0-9_]", "", normalized)

    if (normalized %in% c("beta", "beta_family", "beta_proportion")) {
      return("beta")
    }

    if (normalized %in% c("binomial", "bernoulli")) {
      return("binomial")
    }

    if (normalized == "poisson") {
      return("poisson")
    }

    if (normalized == "gamma") {
      return("gamma")
    }

    if (normalized %in% c("gaussian", "normal")) {
      return("gaussian")
    }

    if (
      normalized %in% c(
        "negbin", "nbinom", "nbinom1", "nbinom2",
        "negative_binomial", "negative_binomial_1",
        "negative_binomial_2"
      ) || grepl("negative.*binomial", normalized)
    ) {
      return("negative_binomial")
    }

    if (normalized %in% c("ordinal", "cumulative")) {
      return("ordinal")
    }

    stop(
      "Unsupported family: '", family_raw, "'. Supported canonical families are: ",
      "beta, binomial, poisson, gamma, gaussian, negative_binomial, and ordinal."
    )
  }

  resolve_input_link <- function(family_input) {
    if (is.list(family_input) && !is.null(family_input$link) &&
        length(family_input$link) >= 1L) {
      link <- tolower(as.character(family_input$link[[1]]))
      if (!is.na(link) && nzchar(link)) return(link)
    }

    if (inherits(family_input, "family") && !is.null(family_input$link)) {
      link <- tolower(as.character(family_input$link[[1]]))
      if (!is.na(link) && nzchar(link)) return(link)
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
      ordinal = "logit",
      stop("Could not determine a default link for family '", family_name, "'.")
    )
  }

  resolve_engine_family <- function(
    family_input,
    family_name,
    model_engine,
    link_name
  ) {
    user_supplied_object <- !is.character(family_input)

    # Preserve compatible user-supplied objects for engines that accept them.
    if (
      user_supplied_object &&
      model_engine %in% c("glm", "glmmTMB", "gam", "gamm", "spamm", "brms")
    ) {
      return(family_input)
    }

    if (model_engine == "glm") {
      return(
        switch(
          family_name,
          gaussian = stats::gaussian(link = link_name),
          poisson = stats::poisson(link = link_name),
          gamma = stats::Gamma(link = link_name),
          binomial = stats::binomial(link = link_name),
          stop(
            "Family '", family_name,
            "' is not supported by `model_engine = 'glm'`."
          )
        )
      )
    }

    if (model_engine == "glmmTMB") {
      if (!requireNamespace("glmmTMB", quietly = TRUE)) {
        stop("Package 'glmmTMB' is required for `model_engine = 'glmmTMB'`.")
      }

      return(
        switch(
          family_name,
          beta = glmmTMB::beta_family(link = link_name),
          gaussian = stats::gaussian(link = link_name),
          poisson = stats::poisson(link = link_name),
          gamma = stats::Gamma(link = link_name),
          binomial = stats::binomial(link = link_name),
          negative_binomial = glmmTMB::nbinom2(link = link_name),
          stop(
            "Family '", family_name,
            "' is not supported by `model_engine = 'glmmTMB'`."
          )
        )
      )
    }

    if (model_engine %in% c("gam", "gamm")) {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required for GAM or GAMM models.")
      }

      if (family_name == "beta") {
        if (model_engine == "gamm") {
          stop(
            "Canonical Beta-family support is not harmonized for ",
            "`model_engine = 'gamm'`. Supply a verified engine-compatible ",
            "family object or use another supported engine."
          )
        }
        return(mgcv::betar(link = link_name))
      }

      return(
        switch(
          family_name,
          gaussian = stats::gaussian(link = link_name),
          poisson = stats::poisson(link = link_name),
          gamma = stats::Gamma(link = link_name),
          binomial = stats::binomial(link = link_name),
          negative_binomial = mgcv::nb(link = link_name),
          stop(
            "Family '", family_name,
            "' is not supported by `model_engine = '", model_engine, "'`."
          )
        )
      )
    }

    if (model_engine == "gls") {
      if (family_name != "gaussian") {
        stop("`model_engine = 'gls'` supports only `family = 'gaussian'`.")
      }
      return(stats::gaussian(link = "identity"))
    }

    if (model_engine == "spamm") {
      if (!requireNamespace("spaMM", quietly = TRUE)) {
        stop("Package 'spaMM' is required for `model_engine = 'spamm'`.")
      }

      return(
        switch(
          family_name,
          gaussian = stats::gaussian(link = link_name),
          poisson = stats::poisson(link = link_name),
          gamma = stats::Gamma(link = link_name),
          binomial = stats::binomial(link = link_name),
          stop(
            "Family '", family_name,
            "' is not mapped automatically for `model_engine = 'spamm'`. ",
            "Supply a verified spaMM-compatible family object."
          )
        )
      )
    }

    if (model_engine == "brms") {
      if (!requireNamespace("brms", quietly = TRUE)) {
        stop("Package 'brms' is required for `model_engine = 'brms'`.")
      }

      return(
        switch(
          family_name,
          beta = brms::Beta(link = link_name),
          gaussian = brms::gaussian(link = link_name),
          poisson = brms::poisson(link = link_name),
          gamma = brms::Gamma(link = link_name),
          binomial = brms::bernoulli(link = link_name),
          negative_binomial = brms::negbinomial(link = link_name),
          ordinal = brms::cumulative(link = link_name),
          stop(
            "Family '", family_name,
            "' is not supported by `model_engine = 'brms'`."
          )
        )
      )
    }

    if (model_engine == "inla") {
      return(
        switch(
          family_name,
          gaussian = "gaussian",
          poisson = "poisson",
          binomial = "binomial",
          gamma = "gamma",
          negative_binomial = "nbinomial",
          stop(
            "Family '", family_name,
            "' is not mapped for `model_engine = 'inla'`."
          )
        )
      )
    }

    if (model_engine == "bdlnm") {
      return(
        switch(
          family_name,
          gaussian = "gaussian",
          poisson = "poisson",
          binomial = "binomial",
          gamma = "gamma",
          negative_binomial = "nbinomial",
          stop(
            "Family '", family_name,
            "' is not mapped for `model_engine = 'bdlnm'`."
          )
        )
      )
    }

    stop("Unsupported `model_engine`: ", model_engine, ".")
  }

  family_name <- resolve_family_name(family)
  link_name <- resolve_input_link(family)
  if (is.null(link_name)) {
    link_name <- default_family_link(family_name)
  }

  if (family_name == "ordinal" && model_engine != "brms") {
    stop(
      "`family = 'ordinal'` is currently supported only with ",
      "`model_engine = 'brms'`."
    )
  }

  engine_family <- resolve_engine_family(
    family_input = family,
    family_name = family_name,
    model_engine = model_engine,
    link_name = link_name
  )

  # =========================================================
  # CHECK CONSISTENCY WITH PREPARE_RESPONSE()
  # =========================================================

  prepared_family_name <- attr(data, "response_family_name")

  if (!is.null(prepared_family_name)) {
    prepared_family_name <- resolve_family_name(prepared_family_name)

    if (!identical(prepared_family_name, family_name)) {
      stop(
        "The response was prepared for family '", prepared_family_name,
        "', but the model is being fitted with family '", family_name, "'."
      )
    }
  }

  # =========================================================
  # NORMALIZE AND VALIDATE EXPOSURE SPECIFICATION
  # =========================================================

  normalize_spec <- function(spec) {
    if (is.null(spec)) return(NULL)

    if (!is.list(spec) || is.null(names(spec)) ||
        anyNA(names(spec)) || any(!nzchar(names(spec)))) {
      stop("`epiexposure_spec` must be a named list.")
    }

    if (anyDuplicated(names(spec))) {
      stop("`epiexposure_spec` must contain unique variable names.")
    }

    for (nm in names(spec)) {
      current <- spec[[nm]]

      if (!is.list(current)) {
        stop("Specification for variable '", nm, "' must be a list.")
      }

      if (is.null(current$lag_max)) {
        stop("Missing `lag_max` for variable '", nm, "'.")
      }

      if (!is.numeric(current$lag_max) || !length(current$lag_max) ||
          any(!is.finite(current$lag_max)) || max(current$lag_max) < 0) {
        stop("Invalid `lag_max` for variable '", nm, "'.")
      }

      if (is.null(current$argvar) || !is.list(current$argvar)) {
        stop("Missing or invalid `argvar` for variable '", nm, "'.")
      }

      if (is.null(current$arglag) || !is.list(current$arglag)) {
        stop("Missing or invalid `arglag` for variable '", nm, "'.")
      }

      current$lag_max <- as.integer(max(current$lag_max))
      spec[[nm]] <- current
    }

    spec
  }

  epiexposure_spec <- normalize_spec(epiexposure_spec)

  # =========================================================
  # ORDER CROSS-BASIS COLUMNS
  # =========================================================

  sort_cb_cols <- function(cols) {
    if (!length(cols)) return(cols)

    vars <- sub("^cb_", "", cols)
    vars <- sub("_[0-9]+$", "", vars)

    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", cols)))
    missing_idx <- is.na(idx)
    idx[missing_idx] <- seq_along(cols)[missing_idx]

    cols[order(vars, idx)]
  }

  cb_cols <- sort_cb_cols(grep("^cb_", names(data), value = TRUE))

  if (!length(cb_cols) && model_engine != "bdlnm") {
    stop("No `cb_*` columns were found in `data`.")
  }

  infer_vars <- function(columns) {
    variables <- sub("^cb_", "", columns)
    variables <- sub("_[0-9]+$", "", variables)
    unique(variables)
  }

  vars_inferred <- if (length(cb_cols)) {
    infer_vars(cb_cols)
  } else if (
    model_engine == "bdlnm" &&
    !is.null(basis_objects) &&
    is.list(basis_objects) &&
    !is.null(names(basis_objects))
  ) {
    names(basis_objects)
  } else {
    character(0)
  }

  if (!length(vars_inferred)) {
    stop("Could not determine the exposure variables used by the fitted model.")
  }

  if (!is.null(epiexposure_spec)) {
    missing_spec_vars <- setdiff(vars_inferred, names(epiexposure_spec))
    if (length(missing_spec_vars)) {
      stop(
        "`epiexposure_spec` is missing variables used by the model: ",
        paste(missing_spec_vars, collapse = ", "), "."
      )
    }
  }

  # =========================================================
  # VALIDATE BDLNM BASIS OBJECTS
  # =========================================================

  if (model_engine == "bdlnm") {
    if (is.null(basis_objects) || !is.list(basis_objects) ||
        is.null(names(basis_objects)) || !length(basis_objects) ||
        anyNA(names(basis_objects)) || any(!nzchar(names(basis_objects))) ||
        anyDuplicated(names(basis_objects))) {
      stop(
        "For `model_engine = 'bdlnm'`, `basis_objects` must be a non-empty ",
        "named list with unique names."
      )
    }
  }

  # =========================================================
  # ATTACH STANDARDIZED METADATA
  # =========================================================

  attach_epiexposure_meta <- function(model_obj) {
    data_template <- data[1, , drop = FALSE]

    attr(model_obj, "epiexposure_engine") <- model_engine
    attr(model_obj, "epiexposure_family_input") <- family
    attr(model_obj, "epiexposure_family") <- engine_family
    attr(model_obj, "epiexposure_family_name") <- family_name
    attr(model_obj, "epiexposure_link") <- link_name
    attr(model_obj, "epiexposure_cb_cols") <- cb_cols
    attr(model_obj, "epiexposure_vars") <- vars_inferred
    attr(model_obj, "epiexposure_data_template") <- data_template
    attr(model_obj, "epiexposure_id_col") <- random_effect
    attr(model_obj, "epiexposure_spec") <- epiexposure_spec
    attr(model_obj, "epiexposure_basis_objects") <- basis_objects

    model_obj
  }

  # =========================================================
  # FORMULA CONSTRUCTION
  # =========================================================

  rhs <- if (length(cb_cols)) paste(cb_cols, collapse = " + ") else ""

  fml_fixed <- if (nzchar(rhs)) {
    paste0("y_model ~ 1 + ", rhs)
  } else {
    "y_model ~ 1"
  }

  fml <- if (
    !is.null(random_effect) &&
    model_engine %in% c("glmmTMB", "gamm", "spamm", "brms")
  ) {
    stats::as.formula(paste0(fml_fixed, " + (1|", random_effect, ")"))
  } else {
    stats::as.formula(fml_fixed)
  }

  # =========================================================
  # FREQUENTIST ENGINES
  # =========================================================

  if (model_engine == "glm") {
    mod <- stats::glm(
      formula = fml,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "glmmTMB") {
    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package 'glmmTMB' is required for `model_engine = 'glmmTMB'`.")
    }

    mod <- glmmTMB::glmmTMB(
      formula = fml,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "gam") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for `model_engine = 'gam'`.")
    }

    mod <- mgcv::gam(
      formula = fml,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "gamm") {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
      stop("Package 'mgcv' is required for `model_engine = 'gamm'`.")
    }

    mod <- mgcv::gamm(
      formula = fml,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "gls") {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop("Package 'nlme' is required for `model_engine = 'gls'`.")
    }

    mod <- nlme::gls(
      model = fml,
      data = data,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  if (model_engine == "spamm") {
    if (!requireNamespace("spaMM", quietly = TRUE)) {
      stop("Package 'spaMM' is required for `model_engine = 'spamm'`.")
    }

    mod <- spaMM::fitme(
      formula = fml,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  # =========================================================
  # BAYESIAN ENGINE: BRMS
  # =========================================================

  if (model_engine == "brms") {
    if (!requireNamespace("brms", quietly = TRUE)) {
      stop("Package 'brms' is required for `model_engine = 'brms'`.")
    }

    mod <- brms::brm(
      formula = fml,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  # =========================================================
  # BAYESIAN ENGINE: INLA
  # =========================================================

  if (model_engine == "inla") {
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("Package 'INLA' is required for `model_engine = 'inla'`.")
    }

    fml_inla <- if (!is.null(random_effect)) {
      stats::as.formula(
        paste0(fml_fixed, " + f(", random_effect, ", model = 'iid')")
      )
    } else {
      stats::as.formula(fml_fixed)
    }

    mod <- INLA::inla(
      formula = fml_inla,
      data = data,
      family = engine_family,
      ...
    )
    return(attach_epiexposure_meta(mod))
  }

  # =========================================================
  # BAYESIAN DLNM ENGINE: BDLNM
  # =========================================================

  if (model_engine == "bdlnm") {
    if (!requireNamespace("bdlnm", quietly = TRUE)) {
      stop("Package 'bdlnm' is required for `model_engine = 'bdlnm'`.")
    }

    basis_names <- names(basis_objects)
    rhs_bdlnm <- paste(basis_names, collapse = " + ")

    fml_bdlnm <- if (!is.null(random_effect)) {
      stats::as.formula(
        paste0(
          "y_model ~ 1 + ", rhs_bdlnm,
          " + f(", random_effect, ", model = 'iid')"
        )
      )
    } else {
      stats::as.formula(paste0("y_model ~ 1 + ", rhs_bdlnm))
    }

    eval_env <- new.env(parent = parent.frame())
    for (nm in basis_names) {
      assign(nm, basis_objects[[nm]], envir = eval_env)
    }
    environment(fml_bdlnm) <- eval_env

    mod <- bdlnm::bdlnm(
      formula = fml_bdlnm,
      data = data,
      family = engine_family,
      ...
    )

    return(attach_epiexposure_meta(mod))
  }

  stop("Unsupported `model_engine`.")
}
