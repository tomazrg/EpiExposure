#' Compute lag-specific decomposition of Exposure Cumulative Impact (ECI)
#'
#' Decomposes the model-weighted Exposure Cumulative Impact (ECI) of one or
#' more exposure profiles into exact lag-specific contributions on the
#' linear-predictor scale.
#'
#' Each evaluated exposure profile is compared with an explicit joint reference
#' profile. The focal exposure is allowed to vary through time while every other
#' fitted exposure remains fixed at its reference value. Lag-specific
#' contributions are obtained by changing one chronological exposure position
#' at a time relative to that same reference profile and reconstructing the
#' corresponding fitted cross-basis through the standard EpiExposure prediction
#' pipeline.
#'
#' The resulting lag contributions are additive: their sum is required to equal
#' the centered cumulative ECI for the complete focal exposure trajectory.
#' `compute_ecilag()` therefore performs an exact decomposition of the fitted
#' DLNM contrast. It does not calculate numerical derivatives, marginal effects,
#' or local sensitivity measures.
#'
#' @param profile Optional exposure-profile input. Supply exactly one of
#'   `profile` or `data`.
#'
#'   For a single focal exposure, `profile` may be one finite numeric vector.
#'   For multiple focal exposures, supply a list containing one numeric vector
#'   per variable. A named list is recommended; if unnamed, its order must match
#'   `vars`.
#'
#'   Exposure profiles must be chronological from the oldest observation to the
#'   most recent observation. Under the EpiExposure exact-profile contract, every
#'   supplied profile must contain exactly `max_lag + 1` observations. The first
#'   value corresponds to lag `max_lag` and the final value to lag 0. Profiles
#'   are never truncated, padded, or silently realigned.
#'
#' @param data Optional non-empty long-format data frame containing observed
#'   exposure profiles. Supply exactly one of `profile` or `data`.
#'
#'   `data` must contain `group`, `time`, and every exposure requested in `vars`.
#'   Each evaluated group must contain exactly `max_lag + 1` observations, with
#'   unique and equally spaced time values. All evaluated groups must use the
#'   same temporal spacing.
#'
#'   When a method-based reference is requested, `data` must additionally contain
#'   every exposure included in the fitted model because the complete joint
#'   reference profile is constructed from those variables.
#'
#' @param group Character scalar naming the independent-profile grouping column
#'   in `data`. Required when `data` is supplied.
#'
#' @param time Character scalar naming the chronological numeric time column in
#'   `data`. Default is `"time"`.
#'
#' @param group_level Optional grouping value or vector of grouping values to
#'   evaluate. If `NULL`, all groups are evaluated in first-occurrence order.
#'
#' @param fit Fitted EpiExposure model returned by the current
#'   `fit_epidlnm()` implementation. The fitted model must satisfy the current
#'   population-level expected-response, central-parameter, uncertainty, and
#'   exact-profile metadata contracts.
#'
#' @param vars Optional character vector naming fitted exposure variables to
#'   decompose.
#'
#'   If the fitted model contains only one exposure, `vars = NULL` uses that
#'   exposure automatically. For multivariable fits, `vars` must be supplied
#'   explicitly.
#'
#'   Each requested variable is decomposed separately. When one focal variable
#'   is evaluated, all other fitted exposures remain fixed at their joint
#'   reference values; contributions from different exposure variables are not
#'   pooled into a single lag decomposition.
#'
#' @param ref Reference exposure specification.
#'
#'   Two forms are supported:
#'
#'   \itemize{
#'     \item Method-based:
#'       `list(method = "median", value = NULL)`,
#'       `list(method = "percentile", value = p)`, or
#'       `list(method = "fixed", value = x)`.
#'     \item Exposure-specific: a named list containing exactly one finite
#'       reference value for every fitted exposure, for example
#'       `list(tmean = 25, rain = 5, wetness = 70)`.
#'   }
#'
#'   Method-based references require `data`. `"median"` uses each fitted
#'   exposure's median, `"percentile"` applies the same percentile probability
#'   to each fitted exposure, and `"fixed"` applies the same numeric value to
#'   every fitted exposure and should therefore be used cautiously when
#'   exposures have different units.
#'
#'   The default is `list(method = "median", value = NULL)`.
#'
#' @param absolute Logical. If `TRUE`, add the absolute lag contribution and,
#'   when uncertainty is requested, its uncertainty summary. The signed
#'   `ECI_weighted` column is always retained.
#'
#'   `ECI_percent` is always based on the absolute magnitude of each
#'   lag contribution:
#'
#'   \deqn{
#'     100
#'     \frac{|C_l|}
#'          {\sum_j |C_j|}.
#'   }
#'
#'   Consequently, positive and negative lag effects do not cancel when relative
#'   lag importance is calculated. If the total absolute contribution is zero,
#'   percentage contributions are undefined and returned as `NA`.
#'
#' @param uncertainty Logical. If `FALSE`, lag contributions are calculated from
#'   the harmonized central fixed/population parameter estimate. If `TRUE`,
#'   joint parameter uncertainty is propagated draw by draw through the complete
#'   lag decomposition.
#'
#' @param output Character. `"summary"` returns deterministic lag contributions
#'   when `uncertainty = FALSE`, or medians, standard deviations, and empirical
#'   intervals when `uncertainty = TRUE`.
#'
#'   `"samples"` additionally returns the lag-specific and total weighted ECI
#'   values for every parameter draw and requires `uncertainty = TRUE`.
#'
#' @param n_samples Positive integer number of parameter draws used when
#'   `uncertainty = TRUE`. At least two are required. Default is 1000.
#'
#' @param interval_probs Numeric vector of length two specifying the empirical
#'   uncertainty interval. The default `c(0.025, 0.975)` gives a 95 percent
#'   interval.
#'
#' @param seed Optional finite integer used for parameter-draw generation or
#'   sampling. When supplied, the caller's random-number state is restored when
#'   the function exits.
#'
#' @param extrapolation Character controlling exposure values outside the range
#'   recorded by the fitted cross-basis: `"error"` (default), `"warn"`, or
#'   `"allow"`. The fitted basis specification is never re-estimated from the
#'   ECI-lag input data.
#'
#' @return A list containing cumulative and lag-specific ECI results.
#'
#'   For a single exposure profile, the principal components are:
#'
#'   \describe{
#'     \item{`ECI_raw`}{Descriptive sum of the focal exposure profile.}
#'     \item{`ECI_raw_centered`}{Descriptive sum of exposure deviations from the
#'       focal reference value.}
#'     \item{`ECI_weighted`}{Centered cumulative exposure contribution on the
#'       linear-predictor scale.}
#'     \item{`reference_value`}{Reference value for the focal exposure.}
#'     \item{`max_lag`}{Common fitted maximum lag.}
#'     \item{`n_exposure_values`}{Number of observations in the profile; always
#'       `max_lag + 1`.}
#'     \item{`by_lag`}{Data frame containing the exact contribution assigned to
#'       each lag.}
#'   }
#'
#'   With `uncertainty = TRUE`, `ECI_weighted` is summarized by its median,
#'   standard deviation, and empirical lower and upper quantiles. The `by_lag`
#'   table similarly contains draw-based summaries for signed contributions,
#'   percentage contributions, and, when `absolute = TRUE`, absolute
#'   contributions.
#'
#'   With `output = "samples"`, `by_lag_samples` and
#'   `ECI_weighted_samples` contain the corresponding draw-level quantities.
#'
#'   When `data` or multiple focal exposure profiles are evaluated, cumulative
#'   results are returned in `eci_summary` and lag-specific results in
#'   `by_lag`, with the grouping column included when applicable.
#'
#'   The result also stores attributes describing the exact ECI-lag
#'   decomposition, common fitted maximum lag, profile length, reference values,
#'   and uncertainty setting.
#'
#' @details
#' ## Centered cumulative impact
#'
#' Let \eqn{X_{ref}} denote the fixed/population design generated when every
#' fitted exposure follows its reference profile. For focal exposure \eqn{v},
#' let \eqn{X_v} denote the design in which only that exposure follows the
#' evaluated complete trajectory. The model-weighted cumulative ECI is
#'
#' \deqn{
#'   ECI_{weighted,v}
#'   =
#'   (X_v-X_{ref})\beta.
#' }
#'
#' This quantity is a centered contrast on the linear-predictor scale. The
#' intercept and all non-focal exposure terms cancel because they are identical
#' in the target and reference designs.
#'
#' ## Exact lag decomposition
#'
#' For each chronological position corresponding to lag \eqn{l}, the function
#' creates an isolated target profile in which only that exposure value differs
#' from the reference trajectory. The resulting design contrast is
#'
#' \deqn{
#'   D_l = X_l-X_{ref}.
#' }
#'
#' Its lag-specific contribution is
#'
#' \deqn{
#'   C_l = D_l\beta.
#' }
#'
#' The implementation explicitly verifies
#'
#' \deqn{
#'   \sum_{l=0}^{L} D_l
#'   =
#'   X_v-X_{ref}
#' }
#'
#' and therefore
#'
#' \deqn{
#'   \sum_{l=0}^{L} C_l
#'   =
#'   ECI_{weighted,v}.
#' }
#'
#' If this additivity condition is not satisfied within numerical tolerance, the
#' function stops rather than returning an approximate decomposition.
#'
#' This is fundamentally different from a derivative-based sensitivity analysis.
#' `compute_ecilag()` assigns the complete centered fitted contrast to exact lag
#' positions; it does not estimate \eqn{\partial\eta/\partial x_l}.
#'
#' ## Chronological and lag indexing
#'
#' Input profiles are chronological. For fitted maximum lag \eqn{L}, a profile
#' of length \eqn{L+1} is interpreted as
#'
#' \deqn{
#'   (x_L, x_{L-1}, \ldots, x_1, x_0),
#' }
#'
#' where the first supplied observation corresponds to lag \eqn{L} and the last
#' to lag 0. Returned `by_lag` tables are ordered by increasing lag.
#'
#' ## Joint reference profile
#'
#' Even when only one exposure is decomposed, the reference condition is defined
#' jointly for all fitted exposures. This keeps the decomposition consistent with
#' the fitted multivariable model. Non-focal exposures remain fixed at their own
#' reference values throughout the entire lag window.
#'
#' ## Parameter uncertainty
#'
#' With `uncertainty = TRUE`, the same joint fixed/population parameter draw is
#' used for all lag contributions within a trajectory. For draw \eqn{s},
#'
#' \deqn{
#'   C_l^{(s)} = D_l\beta^{(s)}
#' }
#'
#' and
#'
#' \deqn{
#'   ECI_{weighted}^{(s)}
#'   =
#'   \sum_l C_l^{(s)}.
#' }
#'
#' Additivity is verified for every draw. Summary statistics are calculated only
#' after the draw-specific contributions and percentage contributions have been
#' obtained.
#'
#' Frequentist engines use joint draws derived from the fitted fixed-effect
#' covariance matrix. Bayesian engines use joint posterior or
#' approximate-posterior fixed/population coefficient draws through the
#' centralized EpiExposure prediction helpers.
#'
#' No residual, observation, process, dispersion, or posterior-predictive noise
#' is added. Uncertainty therefore represents parameter uncertainty in the fitted
#' exposure-lag association.
#'
#' @seealso
#' `compute_eci()`, `summarise_effects()`, `predict_outcomes()`
#'
#' @export
#'
#' @export
compute_ecilag <- function(
    profile = NULL,
    data = NULL,
    group = NULL,
    time = "time",
    group_level = NULL,
    fit,
    vars = NULL,
    ref = list(method = "median", value = NULL),
    absolute = TRUE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    interval_probs = c(0.025, 0.975),
    seed = NULL,
    extrapolation = c("error", "warn", "allow")
) {
  output <- match.arg(output)
  extrapolation <- match.arg(extrapolation)
  stopf <- function(...) stop(..., call. = FALSE)
  scalar_name <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  scalar_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  scalar_int <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x == as.integer(x)
  nearly_equal <- function(x, y, tol = 1e-8) {
    abs(x - y) <= tol * pmax(1, abs(x), abs(y))
  }
  summarize_finite <- function(x) {
    if (!is.numeric(x) || !length(x) || anyNA(x) || any(!is.finite(x)))
      stopf("Internal ECI-lag draw values are not finite.")
    c(estimate = stats::median(x), sd = stats::sd(x),
      lower = stats::quantile(x, interval_probs[1L], names = FALSE),
      upper = stats::quantile(x, interval_probs[2L], names = FALSE))
  }
  summarize_defined <- function(x) {
    z <- x[is.finite(x)]
    if (!length(z)) return(c(estimate=NA, sd=NA, lower=NA, upper=NA, n_defined=0))
    c(estimate=stats::median(z), sd=if(length(z)>1L) stats::sd(z) else NA,
      lower=stats::quantile(z, interval_probs[1L], names=FALSE),
      upper=stats::quantile(z, interval_probs[2L], names=FALSE), n_defined=length(z))
  }

  if (missing(fit) || is.null(fit)) stopf("`fit` must be provided.")
  if (!scalar_flag(absolute) || !scalar_flag(uncertainty))
    stopf("`absolute` and `uncertainty` must be TRUE or FALSE.")
  if (!uncertainty && output == "samples") stopf("`output = 'samples'` requires `uncertainty = TRUE`.")
  if (!scalar_int(n_samples) || n_samples < 1L || (uncertainty && n_samples < 2L))
    stopf("`n_samples` must be a positive integer and at least 2 with uncertainty.")
  n_samples <- as.integer(n_samples)
  if (!is.numeric(interval_probs) || length(interval_probs)!=2L || anyNA(interval_probs) ||
      any(!is.finite(interval_probs)) || interval_probs[1L]<=0 || interval_probs[2L]>=1 ||
      interval_probs[1L]>=interval_probs[2L]) stopf("Invalid `interval_probs`.")
  if (!is.null(seed) && !scalar_int(seed)) stopf("`seed` must be NULL or one finite integer.")

  using_profile <- !is.null(profile); using_data <- !is.null(data)
  if (identical(using_profile, using_data)) stopf("Supply exactly one of `profile` or `data`.")

  metadata <- .get_epiexposure_metadata(fit)
  if (!identical(metadata$prediction_level, "population") ||
      !identical(metadata$prediction_estimand, "expected_response") ||
      !identical(metadata$point_prediction_contract, "central_expected_response") ||
      !identical(metadata$uncertainty_contract, "draw_by_draw_median_quantiles"))
    stopf("The fitted model does not satisfy the current EpiExposure prediction contract.")

  fitted_lags <- vapply(metadata$vars, function(v) as.integer(metadata$spec[[v]]$max_lag), integer(1))
  if (anyNA(fitted_lags) || length(unique(fitted_lags)) != 1L)
    stopf("All fitted exposures must have the same valid `max_lag`.")
  max_lag <- unname(fitted_lags[1L]); history_length <- max_lag + 1L

  if (is.null(vars)) {
    if (length(metadata$vars) != 1L) stopf("Supply `vars` explicitly for a multivariable fit.")
    variables <- metadata$vars
  } else {
    if (!is.character(vars) || !length(vars) || anyNA(vars) || any(!nzchar(vars)) || anyDuplicated(vars))
      stopf("`vars` must contain unique fitted exposure names.")
    unknown <- setdiff(vars, metadata$vars)
    if (length(unknown)) stopf("Unknown fitted exposure(s): ", paste(unknown, collapse=", "), ".")
    variables <- vars
  }

  # Validate input and create chronological profile records (oldest to newest).
  records <- list()
  if (using_profile) {
    if (is.numeric(profile) && is.null(dim(profile))) {
      if (length(variables)!=1L) stopf("A numeric `profile` requires exactly one `vars`.")
      profile <- stats::setNames(list(as.numeric(profile)), variables)
    }
    if (!is.list(profile)) stopf("`profile` must be a numeric vector or list.")
    if (is.null(names(profile))) {
      if (length(profile)!=length(variables)) stopf("Unnamed `profile` must follow `vars` exactly.")
      names(profile) <- variables
    }
    if (!setequal(names(profile), variables)) stopf("`profile` names must match `vars` exactly.")
    profile <- profile[variables]
    for (v in variables) {
      x <- profile[[v]]
      if (!is.numeric(x) || length(x)!=history_length || anyNA(x) || any(!is.finite(x)))
        stopf("Profile for '", v, "' must contain exactly ", history_length, " finite values.")
      records[[length(records)+1L]] <- list(variable=v, chronology=as.numeric(x), group_value=NULL)
    }
  } else {
    if (!is.data.frame(data) || !nrow(data)) stopf("`data` must be a non-empty data.frame.")
    if (!scalar_name(group) || !group %in% names(data)) stopf("`group` must name a column in `data`.")
    if (!scalar_name(time) || !time %in% names(data)) stopf("`time` must name a column in `data`.")
    needed <- unique(c(group, time, variables))
    if (length(setdiff(needed, names(data)))) stopf("Required columns are missing from `data`.")
    .epix_validate_regular_time(data, group_index=as.character(data[[group]]), time_col=time)
    counts <- table(as.character(data[[group]]))
    if (any(counts != history_length)) stopf("Every group must contain exactly ", history_length, " observations.")
    available <- unique(data[[group]])
    levels_use <- if (is.null(group_level)) available else {
      idx <- match(as.character(group_level), as.character(available))
      if (anyNA(idx) || anyDuplicated(idx)) stopf("Invalid or duplicated `group_level`.")
      available[idx]
    }
    for (g in levels_use) {
      d <- data[as.character(data[[group]]) == as.character(g), , drop=FALSE]
      d <- d[order(d[[time]]), , drop=FALSE]
      for (v in variables) {
        x <- d[[v]]
        if (!is.numeric(x) || anyNA(x) || any(!is.finite(x))) stopf("Exposure '",v,"' must be finite numeric.")
        records[[length(records)+1L]] <- list(variable=v, chronology=as.numeric(x), group_value=g)
      }
    }
  }

  # Resolve joint reference values.
  method_ref <- is.list(ref) && !is.null(names(ref)) && setequal(names(ref), c("method","value"))
  if (method_ref) {
    if (using_profile) stopf("Method-based `ref` requires `data`.")
    method <- ref$method
    if (!scalar_name(method) || !method %in% c("median","percentile","fixed")) stopf("Invalid reference method.")
    missing_ref_data <- setdiff(metadata$vars, names(data))
    if (length(missing_ref_data)) stopf("Reference data are missing fitted exposures: ", paste(missing_ref_data, collapse=", "), ".")
    if (method == "median") {
      if (!is.null(ref$value)) stopf("`ref$value` must be NULL for median reference.")
      ref_values <- vapply(metadata$vars, function(v) stats::median(data[[v]]), numeric(1))
    } else if (method == "percentile") {
      p <- ref$value
      if (!is.numeric(p) || length(p)!=1L || is.na(p) || p<0 || p>1) stopf("Invalid reference percentile.")
      ref_values <- vapply(metadata$vars, function(v) stats::quantile(data[[v]], p, names=FALSE), numeric(1))
    } else {
      x <- ref$value
      if (!is.numeric(x) || length(x)!=1L || is.na(x) || !is.finite(x)) stopf("Invalid fixed reference.")
      ref_values <- stats::setNames(rep(as.numeric(x), length(metadata$vars)), metadata$vars)
    }
  } else {
    if (!is.list(ref) || is.null(names(ref)) || !setequal(names(ref), metadata$vars))
      stopf("Exposure-specific `ref` names must match all fitted exposures exactly.")
    ref <- ref[metadata$vars]
    ref_values <- vapply(ref, function(x) {
      if (!is.numeric(x) || length(x)!=1L || is.na(x) || !is.finite(x)) stopf("Each reference must be one finite value.")
      as.numeric(x)
    }, numeric(1))
  }
  names(ref_values) <- metadata$vars

  # Extrapolation checks for reference and evaluated exposure values.
  for (v in metadata$vars) {
    def <- .epix_basis_definition(metadata, v)
    vals <- ref_values[[v]]
    if (v %in% variables) {
      actual <- unlist(
        lapply(
          records[vapply(records, function(z) identical(z$variable, v), logical(1))],
          `[[`,
          "chronology"
        ),
        use.names = FALSE
      )
      vals <- unique(c(vals, actual))
    }
    .epix_handle_extrapolation(vals, def$exposure_range, v, extrapolation)
  }

  reference_profiles <- stats::setNames(lapply(metadata$vars, function(v) rep(ref_values[[v]], history_length)), metadata$vars)
  build_design <- function(profiles, context) {
    nd <- .build_newdata_basis(fit=fit, profiles=profiles, extrapolation="allow")
    X <- .epix_standard_fixed_design(nd, metadata)
    if (!is.matrix(X) || nrow(X)!=1L || anyNA(X) || any(!is.finite(X))) stopf("Invalid fixed design for ", context, ".")
    X
  }
  X_reference <- build_design(reference_profiles, "joint reference")

  isolated_design <- function(variable, position, value, variable_columns) {
    profiles <- reference_profiles
    profiles[[variable]][position] <- value
    X <- build_design(profiles, paste0("exposure '",variable,"', position ",position))
    if (!identical(colnames(X), colnames(X_reference))) stopf("Fixed-design columns changed during lag isolation.")
    z <- as.numeric(X[1L,] - X_reference[1L,]); names(z) <- colnames(X_reference)
    other <- setdiff(names(z), variable_columns)
    tol <- 1e-8 * pmax(1, abs(z[other]))
    if (length(other) && any(abs(z[other]) > tol)) stopf("Isolated-lag construction changed non-focal terms for '",variable,"'.")
    z[variable_columns]
  }

  central <- .extract_central_parameters(fit)
  if (!identical(names(central), colnames(X_reference))) stopf("Central parameters do not align with fixed design.")
  draws <- NULL
  if (uncertainty) {
    if (!is.null(seed)) {
      had_seed <- exists(".Random.seed", .GlobalEnv, inherits=FALSE)
      if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
      on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else if (exists(".Random.seed",.GlobalEnv,FALSE)) rm(".Random.seed",envir=.GlobalEnv), add=TRUE)
      set.seed(as.integer(seed))
    }
    draws <- .extract_parameter_draws(fit, n_samples=n_samples)
    if (!is.matrix(draws) || nrow(draws)!=n_samples || !identical(colnames(draws),colnames(X_reference))) stopf("Invalid parameter draws.")
  }

  process_one <- function(rec) {
    v <- rec$variable; chronology <- rec$chronology
    vc <- .epix_cb_cols_for_var(metadata$cb_cols, v)
    if (!length(vc) || any(!vc %in% colnames(X_reference))) stopf("Invalid cross-basis mapping for '",v,"'.")
    lag_position <- seq.int(max_lag, 0L)
    C <- matrix(
      NA_real_,
      history_length,
      length(vc),
      dimnames = list(NULL, vc)
    )
    for (pos in seq_len(history_length)) {
      C[pos, ] <- isolated_design(
        v,
        pos,
        chronology[pos],
        vc
      )
    }
    full_profiles <- reference_profiles; full_profiles[[v]] <- chronology
    X_full <- build_design(full_profiles, paste0("complete profile '",v,"'"))
    full <- as.numeric(X_full[1L,]-X_reference[1L,]); names(full) <- colnames(X_reference)
    summed <- colSums(C); target <- full[vc]; diff <- summed-target
    tol <- 1e-8*pmax(1,abs(summed),abs(target))
    if (any(abs(diff)>tol)) stopf("Exact lag-specific designs do not sum to the centered trajectory for '",v,
                                  "'. Maximum absolute difference: ",format(max(abs(diff)),scientific=TRUE),".")
    beta <- central[vc]
    cp <- as.numeric(C %*% beta)
    weighted <- as.numeric(sum(full*central))
    if (!nearly_equal(sum(cp),weighted)) stopf("Central lag contributions do not sum to ECI_weighted for '",v,"'.")
    pct <- if (sum(abs(cp))>0) 100*abs(cp)/sum(abs(cp)) else rep(NA_real_,history_length)
    base <- data.frame(
      var = v,
      lag = lag_position,
      exposure = chronology,
      reference_value = ref_values[[v]],
      exposure_minus_reference = chronology - ref_values[[v]],
      ECI_weighted = cp,
      ECI_percent = pct,
      stringsAsFactors = FALSE
    )
    if (absolute) {
      base$ECI_absolute <- abs(cp)
    }
    sample_table <- weighted_table <- NULL
    summary <- data.frame(var=v,reference_value=ref_values[[v]],max_lag=max_lag,n_exposure_values=history_length,
                          ECI_raw=sum(chronology),ECI_raw_centered=sum(chronology-ref_values[[v]]),ECI_weighted=weighted,stringsAsFactors=FALSE)
    if (uncertainty) {
      CD <- draws[, vc, drop = FALSE] %*% t(C)
      WD <- as.numeric(draws %*% full)
      if (any(abs(rowSums(CD)-WD)>1e-8*pmax(1,abs(rowSums(CD)),abs(WD)))) stopf("Draw-level additivity failed for '",v,"'.")
      den <- rowSums(abs(CD)); PD <- matrix(NA_real_,n_samples,history_length); ok <- den>0
      PD[ok,] <- 100*abs(CD[ok,,drop=FALSE])/den[ok]
      cs <- t(vapply(
        seq_len(history_length),
        function(j) summarize_finite(CD[, j]),
        numeric(4)
      ))
      ps <- t(vapply(
        seq_len(history_length),
        function(j) summarize_defined(PD[, j]),
        numeric(5)
      ))
      base$ECI_weighted<-cs[,1]; base$ECI_weighted_sd<-cs[,2]; base$ECI_weighted_lower<-cs[,3]; base$ECI_weighted_upper<-cs[,4]
      base$ECI_percent<-ps[,1]; base$ECI_percent_sd<-ps[,2]; base$ECI_percent_lower<-ps[,3]; base$ECI_percent_upper<-ps[,4]; base$ECI_percent_n<-as.integer(ps[,5])
      if (absolute) {
        acs <- t(vapply(
          seq_len(history_length),
          function(j) summarize_finite(abs(CD[, j])),
          numeric(4)
        ))
        base$ECI_absolute <- acs[, 1]
        base$ECI_absolute_sd <- acs[, 2]
        base$ECI_absolute_lower <- acs[, 3]
        base$ECI_absolute_upper <- acs[, 4]
      }
      ws <- summarize_finite(WD)
      summary$ECI_weighted<-ws[1]; summary$ECI_weighted_sd<-ws[2]; summary$ECI_weighted_lower<-ws[3]; summary$ECI_weighted_upper<-ws[4]
      if (output=="samples") {
        sample_table <- do.call(rbind,lapply(seq_len(n_samples),function(s) {
          z <- data.frame(
            sample = s,
            var = v,
            lag = lag_position,
            exposure = chronology,
            reference_value = ref_values[[v]],
            exposure_minus_reference = chronology - ref_values[[v]],
            ECI_weighted = CD[s, ],
            ECI_percent = PD[s, ],
            stringsAsFactors = FALSE
          )
          if (absolute) {
            z$ECI_absolute <- abs(z$ECI_weighted)
          }
          z
        }))
        weighted_table<-data.frame(sample=seq_len(n_samples),ECI_weighted=WD)
      }
    }
    base<-base[order(base$lag),,drop=FALSE]; rownames(base)<-NULL
    list(summary=summary,by_lag=base,by_lag_samples=sample_table,ECI_weighted_samples=weighted_table)
  }

  results <- lapply(records, process_one)
  if (using_data) for (i in seq_along(results)) {
    g <- records[[i]]$group_value
    for (nm in c("summary","by_lag","by_lag_samples","ECI_weighted_samples")) if (!is.null(results[[i]][[nm]])) {
      results[[i]][[nm]][[group]] <- g
      results[[i]][[nm]] <- results[[i]][[nm]][,c(group,setdiff(names(results[[i]][[nm]]),group)),drop=FALSE]
    }
  }
  if (length(results)==1L && !using_data) {
    z<-results[[1L]]; s<-z$summary
    out<-list(ECI_raw=s$ECI_raw[[1]],ECI_raw_centered=s$ECI_raw_centered[[1]],ECI_weighted=s$ECI_weighted[[1]],
              reference_value=s$reference_value[[1]],max_lag=s$max_lag[[1]],n_exposure_values=s$n_exposure_values[[1]],by_lag=z$by_lag)
    if (uncertainty) { out$ECI_weighted_sd<-s$ECI_weighted_sd[[1]]; out$ECI_weighted_lower<-s$ECI_weighted_lower[[1]]; out$ECI_weighted_upper<-s$ECI_weighted_upper[[1]] }
    if (uncertainty && output=="samples") { out$by_lag_samples<-z$by_lag_samples; out$ECI_weighted_samples<-z$ECI_weighted_samples }
  } else {
    out<-list(eci_summary=do.call(rbind,lapply(results,`[[`,"summary")),by_lag=do.call(rbind,lapply(results,`[[`,"by_lag")))
    if (uncertainty && output=="samples") {
      out$by_lag_samples<-do.call(rbind,lapply(results,`[[`,"by_lag_samples"))
      out$ECI_weighted_samples<-do.call(rbind,lapply(seq_along(results),function(i){z<-results[[i]]$ECI_weighted_samples;z$var<-records[[i]]$variable;z}))
    }
  }
  attr(out,"epiexposure_ecilag_contract") <- "exact_lag_contributions_sum_to_centered_eci_link"
  attr(out,"epiexposure_ecilag_scale") <- "link"
  attr(out,"epiexposure_ecilag_reference_values") <- ref_values
  attr(out,"epiexposure_ecilag_max_lag") <- max_lag
  attr(out,"epiexposure_ecilag_history_length") <- history_length
  attr(out,"epiexposure_uncertainty") <- uncertainty
  out
}
