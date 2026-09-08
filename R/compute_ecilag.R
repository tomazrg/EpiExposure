#' Compute lag-specific decomposition of Exposure Cumulative Impact (ECI)
#'
#' Decomposes a complete exposure trajectory into exact lag-specific centered
#' contributions on the linear-predictor scale. The complete trajectory and
#' every isolated-lag trajectory are constructed through the same EpiExposure
#' basis pipeline.
#'
#' This function performs an additive lag decomposition of centered ECI. It does
#' not calculate numerical derivatives or local sensitivity measures.
#'
#' @export
compute_ecilag <- function(
    profile = NULL,
    data = NULL,
    group = NULL,
    time = "time",
    group_level = NULL,
    fit,
    var = NULL,
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

  if (is.null(var)) {
    if (length(metadata$vars) != 1L) stopf("Supply `var` explicitly for a multivariable fit.")
    variables <- metadata$vars
  } else {
    if (!is.character(var) || !length(var) || anyNA(var) || any(!nzchar(var)) || anyDuplicated(var))
      stopf("`var` must contain unique fitted exposure names.")
    unknown <- setdiff(var, metadata$vars)
    if (length(unknown)) stopf("Unknown fitted exposure(s): ", paste(unknown, collapse=", "), ".")
    variables <- var
  }

  # Validate input and create chronological profile records (oldest to newest).
  records <- list()
  if (using_profile) {
    if (is.numeric(profile) && is.null(dim(profile))) {
      if (length(variables)!=1L) stopf("A numeric `profile` requires exactly one `var`.")
      profile <- stats::setNames(list(as.numeric(profile)), variables)
    }
    if (!is.list(profile)) stopf("`profile` must be a numeric vector or list.")
    if (is.null(names(profile))) {
      if (length(profile)!=length(variables)) stopf("Unnamed `profile` must follow `var` exactly.")
      names(profile) <- variables
    }
    if (!setequal(names(profile), variables)) stopf("`profile` names must match `var` exactly.")
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

  # This is the key correction: every isolated-lag contrast uses the same basis
  # pipeline as the complete trajectory. No crosspred-based parallel route.
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
      contribution = cp,
      percent_contribution = pct,
      stringsAsFactors = FALSE
    )
    if (absolute) {
      base$abs_contribution <- abs(cp)
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
      base$contribution<-cs[,1]; base$contribution_sd<-cs[,2]; base$contribution_lower<-cs[,3]; base$contribution_upper<-cs[,4]
      base$percent_contribution<-ps[,1]; base$percent_contribution_sd<-ps[,2]; base$percent_contribution_lower<-ps[,3]; base$percent_contribution_upper<-ps[,4]; base$percent_contribution_n_defined<-as.integer(ps[,5])
      if (absolute) {
        acs <- t(vapply(
          seq_len(history_length),
          function(j) summarize_finite(abs(CD[, j])),
          numeric(4)
        ))
        base$abs_contribution <- acs[, 1]
        base$abs_contribution_sd <- acs[, 2]
        base$abs_contribution_lower <- acs[, 3]
        base$abs_contribution_upper <- acs[, 4]
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
            contribution = CD[s, ],
            percent_contribution = PD[s, ],
            stringsAsFactors = FALSE
          )
          if (absolute) {
            z$abs_contribution <- abs(z$contribution)
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
