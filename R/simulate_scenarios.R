#' Simulate epidemiological DLNM scenarios
#'
#' Generates predicted outcomes under user-defined epidemiological scenarios
#' using distributed lag nonlinear models (DLNMs).
#'
#' The function combines exposure conditions specified for one or more lag
#' periods, assembles complete exposure histories, and evaluates the resulting
#' outcome predictions through `predict_outcome()`.
#'
#' Scenarios can be supplied either as:
#'
#' - a structured object returned by `simulate_range()`, or
#' - a named list of scenario definitions.
#'
#' Exposure profiles are assembled in chronological order:
#'
#' - the first profile element represents the oldest exposure;
#' - the last profile element represents the most recent exposure.
#'
#' Internal conversion to retrospective lag indexing is handled
#' automatically by the package and is never required from the user.
#'
#' Reference exposure values can be supplied explicitly through
#' `ref_vals`. When omitted, the function attempts to determine
#' appropriate reference values from:
#'
#' 1. the median of each variable in `data`;
#' 2. centering values stored in the fitted DLNM specification.
#'
#' When `uncertainty = FALSE`, predictions are computed using the
#' central coefficient estimates stored in the fitted model.
#'
#' When `uncertainty = TRUE`, uncertainty is propagated through
#' `predict_outcome()` using model-specific coefficient draws:
#'
#' - posterior draws for Bayesian models;
#' - simulated coefficient draws based on the asymptotic covariance
#'   matrix for frequentist models.
#'
#' If `output = "summary"`, simulated predictions are summarized by:
#'
#' - median prediction;
#' - empirical standard deviation;
#' - empirical 95% uncertainty interval.
#'
#' If `output = "samples"`, all simulated predictions are returned
#' without aggregation.
#'
#' For mixed-effects models, `pop_level = TRUE` produces population-level
#' predictions when supported by the fitted engine, whereas
#' `pop_level = FALSE` includes conditional predictions.
#'
#' The function is intended for scenario analysis, intervention
#' assessment, counterfactual simulations, and exposure profile
#' comparisons across multiple lag periods.
#'
#' Conceptually, each scenario represents a complete exposure history.
#' The prediction associated with a scenario therefore reflects the
#' cumulative effect of the entire exposure profile rather than the
#' effect of a single lag or a single exposure value.
#'
#' @export
simulate_scenarios <- function(
    fit,
    scenarios,
    data = NULL,
    periods = NULL,
    ref_vals = NULL,
    pop_level = TRUE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    seed = NULL
) {

  output <- match.arg(output)

  `%||%` <- function(a,b) if(!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if(length(x) <= 1L) return(0)
    stats::sd(x)
  }

  is_whole_number <- function(x){
    is.numeric(x) && length(x)==1L && is.finite(x) &&
      abs(x-round(x)) < sqrt(.Machine$double.eps)
  }

  if (is.null(fit)) stop("`fit` cannot be NULL.")

  if (!is.logical(pop_level) || length(pop_level)!=1L || is.na(pop_level)) {
    stop("`pop_level` must be TRUE or FALSE.")
  }

  if (!is.logical(uncertainty) || length(uncertainty)!=1L || is.na(uncertainty)) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }

  if (!is_whole_number(n_samples) || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  if (uncertainty && n_samples < 2L) {
    stop("`n_samples` must be at least 2 when `uncertainty = TRUE`.")
  }

  if (!is.null(seed)) {
    if (!is_whole_number(seed)) stop("`seed` must be NULL or a finite integer.")
    set.seed(as.integer(seed))
  }

  fit_spec <- attr(fit, "epiexposure_spec")
  fit_vars <- attr(fit, "epiexposure_vars")
  dat_template <- attr(fit, "epiexposure_data_template")

  if (is.null(dat_template) || !is.data.frame(dat_template) || nrow(dat_template) < 1L) {
    stop("The fitted model does not contain `epiexposure_data_template`.")
  }

  if (is.null(fit_spec) || !is.list(fit_spec)) {
    stop("The fitted model does not contain `epiexposure_spec`.")
  }

  scenario_info <- NULL
  periods_from_scenarios <- NULL

  if (is.list(scenarios) && "scenarios" %in% names(scenarios)) {
    scenario_info <- scenarios$info %||% NULL
    periods_from_scenarios <- scenarios$periods %||% NULL
    scenarios <- scenarios$scenarios
  }

  if (!is.null(scenario_info) && "profile_order" %in% names(scenario_info)) {
    if (any(scenario_info$profile_order != "chronological")) {
      stop("All profile-based scenarios must be stored in chronological order.")
    }
  }

  if (!is.list(scenarios) || is.null(names(scenarios)) || any(names(scenarios)=="")) {
    stop("`scenarios` must be a named list.")
  }

  periods <- periods %||% periods_from_scenarios

  if (is.null(periods)) stop("`periods` is missing.")

  req_cols <- c("period","lag_start","lag_end")
  if (!all(req_cols %in% names(periods))) {
    stop("`periods` must contain columns: period, lag_start, lag_end.")
  }

  if (anyDuplicated(periods$period)) stop("`periods$period` must contain unique labels.")

  if (any(!is.finite(periods$lag_start)) || any(!is.finite(periods$lag_end))) {
    stop("`lag_start` and `lag_end` must be finite.")
  }

  if (any(abs(periods$lag_start-round(periods$lag_start)) > 0) ||
      any(abs(periods$lag_end-round(periods$lag_end)) > 0)) {
    stop("`lag_start` and `lag_end` must be integers.")
  }

  if (any(periods$lag_start > periods$lag_end)) {
    stop("Every period must satisfy `lag_start <= lag_end`.")
  }

  vars <- if (!is.null(fit_vars) && length(fit_vars)>0L) fit_vars else unique(unlist(lapply(scenarios,function(s) unique(unlist(lapply(s,names),use.names=FALSE))),use.names=FALSE))

  get_var_lagmax <- function(v){
    if(!is.null(fit_spec[[v]]$max_lag)) return(as.integer(max(fit_spec[[v]]$max_lag)))
    stop("Could not determine max_lag for variable: ",v)
  }

  var_lagmax <- setNames(lapply(vars,get_var_lagmax),vars)

  max_period_lag <- max(periods$lag_end)
  for(v in vars){
    if(max_period_lag > var_lagmax[[v]]){
      stop("Period definitions exceed max_lag for variable '",v,"'.")
    }
  }

  if (is.null(ref_vals)) {
    ref_vals <- setNames(vector("list", length(vars)), vars)
    for(v in vars){
      if(!is.null(data) && is.data.frame(data) && v %in% names(data)){
        ref_vals[[v]] <- as.numeric(stats::median(data[[v]], na.rm=TRUE)); next
      }
      if(!is.null(fit_spec[[v]]$argvar$cen)){
        ref_vals[[v]] <- as.numeric(fit_spec[[v]]$argvar$cen); next
      }
      stop("Could not determine default reference value for variable '",v,"'.")
    }
  }

  lag_to_idx <- function(lags, N){
    idx <- N - lags
    idx[idx >= 1L & idx <= N]
  }

  build_profile <- function(var, scen_values){
    N_var <- var_lagmax[[var]] + 1L
    x <- rep(ref_vals[[var]], N_var)

    for(i in seq_len(nrow(periods))){
      p_id <- periods$period[i]
      if(!is.null(scen_values[[p_id]]) && !is.null(scen_values[[p_id]][[var]])){
        lags <- seq(periods$lag_start[i], periods$lag_end[i])
        idx <- lag_to_idx(lags, N_var)
        val <- scen_values[[p_id]][[var]]
        if(!(length(val) %in% c(1L,length(idx)))){
          stop("Scenario '",p_id,"', variable '",var,"' has invalid length.")
        }
        x[idx] <- val
      }
    }
    x
  }

  expand_scenario_points <- function(scen_values){
    lens <- unlist(lapply(scen_values,function(block) sapply(block,length)), use.names=FALSE)
    lens_gt1 <- unique(lens[lens > 1L])
    if(length(lens_gt1)==0L) return(list(scen_values))
    if(length(lens_gt1)>1L) stop("All varying scenario vectors must have the same length.")
    n_pts <- lens_gt1[1]
    lapply(seq_len(n_pts), function(i){
      lapply(scen_values, function(block){
        lapply(block, function(val) if(length(val)>1L) val[i] else val)
      })
    })
  }

  attach_scenario_info <- function(name, point_index=1L, n_rows=1L, total_points=1L){
    if(is.null(scenario_info)) return(data.frame(scenario=rep(name,n_rows), stringsAsFactors=FALSE))
    df_info <- scenario_info[scenario_info$scenario==name,,drop=FALSE]
    if(nrow(df_info)==0L) return(data.frame(scenario=rep(name,n_rows), stringsAsFactors=FALSE))
    if(nrow(df_info)==1L){ df_info <- df_info[rep(1L,n_rows),,drop=FALSE]; rownames(df_info)<-NULL; return(df_info)}
    if(nrow(df_info)==total_points){ df_info <- df_info[point_index,,drop=FALSE]; df_info <- df_info[rep(1L,n_rows),,drop=FALSE]; rownames(df_info)<-NULL; return(df_info)}
    stop("scenario_info rows mismatch.")
  }

  re_mode <- if(isTRUE(pop_level)) "population" else "conditional"
  pred_output <- if(uncertainty) "samples" else output

  out_list <- purrr::imap(scenarios, function(scen,name){
    pts <- expand_scenario_points(scen)
    pred_list <- vector("list", length(pts))

    for(j in seq_along(pts)){
      profiles_j <- lapply(vars, function(v) build_profile(v, pts[[j]]))
      names(profiles_j) <- vars

      pred_j <- predict_outcome(
        fit = fit,
        profiles = profiles_j,
        re = re_mode,
        allow_new_levels = TRUE,
        type = "response",
        uncertainty = uncertainty,
        output = pred_output,
        n_samples = n_samples,
        seed = seed
      )

      meta_j <- attach_scenario_info(name,j,1L,length(pts))

      if(!uncertainty){
        meta_j$prediction <- pred_j$prediction[1]
        pred_list[[j]] <- meta_j
      } else if(output == "summary"){
        meta_j$prediction <- stats::median(pred_j$prediction, na.rm=TRUE)
        meta_j$sd <- safe_sd(pred_j$prediction)
        meta_j$lower <- safe_quantile(pred_j$prediction)[1]
        meta_j$upper <- safe_quantile(pred_j$prediction)[2]
        pred_list[[j]] <- meta_j
      } else {
        tmp <- meta_j[rep(1L,nrow(pred_j)),,drop=FALSE]
        tmp$sample <- pred_j$sample
        tmp$prediction <- pred_j$prediction
        pred_list[[j]] <- tmp
      }
    }
    do.call(rbind,pred_list)
  })

  out <- do.call(rbind,out_list)
  rownames(out) <- NULL
  out
}
