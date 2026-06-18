#' Find the best DLNM model structure using LOOCV and Lin's CCC
#'
#' Tests multiple combinations of exposure variables, exposure-basis degrees of
#' freedom, and lag-basis degrees of freedom using leave-one-out cross-validation
#' (LOOCV). Candidate models are ranked using Lin's concordance correlation
#' coefficient (CCC), alongside the bias correction factor (`Cb`), Pearson
#' correlation (`rho`), RMSE, and MAE.
#'
#' This function is designed to work directly with the same long-format dataset
#' used by the EpiExposure workflow:
#'
#' \enumerate{
#'   \item `define_exposure()`
#'   \item `build_design()`
#'   \item `prepare_response()`
#'   \item `fit_epidlnm()`
#' }
#'
#' The input data should contain one row per time point within each epidemic,
#' with exposure variables measured over time and the response repeated for each
#' row within the same epidemic.
#'
#' @param epi_data Long-format data frame containing the response, time index,
#'   epidemic/group identifier, and daily exposure variables.
#'
#' @param response Character. Name of the response variable in `epi_data`.
#'   The response may be repeated across rows within each epidemic/group.
#'
#' @param group Character. Name of the grouping variable used for LOOCV.
#'   Default is `"epi_id"`.
#'
#' @param time Character. Name of the temporal ordering variable.
#'   Default is `"dpp"`.
#'
#' @param vars Character vector of candidate exposure variables to test.
#'
#' @param lag_max Integer. Maximum lag used to build the DLNM crossbasis.
#'
#' @param df_var_grid Numeric vector. Candidate degrees of freedom for the
#'   exposure-response basis.
#'
#' @param df_lag_grid Numeric vector. Candidate degrees of freedom for the
#'   lag-response basis.
#'
#' @param min_vars Integer. Minimum number of exposure variables per candidate
#'   model. Default is `1`.
#'
#' @param max_vars Integer or `NULL`. Maximum number of exposure variables per
#'   candidate model. If `NULL`, all sizes up to `length(vars)` are tested.
#'
#' @param var_sets Optional list of character vectors. If supplied, these exact
#'   variable combinations are tested and `min_vars`/`max_vars` are ignored.
#'
#' @param fun_var Character. Exposure-response basis function passed to
#'   `define_exposure()`. Default is `"ns"`.
#'
#' @param fun_lag Character. Lag-response basis function passed to
#'   `define_exposure()`. Default is `"ns"`.
#'
#' @param model_engine Character. Model engine passed to `fit_epidlnm()`.
#'   Supported values are those supported by `fit_epidlnm()`: `"glm"`,
#'   `"glmmTMB"`, `"gam"`, `"gamm"`, `"gls"`, `"spamm"`, `"brms"`, `"inla"`,
#'   and `"bdlnm"`.
#'
#' @param family Character or family object. Model family passed to
#'   `prepare_response()` and `fit_epidlnm()`. For example: `"beta"`,
#'   `"gaussian"`, `"poisson"`, `"gamma"`, `"binomial"`, `"negbin"`,
#'   `"negative_binomial"`, or a supported family object when appropriate.
#'
#' @param random_effect Optional character. Name of the random-intercept
#'   grouping variable used by `fit_epidlnm()`; for example `"siteyear"`.
#'
#' @param min_success Integer. Minimum number of successful LOOCV folds required
#'   for a candidate model to be ranked. Default is `2`.
#'
#' @param top_n Integer or `Inf`. Number of top-ranked models to return.
#'   Default is `Inf`, returning all successful candidates.
#'
#' @param verbose Logical. If `TRUE`, prints progress messages.
#'
#' @param ... Additional arguments passed to `fit_epidlnm()`.
#'
#' @return A data.frame ranked by decreasing CCC. Columns include:
#' \itemize{
#'   \item `rank`
#'   \item `model_id`
#'   \item `df_var`
#'   \item `df_lag`
#'   \item `vars`
#'   \item `n_vars`
#'   \item `CCC`
#'   \item `Cb`
#'   \item `rho`
#'   \item `RMSE`
#'   \item `MAE`
#'   \item `n_folds`
#'   \item `n_success`
#'   \item `n_failed`
#' }
#'
#' The returned object includes two attributes:
#' \itemize{
#'   \item `"predictions"`: fold-level observed and predicted values.
#'   \item `"failures"`: candidate/fold-level failure information, if any.
#' }
#'
#' @details
#' **Input structure.** This function expects a single long-format dataset,
#' similar to the object passed through `define_exposure()`, `build_design()`,
#' `prepare_response()`, and `fit_epidlnm()`. The response can be repeated
#' across daily rows within each epidemic/group, but it must be unique within
#' each group after removing duplicates.
#'
#' **Temporal dimension.** The data are ordered internally by `time` within
#' each `group` during basis construction. Each group must contain at least
#' `lag_max + 1` observations for each tested exposure variable.
#'
#' **Cross-validation.** LOOCV is performed by leaving out one `group` level at
#' a time. DLNM basis templates are built on the training data only and then
#' applied to the held-out group.
#'
#' **Prediction across engines.** The function first attempts engine-specific
#' prediction on the response scale. If direct prediction is unavailable for an
#' engine, it falls back to a fixed-effect linear predictor using the model
#' coefficients and applies an inverse-link transformation based on `family`.
#'
#' **Checklist notes.**
#' \itemize{
#'   \item Temporal dimension is controlled through `lag_max` and checked per group.
#'   \item The DLNM API is delegated to `define_exposure()` and `build_design()`.
#'   \item Coefficient extraction is used only as a prediction fallback.
#'   \item Coefficient ordering follows the stable `cb_<var>_<index>` convention
#'     created by `build_design()`.
#'   \item `b_cb_` standardization is only relevant for posterior-draw functions;
#'     this LOOCV ranking function does not summarize posterior surfaces.
#'   \item Surface summarization is not performed; models are ranked by predictive
#'     agreement via Lin's CCC.
#' }
#'
#' @export
find_bestfit <- function(
    epi_data,
    response = "y",
    group = "epi_id",
    time = "dpp",
    vars,
    lag_max,
    df_var_grid = c(3, 4, 5),
    df_lag_grid = c(3, 4, 5),
    min_vars = 1,
    max_vars = NULL,
    var_sets = NULL,
    fun_var = "ns",
    fun_lag = "ns",
    model_engine = "glmmTMB",
    family = "beta",
    random_effect = NULL,
    min_success = 2,
    top_n = Inf,
    verbose = TRUE,
    ...
) {
  
  # ------------------------------------------------------------
  # Basic validations
  # ------------------------------------------------------------
  if (!is.data.frame(epi_data)) {
    stop("`epi_data` must be a data.frame.")
  }
  
  if (!is.character(response) || length(response) != 1L) {
    stop("`response` must be a single character string.")
  }
  
  if (!response %in% names(epi_data)) {
    stop("`response` ('", response, "') was not found in `epi_data`.")
  }
  
  if (!is.character(group) || length(group) != 1L) {
    stop("`group` must be a single character string.")
  }
  
  if (!group %in% names(epi_data)) {
    stop("`group` ('", group, "') was not found in `epi_data`.")
  }
  
  if (!is.character(time) || length(time) != 1L) {
    stop("`time` must be a single character string.")
  }
  
  if (!time %in% names(epi_data)) {
    stop("`time` ('", time, "') was not found in `epi_data`.")
  }
  
  if (!is.character(vars) || length(vars) < 1L) {
    stop("`vars` must be a non-empty character vector.")
  }
  
  missing_vars <- setdiff(vars, names(epi_data))
  if (length(missing_vars) > 0) {
    stop(
      "The following variables in `vars` were not found in `epi_data`: ",
      paste(missing_vars, collapse = ", "),
      "."
    )
  }
  
  if (!is.numeric(lag_max) || length(lag_max) < 1L ||
      any(!is.finite(lag_max)) || max(lag_max) < 0) {
    stop("`lag_max` must be a non-negative integer or numeric vector.")
  }
  
  lag_max <- as.integer(max(lag_max))
  
  if (!is.numeric(df_var_grid) || length(df_var_grid) < 1L) {
    stop("`df_var_grid` must be a non-empty numeric vector.")
  }
  
  if (!is.numeric(df_lag_grid) || length(df_lag_grid) < 1L) {
    stop("`df_lag_grid` must be a non-empty numeric vector.")
  }
  
  if (!is.numeric(min_vars) || length(min_vars) != 1L || min_vars < 1) {
    stop("`min_vars` must be a positive integer.")
  }
  
  min_vars <- as.integer(min_vars)
  
  if (is.null(max_vars)) {
    max_vars <- length(vars)
  }
  
  if (!is.numeric(max_vars) || length(max_vars) != 1L ||
      max_vars < min_vars || max_vars > length(vars)) {
    stop("`max_vars` must be between `min_vars` and length(`vars`).")
  }
  
  max_vars <- as.integer(max_vars)
  
  if (!is.null(random_effect)) {
    if (!is.character(random_effect) || length(random_effect) != 1L) {
      stop("`random_effect` must be NULL or a single character string.")
    }
    
    if (!random_effect %in% names(epi_data)) {
      stop("`random_effect` ('", random_effect, "') was not found in `epi_data`.")
    }
  }
  
  if (!is.numeric(min_success) || length(min_success) != 1L ||
      !is.finite(min_success) || min_success < 2) {
    stop("`min_success` must be an integer >= 2.")
  }
  
  min_success <- as.integer(min_success)
  
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("`verbose` must be TRUE or FALSE.")
  }
  
  model_engine <- match.arg(
    model_engine,
    choices = c("glm", "glmmTMB", "gam", "gamm", "gls", "spamm", "brms", "inla", "bdlnm")
  )
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  # ------------------------------------------------------------
  # Internal standardization to EpiExposure pipeline names
  # ------------------------------------------------------------
  dat_long <- epi_data
  
  dat_long$epi_id <- dat_long[[group]]
  dat_long$dpp <- dat_long[[time]]
  dat_long$y <- dat_long[[response]]
  
  if (!is.null(random_effect)) {
    dat_long[[random_effect]] <- epi_data[[random_effect]]
  }
  
  # Validate a single response per group
  y_check <- unique(dat_long[, c("epi_id", "y"), drop = FALSE])
  
  if (any(duplicated(y_check$epi_id))) {
    stop(
      "Multiple distinct response values were found within at least one `group`. ",
      "The response must be constant/repeated within each group."
    )
  }
  
  # ------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------
  ccc_lins <- function(obs, pred) {
    
    ok <- is.finite(obs) & is.finite(pred)
    obs <- obs[ok]
    pred <- pred[ok]
    
    if (length(obs) < 2L) {
      return(list(
        CCC = NA_real_,
        Cb = NA_real_,
        rho = NA_real_,
        RMSE = NA_real_,
        MAE = NA_real_
      ))
    }
    
    mx <- mean(obs)
    my <- mean(pred)
    
    vx <- stats::var(obs)
    vy <- stats::var(pred)
    sxy <- stats::cov(obs, pred)
    
    rho <- suppressWarnings(stats::cor(obs, pred))
    
    CCC <- (2 * sxy) / (vx + vy + (mx - my)^2)
    
    Cb <- if (is.finite(rho) && abs(rho) > .Machine$double.eps) {
      CCC / rho
    } else {
      NA_real_
    }
    
    RMSE <- sqrt(mean((obs - pred)^2))
    MAE <- mean(abs(obs - pred))
    
    list(
      CCC = as.numeric(CCC),
      Cb = as.numeric(Cb),
      rho = as.numeric(rho),
      RMSE = as.numeric(RMSE),
      MAE = as.numeric(MAE)
    )
  }
  
  get_linkinv_from_family <- function(family) {
    
    if (inherits(family, "family") && !is.null(family$linkinv)) {
      return(family$linkinv)
    }
    
    if (is.character(family)) {
      if (family %in% c("beta", "binomial")) return(stats::plogis)
      if (family %in% c("poisson", "gamma", "negbin", "negative_binomial")) return(exp)
      if (family %in% c("gaussian")) return(identity)
    }
    
    identity
  }
  
  extract_fixed_coef <- function(model) {
    
    if (inherits(model, "glmmTMB")) {
      return(glmmTMB::fixef(model)$cond)
    }
    
    if (inherits(model, "merMod")) {
      return(lme4::fixef(model))
    }
    
    if (inherits(model, "gamm")) {
      return(stats::coef(model$gam))
    }
    
    if (inherits(model, "lme")) {
      return(nlme::fixef(model))
    }
    
    if (inherits(model, "gls")) {
      return(stats::coef(model))
    }
    
    if (inherits(model, "gam")) {
      return(stats::coef(model))
    }
    
    if (inherits(model, "HLfit")) {
      return(spaMM::fixef(model))
    }
    
    if (inherits(model, "brmsfit")) {
      fe <- brms::fixef(model)
      beta <- fe[, "Estimate"]
      names(beta) <- rownames(fe)
      return(beta)
    }
    
    if (inherits(model, "inla")) {
      return(model$summary.fixed$mean)
    }
    
    if (inherits(model, "bdlnm")) {
      if (!is.null(model$coefficients.summary)) {
        return(model$coefficients.summary[, "mean"])
      }
    }
    
    beta <- tryCatch(stats::coef(model), error = function(e) NULL)
    
    if (is.null(beta) || !is.numeric(beta)) {
      stop("Could not extract fixed-effect coefficients for prediction fallback.")
    }
    
    beta
  }
  
  predict_response_engine <- function(fit, newdata, family, model_engine) {
    
    # First try engine-specific prediction
    pred_try <- tryCatch({
      
      if (inherits(fit, "glmmTMB")) {
        stats::predict(fit, newdata = newdata, type = "response", allow.new.levels = TRUE)
        
      } else if (inherits(fit, "brmsfit")) {
        if (!requireNamespace("brms", quietly = TRUE)) {
          stop("Package 'brms' is required for brms prediction.")
        }
        ep <- brms::posterior_epred(
          fit,
          newdata = newdata,
          re_formula = NA,
          allow_new_levels = TRUE
        )
        colMeans(ep)
        
      } else if (inherits(fit, "gamm")) {
        stats::predict(fit$gam, newdata = newdata, type = "response")
        
      } else if (inherits(fit, "gam")) {
        stats::predict(fit, newdata = newdata, type = "response")
        
      } else if (inherits(fit, "glm")) {
        stats::predict(fit, newdata = newdata, type = "response")
        
      } else if (inherits(fit, "gls")) {
        stats::predict(fit, newdata = newdata)
        
      } else if (inherits(fit, "lme")) {
        stats::predict(fit, newdata = newdata, level = 0)
        
      } else if (inherits(fit, "HLfit")) {
        stats::predict(fit, newdata = newdata, type = "response")
        
      } else {
        stop("No direct prediction method used.")
      }
      
    }, error = function(e) NULL)
    
    if (!is.null(pred_try) && all(is.finite(as.numeric(pred_try)))) {
      return(as.numeric(pred_try))
    }
    
    # Fallback: fixed-effect linear predictor
    beta <- extract_fixed_coef(fit)
    
    nm <- names(beta)
    
    if (is.null(nm)) {
      stop("Prediction fallback failed: coefficients have no names.")
    }
    
    eta <- rep(0, nrow(newdata))
    
    if ("(Intercept)" %in% nm) {
      eta <- eta + as.numeric(beta["(Intercept)"])
    }
    
    fixed_terms <- intersect(setdiff(nm, "(Intercept)"), names(newdata))
    
    if (length(fixed_terms) == 0) {
      stop("Prediction fallback failed: no matching fixed-effect columns found in `newdata`.")
    }
    
    X <- as.matrix(newdata[, fixed_terms, drop = FALSE])
    b <- as.numeric(beta[fixed_terms])
    
    eta <- eta + as.numeric(X %*% b)
    
    linkinv <- get_linkinv_from_family(family)
    
    as.numeric(linkinv(eta))
  }
  
  check_temporal_coverage <- function(dat, vars_use, lag_max) {
    
    bad <- dat |>
      dplyr::group_by(epi_id) |>
      dplyr::summarise(n_days = dplyr::n_distinct(dpp), .groups = "drop") |>
      dplyr::filter(n_days < lag_max + 1L)
    
    if (nrow(bad) > 0) {
      stop(
        "Some groups do not have enough temporal coverage for lag_max.\n",
        "Required days per group: ", lag_max + 1L, "\n",
        "Example problematic group(s): ",
        paste(utils::head(bad$epi_id, 5), collapse = ", ")
      )
    }
    
    for (v in vars_use) {
      bad_v <- dat |>
        dplyr::group_by(epi_id) |>
        dplyr::summarise(n_finite = sum(is.finite(.data[[v]])), .groups = "drop") |>
        dplyr::filter(n_finite < lag_max + 1L)
      
      if (nrow(bad_v) > 0) {
        stop(
          "Some groups do not have enough finite observations for variable '", v, "'.\n",
          "Required finite observations per group: ", lag_max + 1L
        )
      }
    }
    
    invisible(TRUE)
  }
  
  # ------------------------------------------------------------
  # Candidate variable sets
  # ------------------------------------------------------------
  if (is.null(var_sets)) {
    
    sizes <- seq.int(min_vars, max_vars)
    
    var_sets <- unlist(
      lapply(sizes, function(k) {
        utils::combn(vars, k, simplify = FALSE)
      }),
      recursive = FALSE
    )
    
  } else {
    
    if (!is.list(var_sets) || length(var_sets) == 0) {
      stop("`var_sets` must be NULL or a non-empty list of character vectors.")
    }
    
    bad_sets <- unique(unlist(lapply(var_sets, function(s) setdiff(s, vars))))
    
    if (length(bad_sets) > 0) {
      stop(
        "Some variables in `var_sets` are not present in `vars`: ",
        paste(bad_sets, collapse = ", "),
        "."
      )
    }
  }
  
  # ------------------------------------------------------------
  # Common checks
  # ------------------------------------------------------------
  check_temporal_coverage(dat_long, vars, lag_max)
  
  fold_ids <- unique(dat_long$epi_id)
  
  if (length(fold_ids) < 2L) {
    stop("At least two groups are required for LOOCV.")
  }
  
  # ------------------------------------------------------------
  # Main loop
  # ------------------------------------------------------------
  results <- list()
  predictions <- list()
  failures <- list()
  
  model_id <- 1L
  pred_id <- 1L
  fail_id <- 1L
  
  total_candidates <- length(df_var_grid) * length(df_lag_grid) * length(var_sets)
  candidate_id <- 1L
  
  for (df_var in df_var_grid) {
    for (df_lag in df_lag_grid) {
      for (vars_use in var_sets) {
        
        if (verbose) {
          message(
            "[", candidate_id, "/", total_candidates, "] ",
            "df_var = ", df_var,
            ", df_lag = ", df_lag,
            ", vars = ", paste(vars_use, collapse = " + ")
          )
        }
        
        candidate_id <- candidate_id + 1L
        
        obs_all <- rep(NA_real_, length(fold_ids))
        pred_all <- rep(NA_real_, length(fold_ids))
        success <- rep(FALSE, length(fold_ids))
        
        for (fold_i in seq_along(fold_ids)) {
          
          test_group <- fold_ids[fold_i]
          train_groups <- setdiff(fold_ids, test_group)
          
          dat_train_long <- dat_long[dat_long$epi_id %in% train_groups, , drop = FALSE]
          dat_test_long  <- dat_long[dat_long$epi_id == test_group, , drop = FALSE]
          
          fold_res <- tryCatch({
            
            cb_templates <- define_exposure(
              wx_long = dat_train_long,
              vars = vars_use,
              lag_max = lag_max,
              df_var = df_var,
              df_lag = df_lag,
              fun_var = fun_var,
              fun_lag = fun_lag
            )
            
            X_train <- build_design(
              wx_long = dat_train_long,
              cb_templates = cb_templates,
              lag_max = lag_max,
              include_response = TRUE
            )
            
            X_test <- build_design(
              wx_long = dat_test_long,
              cb_templates = cb_templates,
              lag_max = lag_max,
              include_response = TRUE
            )
            
            # Add metadata columns required by the model, such as random effects
            if (!is.null(random_effect)) {
              
              meta_train <- unique(dat_train_long[, c("epi_id", random_effect), drop = FALSE])
              meta_test  <- unique(dat_test_long[,  c("epi_id", random_effect), drop = FALSE])
              
              X_train <- merge(X_train, meta_train, by = "epi_id", all.x = TRUE)
              X_test  <- merge(X_test,  meta_test,  by = "epi_id", all.x = TRUE)
            }
            
            dat_train <- prepare_response(
              dat = X_train,
              y_var = "y",
              family_choice = if (is.character(family)) family else "gaussian"
            )
            
            dat_test <- prepare_response(
              dat = X_test,
              y_var = "y",
              family_choice = if (is.character(family)) family else "gaussian"
            )
            
            fit <- fit_epidlnm(
              dat = dat_train,
              model_engine = model_engine,
              family = family,
              random_effect = random_effect,
              epiexposure_spec = attr(cb_templates, "spec"),
              basis_objects = cb_templates,
              ...
            )
            
            pred <- predict_response_engine(
              fit = fit,
              newdata = dat_test,
              family = family,
              model_engine = model_engine
            )
            
            list(
              obs = dat_test$y_model[1],
              pred = as.numeric(pred[1])
            )
            
          }, error = function(e) {
            
            failures[[fail_id]] <<- data.frame(
              model_id = model_id,
              fold = fold_i,
              group = as.character(test_group),
              df_var = df_var,
              df_lag = df_lag,
              vars = paste(vars_use, collapse = " + "),
              error = conditionMessage(e),
              stringsAsFactors = FALSE
            )
            
            fail_id <<- fail_id + 1L
            NULL
          })
          
          if (!is.null(fold_res)) {
            
            obs_all[fold_i] <- fold_res$obs
            pred_all[fold_i] <- fold_res$pred
            success[fold_i] <- TRUE
            
            predictions[[pred_id]] <- data.frame(
              model_id = model_id,
              fold = fold_i,
              group = as.character(test_group),
              df_var = df_var,
              df_lag = df_lag,
              vars = paste(vars_use, collapse = " + "),
              observed = fold_res$obs,
              predicted = fold_res$pred,
              stringsAsFactors = FALSE
            )
            
            pred_id <- pred_id + 1L
          }
        }
        
        n_success <- sum(success)
        n_failed <- length(fold_ids) - n_success
        
        if (n_success >= min_success) {
          
          metrics <- ccc_lins(
            obs = obs_all[success],
            pred = pred_all[success]
          )
          
          results[[length(results) + 1L]] <- data.frame(
            model_id = model_id,
            df_var = df_var,
            df_lag = df_lag,
            vars = paste(vars_use, collapse = " + "),
            n_vars = length(vars_use),
            CCC = metrics$CCC,
            Cb = metrics$Cb,
            rho = metrics$rho,
            RMSE = metrics$RMSE,
            MAE = metrics$MAE,
            n_folds = length(fold_ids),
            n_success = n_success,
            n_failed = n_failed,
            stringsAsFactors = FALSE
          )
        }
        
        model_id <- model_id + 1L
      }
    }
  }
  
  if (length(results) == 0) {
    stop("No candidate model produced enough successful LOOCV predictions.")
  }
  
  results_df <- do.call(rbind, results)
  
  results_df <- results_df[order(
    -results_df$CCC,
    -results_df$Cb,
    -results_df$rho,
    results_df$RMSE,
    results_df$MAE
  ), , drop = FALSE]
  
  results_df$rank <- seq_len(nrow(results_df))
  results_df <- results_df[, c("rank", setdiff(names(results_df), "rank")), drop = FALSE]
  
  if (is.finite(top_n)) {
    results_df <- utils::head(results_df, top_n)
  }
  
  pred_df <- if (length(predictions) > 0) {
    do.call(rbind, predictions)
  } else {
    data.frame()
  }
  
  fail_df <- if (length(failures) > 0) {
    do.call(rbind, failures)
  } else {
    data.frame()
  }
  
  attr(results_df, "predictions") <- pred_df
  attr(results_df, "failures") <- fail_df
  
  results_df
}