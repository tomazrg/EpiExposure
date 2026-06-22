#' Create prediction and lag-contribution ensembles from best-fit DLNM candidate models
#'
#' Builds ensemble predictions from the out-of-fold predictions produced by
#' `find_bestfit()` and, optionally, builds an ensemble of lag-specific
#' contributions using lag-decomposition outputs computed for the selected
#' candidate models.
#'
#' The function supports three ensemble scopes:
#'
#' \itemize{
#'   \item `"model"`: ensemble of model predictions only.
#'   \item `"lag"`: ensemble of lag-specific contributions only.
#'   \item `"both"`: both prediction ensemble and lag-contribution ensemble.
#' }
#'
#' For prediction ensembles, the function uses the out-of-fold predictions
#' stored in `attr(bestfit, "predictions")` or supplied through `predictions`.
#' These predictions are expected to contain the standard columns:
#' `model_id`, `group`, `fold`, `observed`, and `predicted`.
#'
#' For lag ensembles, the function can use either:
#'
#' \itemize{
#'   \item a user-supplied `lag_data` data frame; or
#'   \item fitted models stored in `attr(bestfit, "fits")`, created by
#'   `find_bestfit(..., keep_fits = TRUE)`, together with `epi_data`.
#' }
#'
#' When `lag_data = NULL` and `ensemble_scope` is `"lag"` or `"both"`, the
#' function automatically computes lag-specific decompositions by calling
#' `compute_ecilag()` for each selected fitted model.
#'
#' @param bestfit A data.frame returned by `find_bestfit()`. It must contain
#'   model-level performance metrics and should have attributes named
#'   `"predictions"` and, for automatic lag ensembles, `"fits"`.
#'
#' @param predictions Optional data.frame of fold-level predictions. If `NULL`,
#'   the function uses `attr(bestfit, "predictions")`. The expected standard
#'   columns are `model_id`, `group`, `fold`, `observed`, and `predicted`.
#'
#' @param lag_data Optional data.frame containing lag-specific contributions by
#'   model. If supplied, it is used directly for lag ensembles. If `NULL` and
#'   `ensemble_scope` is `"lag"` or `"both"`, the function attempts to compute
#'   it automatically using fitted models stored in `attr(bestfit, "fits")`.
#'   The expected standard columns are `model_id`, `var`, `lag`,
#'   `contribution`, and optionally `weight`.
#'
#' @param fit_list Optional named list of fitted models. If `NULL`, the function
#'   uses `attr(bestfit, "fits")`. Names must correspond to model IDs.
#'
#' @param epi_data Optional long-format data frame used to compute lag-specific
#'   decompositions automatically through `compute_ecilag()`. Required when
#'   `lag_data = NULL` and `ensemble_scope` is `"lag"` or `"both"`.
#'
#' @param group Character. Column in `epi_data` identifying groups, such as
#'   `"epi_id"`. Used only when automatic lag decomposition is requested.
#'
#' @param var Optional character vector indicating which exposure variables
#'   should be used for lag decomposition. If `NULL`, variables are resolved
#'   from each fitted model according to `lg_strategy`.
#'
#' @param reverse Logical. Passed to `compute_ecilag()` when automatic lag
#'   decomposition is requested. Use `TRUE` when the profiles in `epi_data` are
#'   ordered from oldest to most recent and need to be reversed so that the first
#'   value corresponds to lag 0.
#'
#' @param compute_ecilag_args Optional named list of additional arguments passed
#'   to `compute_ecilag()`, such as `uncertainty`, `output`, `n_samples`, `eps`,
#'   `center`, or `absolute`.
#'
#' @param ensemble_scope Character. Which ensemble to compute. Options are
#'   `"model"`, `"lag"`, and `"both"`. Default is `"model"`.
#'
#' @param method Character. Ensemble strategy. Options are `"best"`,
#'   `"hard_voting"`, `"unweighted"`, `"weighted"`, and `"stacked"`.
#'
#'   The same method is used for model prediction ensembles and lag-contribution
#'   ensembles. However, `"stacked"` is only supported when
#'   `ensemble_scope = "model"`, because lag-specific contributions do not have
#'   an observed response target for stacking.
#'
#' @param top_n Integer. Number of top-ranked models to include in the ensemble.
#'   Ignored when `model_ids` is supplied. Default is `3`.
#'
#' @param model_ids Optional vector of model IDs to use in the ensemble. If
#'   supplied, `top_n` is ignored.
#'
#' @param weight_metric Character. Column in `bestfit` used to compute weights
#'   for `"weighted"` ensembles. Default is `"CCC"`.
#'
#' @param weight_transform Character. Weighting rule used when
#'   `method = "weighted"`. Options are:
#'   \describe{
#'     \item{"softmax"}{Uses a softmax transformation of `weight_metric`.}
#'     \item{"positive"}{Uses positive metric values normalized to sum to one.}
#'     \item{"rank_inverse"}{Uses inverse rank weights.}
#'     \item{"uniform"}{Uses equal weights.}
#'   }
#'
#' @param stacking_model Character. Meta-model used when `method = "stacked"`.
#'   Options are `"ridge"` and `"lm"`.
#'
#' @param stack_objective Character. Objective used to translate stacking into
#'   model weights. Options are:
#'   \describe{
#'     \item{"regularized"}{Default. Uses ridge-stacking weights. This is stable
#'     when candidate models are correlated.}
#'     \item{"min_error"}{Uses least-squares stacking weights that minimize
#'     prediction error on out-of-fold predictions.}
#'     \item{"constrained"}{Projects stacking weights to the non-negative simplex,
#'     forcing weights to be non-negative and sum to one.}
#'   }
#'
#' @param lambda Numeric. Ridge penalty used when `stacking_model = "ridge"` or
#'   when `stack_objective = "regularized"`. Default is `1e-6`.
#'
#' @param stack_intercept Logical. If `TRUE`, includes an intercept in the
#'   stacking model. Default is `TRUE`.
#'
#' @param model_col Character. Name of the model ID column in `bestfit`,
#'   `predictions`, and `lag_data`. Default is `"model_id"`.
#'
#' @param id_cols Character vector identifying each out-of-fold prediction row.
#'   Default is `c("group", "fold")`, matching the standard output of
#'   `find_bestfit()`.
#'
#' @param lag_group_cols Character vector of columns used to align lag
#'   contributions before averaging. Default is `c("var", "lag")`.
#'   If `lag_data` contains a grouping column such as `"epi_id"` or `"group"`,
#'   include it here to compute lag ensembles separately by group.
#'
#' @param lg_strategy Character. Strategy used to decide which variables are
#'   decomposed for each selected model when `lag_data` is generated
#'   automatically. Options are:
#'   \describe{
#'     \item{"requested_available"}{Default. Uses the intersection between
#'     `var` and the variables actually present in each fitted model. Models
#'     with none of the requested variables are skipped for lag decomposition.}
#'     \item{"model"}{Ignores `var` and uses all exposure variables stored in
#'     each fitted model.}
#'     \item{"strict"}{Requires all variables in `var` to be present in every
#'     selected fitted model; otherwise an error is raised.}
#'   }
#'
#' @param rl_weights Logical. If `TRUE`, model weights are renormalized within
#'   each lag-ensemble grouping unit defined by `lag_group_cols`. This is useful
#'   when not all selected models contain the same variables, because the weights
#'   are rescaled over the models that actually contribute to each variable/lag
#'   combination. Default is `TRUE`.
#'
#' @param hb Logical. Higher is better. If `TRUE`, larger values of
#'   `weight_metric` indicate better models. Default is `TRUE`.
#'
#' @param verbose Logical. If `TRUE`, prints informative messages.
#'
#' @return A list containing:
#' \itemize{
#'   \item `ensemble_summary`: performance metrics for model-level ensemble
#'   predictions, when requested.
#'   \item `ensemble_predictions`: observed and ensemble-predicted values,
#'   when model-level ensemble is requested.
#'   \item `ensemble_by_lag`: lag-level ensemble contributions, when requested.
#'   \item `lag_data`: lag-specific contribution table used to build the lag
#'   ensemble.
#'   \item `model_weights`: model-level and lag-level weights assigned to
#'   selected models.
#'   \item `selected_models`: subset of `bestfit` used in the ensemble.
#'   \item `method`: ensemble method used.
#'   \item `ensemble_scope`: requested ensemble scope.
#'   \item `lg_strategy`: variable-resolution strategy used for lag ensemble.
#'   \item `rl_weights`: whether lag weights were renormalized.
#' }
#'
#' @details
#' `compute_ecilag()` estimates lag-specific contributions within a single
#' model. Those contributions are already expressed on the lag scale and are
#' therefore comparable across models. In contrast, raw DLNM coefficients are
#' not directly comparable across models because different models may use
#' different basis dimensions or smoothing structures.
#'
#' The lag ensemble combines lag-specific contributions, not raw coefficients:
#'
#' \deqn{
#' C_{\mathrm{ens},l} = \sum_m \alpha_m C_{m,l}
#' }
#'
#' where \eqn{C_{m,l}} is the lag-specific contribution from model \eqn{m} at
#' lag \eqn{l}, and \eqn{\alpha_m} is the model weight.
#'
#' If `rl_weights = TRUE`, the weights are renormalized within each grouping
#' unit used for the lag ensemble:
#'
#' \deqn{
#' \tilde{\alpha}_m =
#' \frac{\alpha_m}{\sum_{j \in \mathcal{M}_{v,l}} \alpha_j}
#' }
#'
#' where \eqn{\mathcal{M}_{v,l}} is the set of selected models contributing to
#' a given variable/lag combination.
#'
#' The relative lag contribution is computed as:
#'
#' \deqn{
#' RI_{\mathrm{ens},l} =
#' 100 \times
#' \frac{|C_{\mathrm{ens},l}|}
#' {\sum_l |C_{\mathrm{ens},l}|}
#' }
#'
#' This quantity sums to approximately 100 percent across lags within each
#' variable/group combination.
#'
#' For continuous responses, classical hard voting is not directly applicable.
#' Therefore, `"hard_voting"` is implemented as hard model selection, equivalent
#' to using the single best-ranked model.
#'
#' Stacking is only supported for model-level prediction ensembles. It is not
#' supported for lag ensembles because lag-specific contributions do not have an
#' observed response target that can be used to train a stacking model.
#'
#' @export
ensemble_bestfit <- function(
    bestfit,
    predictions = NULL,
    lag_data = NULL,
    fit_list = NULL,
    epi_data = NULL,
    group = "epi_id",
    var = NULL,
    reverse = FALSE,
    compute_ecilag_args = list(),
    ensemble_scope = c("model", "lag", "both"),
    method = c("unweighted", "weighted", "hard_voting", "best", "stacked"),
    top_n = 3,
    model_ids = NULL,
    weight_metric = "CCC",
    weight_transform = c("softmax", "positive", "rank_inverse", "uniform"),
    stacking_model = c("ridge", "lm"),
    stack_objective = c("regularized", "min_error", "constrained"),
    lambda = 1e-6,
    stack_intercept = TRUE,
    model_col = "model_id",
    id_cols = c("group", "fold"),
    lag_group_cols = c("var", "lag"),
    lg_strategy = c("requested_available", "model", "strict"),
    rl_weights = TRUE,
    hb = TRUE,
    verbose = TRUE
) {
  
  # ------------------------------------------------------------
  # Internal standard column names
  # ------------------------------------------------------------
  response_col <- "observed"
  prediction_col <- "predicted"
  lag_col <- "lag"
  lag_contribution_col <- "contribution"
  lag_weight_col <- "weight"
  
  ensemble_scope <- match.arg(ensemble_scope)
  method <- match.arg(method)
  weight_transform <- match.arg(weight_transform)
  stacking_model <- match.arg(stacking_model)
  stack_objective <- match.arg(stack_objective)
  lg_strategy <- match.arg(lg_strategy)
  
  # ------------------------------------------------------------
  # Block conceptually invalid combinations
  # ------------------------------------------------------------
  if (ensemble_scope %in% c("lag", "both") && method == "stacked") {
    stop(
      "`method = 'stacked'` is only supported when `ensemble_scope = 'model'`. ",
      "Stacking requires observed responses for training the meta-model, and ",
      "lag-specific contributions do not have an observed response target. ",
      "Use `method = 'weighted'`, `method = 'unweighted'`, or `method = 'best'` ",
      "for lag-based ensembles."
    )
  }
  
  # ------------------------------------------------------------
  # Validations
  # ------------------------------------------------------------
  if (!is.data.frame(bestfit)) {
    stop("`bestfit` must be a data.frame returned by `find_bestfit()`.")
  }
  
  if (!model_col %in% names(bestfit)) {
    stop("`model_col` ('", model_col, "') was not found in `bestfit`.")
  }
  
  if (!weight_metric %in% names(bestfit)) {
    stop("`weight_metric` ('", weight_metric, "') was not found in `bestfit`.")
  }
  
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
    stop("`top_n` must be a positive integer.")
  }
  
  top_n <- as.integer(top_n)
  
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be a non-negative numeric scalar.")
  }
  
  if (!is.logical(stack_intercept) || length(stack_intercept) != 1L) {
    stop("`stack_intercept` must be TRUE or FALSE.")
  }
  
  if (!is.logical(reverse) || length(reverse) != 1L) {
    stop("`reverse` must be TRUE or FALSE.")
  }
  
  if (!is.list(compute_ecilag_args)) {
    stop("`compute_ecilag_args` must be a named list.")
  }
  
  if (!is.logical(rl_weights) || length(rl_weights) != 1L) {
    stop("`rl_weights` must be TRUE or FALSE.")
  }
  
  if (!is.logical(hb) || length(hb) != 1L) {
    stop("`hb` must be TRUE or FALSE.")
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("`verbose` must be TRUE or FALSE.")
  }
  
  # ------------------------------------------------------------
  # Helper: Lin CCC and error metrics
  # ------------------------------------------------------------
  ccc_lins <- function(obs, pred) {
    
    ok <- is.finite(obs) & is.finite(pred)
    obs <- obs[ok]
    pred <- pred[ok]
    
    if (length(obs) < 2L) {
      return(data.frame(
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
    
    data.frame(
      CCC = as.numeric(CCC),
      Cb = as.numeric(Cb),
      rho = as.numeric(rho),
      RMSE = as.numeric(RMSE),
      MAE = as.numeric(MAE)
    )
  }
  
  # ------------------------------------------------------------
  # Helper: project vector onto non-negative simplex
  # ------------------------------------------------------------
  project_simplex <- function(v) {
    
    v <- as.numeric(v)
    
    if (all(!is.finite(v))) {
      return(rep(1 / length(v), length(v)))
    }
    
    v[!is.finite(v)] <- 0
    
    n <- length(v)
    u <- sort(v, decreasing = TRUE)
    cssv <- cumsum(u)
    
    rho <- max(which(u + (1 - cssv) / seq_along(u) > 0))
    
    theta <- (cssv[rho] - 1) / rho
    
    w <- pmax(v - theta, 0)
    
    if (sum(w) <= 0) {
      w <- rep(1 / n, n)
    } else {
      w <- w / sum(w)
    }
    
    w
  }
  
  # ------------------------------------------------------------
  # Select models
  # ------------------------------------------------------------
  bestfit_use <- bestfit
  
  if (!is.null(model_ids)) {
    
    missing_models <- setdiff(model_ids, bestfit_use[[model_col]])
    
    if (length(missing_models) > 0) {
      stop(
        "The following `model_ids` were not found in `bestfit`: ",
        paste(missing_models, collapse = ", "),
        "."
      )
    }
    
    selected_models <- bestfit_use[bestfit_use[[model_col]] %in% model_ids, , drop = FALSE]
    
  } else {
    
    ord_metric <- bestfit_use[[weight_metric]]
    
    if (!hb) {
      ord_metric <- -ord_metric
    }
    
    ord <- order(-ord_metric, na.last = NA)
    
    selected_models <- bestfit_use[ord, , drop = FALSE]
    selected_models <- utils::head(selected_models, top_n)
  }
  
  selected_ids <- selected_models[[model_col]]
  selected_ids_chr <- as.character(selected_ids)
  
  if (length(selected_ids) == 0) {
    stop("No models were selected for the ensemble.")
  }
  
  if (verbose) {
    message(
      "Selected ", length(selected_ids), " model(s): ",
      paste(selected_ids, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------
  # Helper: compute model weights for non-stacked methods
  # ------------------------------------------------------------
  compute_model_weights <- function(selected_models, selected_ids, ensemble_method) {
    
    selected_ids_chr <- as.character(selected_ids)
    
    if (ensemble_method == "unweighted") {
      w <- rep(1 / length(selected_ids), length(selected_ids))
      names(w) <- selected_ids_chr
      return(w)
    }
    
    if (ensemble_method %in% c("best", "hard_voting")) {
      w <- rep(0, length(selected_ids))
      w[1] <- 1
      names(w) <- selected_ids_chr
      return(w)
    }
    
    metric <- selected_models[[weight_metric]]
    
    if (!hb) {
      metric <- -metric
    }
    
    if (weight_transform == "uniform") {
      w <- rep(1 / length(selected_ids), length(selected_ids))
      
    } else if (weight_transform == "positive") {
      
      metric_pos <- pmax(metric, 0)
      
      if (sum(metric_pos, na.rm = TRUE) <= 0) {
        w <- rep(1 / length(selected_ids), length(selected_ids))
      } else {
        w <- metric_pos / sum(metric_pos, na.rm = TRUE)
      }
      
    } else if (weight_transform == "rank_inverse") {
      
      ranks <- seq_along(selected_ids)
      w <- 1 / ranks
      w <- w / sum(w)
      
    } else if (weight_transform == "softmax") {
      
      metric_centered <- metric - max(metric, na.rm = TRUE)
      exp_metric <- exp(metric_centered)
      w <- exp_metric / sum(exp_metric, na.rm = TRUE)
      
    } else {
      stop("Unsupported `weight_transform`.")
    }
    
    names(w) <- selected_ids_chr
    w
  }
  
  # ------------------------------------------------------------
  # Helper: build prediction matrix
  # ------------------------------------------------------------
  build_prediction_matrix <- function(predictions, selected_ids) {
    
    pred_use <- predictions[predictions[[model_col]] %in% selected_ids, , drop = FALSE]
    
    if (nrow(pred_use) == 0) {
      stop("No predictions were found for the selected models.")
    }
    
    make_key <- function(dat, cols) {
      apply(dat[, cols, drop = FALSE], 1, paste, collapse = "||")
    }
    
    pred_use$.obs_key <- make_key(pred_use, id_cols)
    
    obs_df <- pred_use[, c(".obs_key", id_cols, response_col), drop = FALSE]
    obs_df <- obs_df[!duplicated(obs_df$.obs_key), , drop = FALSE]
    
    obs_keys <- obs_df$.obs_key
    model_names <- paste0("model_", selected_ids)
    
    pred_mat <- matrix(
      NA_real_,
      nrow = length(obs_keys),
      ncol = length(selected_ids),
      dimnames = list(obs_keys, model_names)
    )
    
    for (j in seq_along(selected_ids)) {
      
      mid <- selected_ids[j]
      
      tmp <- pred_use[pred_use[[model_col]] == mid, , drop = FALSE]
      
      key_j <- tmp$.obs_key
      val_j <- tmp[[prediction_col]]
      
      idx <- match(key_j, obs_keys)
      
      pred_mat[idx, j] <- val_j
    }
    
    complete_rows <- stats::complete.cases(pred_mat)
    
    if (sum(complete_rows) < 2L) {
      stop(
        "Fewer than two complete out-of-fold observations are available across selected models."
      )
    }
    
    pred_mat <- pred_mat[complete_rows, , drop = FALSE]
    obs_df <- obs_df[complete_rows, , drop = FALSE]
    obs <- obs_df[[response_col]]
    
    list(
      pred_mat = pred_mat,
      obs_df = obs_df,
      obs = obs
    )
  }
  
  # ------------------------------------------------------------
  # Helper: stacking weights
  # ------------------------------------------------------------
  compute_stacking <- function(pred_mat, obs) {
    
    X <- pred_mat
    y <- obs
    
    fit_stack <- function(X_train, y_train) {
      
      if (stack_intercept) {
        X_design <- cbind("(Intercept)" = 1, X_train)
      } else {
        X_design <- X_train
      }
      
      if (stacking_model == "lm" || stack_objective == "min_error") {
        
        fit <- stats::lm.fit(x = X_design, y = y_train)
        coef <- fit$coefficients
        coef[!is.finite(coef)] <- 0
        return(coef)
      }
      
      if (stacking_model == "ridge" || stack_objective == "regularized") {
        
        p <- ncol(X_design)
        penalty <- diag(lambda, p)
        
        if (stack_intercept) {
          penalty[1, 1] <- 0
        }
        
        XtX <- crossprod(X_design)
        Xty <- crossprod(X_design, y_train)
        
        coef <- tryCatch(
          as.numeric(solve(XtX + penalty, Xty)),
          error = function(e) rep(0, p)
        )
        
        return(coef)
      }
      
      stop("Unsupported stacking configuration.")
    }
    
    ensemble_pred <- numeric(nrow(X))
    
    for (i in seq_len(nrow(X))) {
      
      train_idx <- setdiff(seq_len(nrow(X)), i)
      
      coef_i <- fit_stack(
        X_train = X[train_idx, , drop = FALSE],
        y_train = y[train_idx]
      )
      
      if (stack_intercept) {
        x_i <- c(1, X[i, ])
      } else {
        x_i <- X[i, ]
      }
      
      ensemble_pred[i] <- sum(x_i * coef_i)
    }
    
    coef_final <- fit_stack(X_train = X, y_train = y)
    
    if (stack_intercept) {
      stack_intercept_value <- coef_final[1]
      stack_weights <- coef_final[-1]
    } else {
      stack_intercept_value <- 0
      stack_weights <- coef_final
    }
    
    if (stack_objective == "constrained") {
      stack_weights <- project_simplex(stack_weights)
      stack_intercept_value <- 0
      ensemble_pred <- as.numeric(X %*% stack_weights)
    }
    
    names(stack_weights) <- selected_ids_chr
    
    list(
      ensemble_pred = ensemble_pred,
      weights = stack_weights,
      intercept = stack_intercept_value
    )
  }
  
  # ------------------------------------------------------------
  # Load predictions when needed
  # ------------------------------------------------------------
  predictions_needed <- ensemble_scope %in% c("model", "both")
  
  if (predictions_needed) {
    
    if (is.null(predictions)) {
      predictions <- attr(bestfit, "predictions")
    }
    
    if (is.null(predictions) || !is.data.frame(predictions)) {
      stop(
        "`predictions` must be supplied or available as attr(bestfit, 'predictions') ",
        "when model-level ensemble is requested."
      )
    }
    
    required_pred <- c(model_col, id_cols, response_col, prediction_col)
    missing_pred <- setdiff(required_pred, names(predictions))
    
    if (length(missing_pred) > 0) {
      stop(
        "The following required column(s) are missing from `predictions`: ",
        paste(missing_pred, collapse = ", "),
        "."
      )
    }
  }
  
  # ------------------------------------------------------------
  # Model-level ensemble
  # ------------------------------------------------------------
  ensemble_summary <- NULL
  ensemble_predictions <- NULL
  prediction_weights <- NULL
  lag_weights <- NULL
  stack_intercepts <- NULL
  
  if (ensemble_scope %in% c("model", "both")) {
    
    pm <- build_prediction_matrix(predictions, selected_ids)
    
    pred_mat <- pm$pred_mat
    obs_df <- pm$obs_df
    obs <- pm$obs
    
    if (method == "stacked") {
      
      stack <- compute_stacking(pred_mat, obs)
      
      ensemble_pred <- stack$ensemble_pred
      prediction_weights <- stack$weights
      stack_intercepts <- stack$intercept
      
    } else {
      
      prediction_weights <- compute_model_weights(
        selected_models = selected_models,
        selected_ids = selected_ids,
        ensemble_method = method
      )
      
      ensemble_pred <- as.numeric(pred_mat %*% as.numeric(prediction_weights))
    }
    
    ensemble_predictions <- obs_df[, id_cols, drop = FALSE]
    ensemble_predictions$observed <- obs
    ensemble_predictions$predicted_ensemble <- ensemble_pred
    
    pred_base_df <- as.data.frame(pred_mat, stringsAsFactors = FALSE)
    ensemble_predictions <- cbind(ensemble_predictions, pred_base_df)
    
    ensemble_metrics <- ccc_lins(
      obs = ensemble_predictions$observed,
      pred = ensemble_predictions$predicted_ensemble
    )
    
    ensemble_summary <- data.frame(
      method = method,
      top_n = length(selected_ids),
      stacking_model = if (method == "stacked") stacking_model else NA_character_,
      stack_objective = if (method == "stacked") stack_objective else NA_character_,
      weight_metric = if (method == "weighted") weight_metric else NA_character_,
      weight_transform = if (method == "weighted") weight_transform else NA_character_,
      ensemble_metrics,
      stringsAsFactors = FALSE
    )
  }
  
  # ------------------------------------------------------------
  # Lag weights
  # ------------------------------------------------------------
  if (ensemble_scope %in% c("lag", "both")) {
    lag_weights <- compute_model_weights(
      selected_models = selected_models,
      selected_ids = selected_ids,
      ensemble_method = method
    )
  }
  
  # ------------------------------------------------------------
  # Automatic lag_data generation
  # ------------------------------------------------------------
  if (ensemble_scope %in% c("lag", "both") && is.null(lag_data)) {
    
    if (is.null(fit_list)) {
      fit_list <- attr(bestfit, "fits")
    }
    
    if (is.null(fit_list) || !is.list(fit_list) || length(fit_list) == 0) {
      stop(
        "`lag_data` was not supplied and no fitted models were found. ",
        "Use `find_bestfit(..., keep_fits = TRUE)` or provide `fit_list`."
      )
    }
    
    if (is.null(epi_data) || !is.data.frame(epi_data)) {
      stop(
        "`epi_data` is required to compute lag-based ensembles automatically. ",
        "Lag contributions depend on observed exposure profiles. ",
        "Provide `epi_data` or supply precomputed `lag_data`."
      )
    }
    
    if (is.null(group) || !is.character(group) || length(group) != 1L) {
      stop("`group` must be a single character string when computing lag_data automatically.")
    }
    
    if (!group %in% names(epi_data)) {
      stop("`group` ('", group, "') was not found in `epi_data`.")
    }
    
    lag_list <- list()
    lag_id <- 1L
    
    for (mid in selected_ids_chr) {
      
      fit_i <- fit_list[[mid]]
      
      if (is.null(fit_i)) {
        stop(
          "Could not find fitted model for model_id = ", mid,
          ". Check `attr(bestfit, 'fits')` or provide `fit_list`."
        )
      }
      
      vars_fit_i <- attr(fit_i, "epiexposure_vars")
      
      if (is.null(vars_fit_i) || length(vars_fit_i) == 0) {
        stop(
          "Fitted model for model_id = ", mid,
          " does not contain `epiexposure_vars` metadata."
        )
      }
      
      if (lg_strategy == "model") {
        
        var_i <- vars_fit_i
        
      } else if (lg_strategy == "requested_available") {
        
        if (is.null(var)) {
          var_i <- vars_fit_i
        } else {
          var_i <- intersect(var, vars_fit_i)
        }
        
        if (length(var_i) == 0) {
          if (verbose) {
            warning(
              "Skipping model_id = ", mid,
              " because none of the requested variables are present in this model. ",
              "Model variables: ", paste(vars_fit_i, collapse = ", "),
              call. = FALSE
            )
          }
          next
        }
        
      } else if (lg_strategy == "strict") {
        
        if (is.null(var)) {
          var_i <- vars_fit_i
        } else {
          missing_vars_i <- setdiff(var, vars_fit_i)
          
          if (length(missing_vars_i) > 0) {
            stop(
              "For model_id = ", mid,
              ", `var` contains name(s) not found in the fitted model: ",
              paste(missing_vars_i, collapse = ", "),
              ". Use exactly the same variable name(s) stored in the model: ",
              paste(vars_fit_i, collapse = ", "),
              "."
            )
          }
          
          var_i <- var
        }
      }
      
      args_i <- c(
        list(
          epi_data = epi_data,
          group = group,
          fit = fit_i,
          var = var_i,
          reverse = reverse
        ),
        compute_ecilag_args
      )
      
      lag_res <- tryCatch(
        do.call(compute_ecilag, args_i),
        error = function(e) {
          stop(
            "Automatic lag decomposition failed for model_id = ", mid, ": ",
            conditionMessage(e)
          )
        }
      )
      
      if (!is.list(lag_res) || is.null(lag_res$by_lag) || !is.data.frame(lag_res$by_lag)) {
        stop(
          "`compute_ecilag()` did not return a valid `by_lag` data.frame ",
          "for model_id = ", mid, "."
        )
      }
      
      tmp <- lag_res$by_lag
      tmp[[model_col]] <- suppressWarnings(type.convert(mid, as.is = TRUE))
      
      lag_list[[lag_id]] <- tmp
      lag_id <- lag_id + 1L
    }
    
    if (length(lag_list) == 0) {
      stop(
        "No lag decompositions were generated. This may happen if none of the ",
        "selected models contain the requested variable(s)."
      )
    }
    
    lag_data <- do.call(rbind, lag_list)
    
    if (verbose) {
      message("Lag data were generated automatically using stored fitted models.")
    }
  }
  
  # ------------------------------------------------------------
  # Validate lag_data after auto-generation or user supply
  # ------------------------------------------------------------
  if (ensemble_scope %in% c("lag", "both")) {
    
    if (is.null(lag_data) || !is.data.frame(lag_data)) {
      stop("`lag_data` must be a data.frame.")
    }
    
    required_lag <- c(model_col, lag_col, lag_contribution_col)
    missing_lag <- setdiff(required_lag, names(lag_data))
    
    if (length(missing_lag) > 0) {
      stop(
        "The following required column(s) are missing from `lag_data`: ",
        paste(missing_lag, collapse = ", "),
        "."
      )
    }
    
    missing_lag_group_cols <- setdiff(lag_group_cols, names(lag_data))
    
    if (length(missing_lag_group_cols) > 0) {
      stop(
        "The following `lag_group_cols` are missing from `lag_data`: ",
        paste(missing_lag_group_cols, collapse = ", "),
        "."
      )
    }
  }
  
  # ------------------------------------------------------------
  # Lag-level ensemble
  # ------------------------------------------------------------
  ensemble_by_lag <- NULL
  
  if (ensemble_scope %in% c("lag", "both")) {
    
    lag_use <- lag_data[as.character(lag_data[[model_col]]) %in% selected_ids_chr, , drop = FALSE]
    
    if (nrow(lag_use) == 0) {
      stop("No lag contributions were found for the selected models.")
    }
    
    lag_use$.model_weight <- as.numeric(
      lag_weights[match(as.character(lag_use[[model_col]]), names(lag_weights))]
    )
    
    if (any(!is.finite(lag_use$.model_weight))) {
      stop("Could not match model weights to `lag_data`.")
    }
    
    split_key <- apply(lag_use[, lag_group_cols, drop = FALSE], 1, paste, collapse = "||")
    split_ids <- unique(split_key)
    
    lag_out <- lapply(split_ids, function(key_i) {
      
      tmp <- lag_use[split_key == key_i, , drop = FALSE]
      
      base <- tmp[1, lag_group_cols, drop = FALSE]
      
      w_tmp <- tmp$.model_weight
      
      if (isTRUE(rl_weights)) {
        if (sum(w_tmp, na.rm = TRUE) > 0) {
          w_tmp <- w_tmp / sum(w_tmp, na.rm = TRUE)
        }
      }
      
      contribution_ens <- sum(
        tmp[[lag_contribution_col]] * w_tmp,
        na.rm = TRUE
      )
      
      out <- base
      out$ensemble_contribution <- contribution_ens
      
      if (!is.null(lag_weight_col) && lag_weight_col %in% names(tmp)) {
        out$ensemble_weight <- sum(
          tmp[[lag_weight_col]] * w_tmp,
          na.rm = TRUE
        )
      }
      
      out$n_models <- length(unique(tmp[[model_col]]))
      out$models_used <- paste(sort(unique(tmp[[model_col]])), collapse = ", ")
      
      out
    })
    
    ensemble_by_lag <- do.call(rbind, lag_out)
    
    group_for_percent <- setdiff(lag_group_cols, lag_col)
    
    if (length(group_for_percent) == 0) {
      
      denom <- sum(abs(ensemble_by_lag$ensemble_contribution), na.rm = TRUE)
      
      ensemble_by_lag$ensemble_percent_contribution <- if (denom > 0) {
        100 * abs(ensemble_by_lag$ensemble_contribution) / denom
      } else {
        NA_real_
      }
      
    } else {
      
      group_key <- apply(
        ensemble_by_lag[, group_for_percent, drop = FALSE],
        1,
        paste,
        collapse = "||"
      )
      
      group_ids <- unique(group_key)
      ensemble_by_lag$ensemble_percent_contribution <- NA_real_
      
      for (g in group_ids) {
        idx <- which(group_key == g)
        denom <- sum(abs(ensemble_by_lag$ensemble_contribution[idx]), na.rm = TRUE)
        
        ensemble_by_lag$ensemble_percent_contribution[idx] <- if (denom > 0) {
          100 * abs(ensemble_by_lag$ensemble_contribution[idx]) / denom
        } else {
          NA_real_
        }
      }
    }
  }
  
  # ------------------------------------------------------------
  # Model weights table
  # ------------------------------------------------------------
  model_weights_df <- data.frame(
    model_id = selected_ids,
    prediction_weight = if (!is.null(prediction_weights)) {
      as.numeric(prediction_weights[match(selected_ids_chr, names(prediction_weights))])
    } else {
      NA_real_
    },
    lag_weight = if (!is.null(lag_weights)) {
      as.numeric(lag_weights[match(selected_ids_chr, names(lag_weights))])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
  
  if (!is.null(stack_intercepts)) {
    attr(model_weights_df, "stack_intercept") <- stack_intercepts
  }
  
  # ------------------------------------------------------------
  # Return
  # ------------------------------------------------------------
  out <- list(
    ensemble_summary = ensemble_summary,
    ensemble_predictions = ensemble_predictions,
    ensemble_by_lag = ensemble_by_lag,
    lag_data = lag_data,
    model_weights = model_weights_df,
    selected_models = selected_models,
    method = method,
    ensemble_scope = ensemble_scope,
    lg_strategy = lg_strategy,
    rl_weights = rl_weights
  )
  
  out
}