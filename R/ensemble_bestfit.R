#' Create ensembles from best-fit DLNM candidate models
#'
#' Builds ensemble predictions from the out-of-fold predictions produced by
#' `find_bestfit()`. The function supports several ensemble strategies:
#'
#' \itemize{
#'   \item `"best"`: uses the single best-ranked model.
#'   \item `"hard_voting"`: alias for `"best"` for continuous responses.
#'   \item `"unweighted"`: averages predictions from selected models with equal weights.
#'   \item `"weighted"`: averages predictions using model-performance weights.
#'   \item `"stacked"`: learns ensemble weights from out-of-fold predictions using
#'   a linear or ridge meta-model.
#' }
#'
#' @param bestfit A data.frame returned by `find_bestfit()`. It must contain
#'   model-level performance metrics and should have an attribute named
#'   `"predictions"` containing fold-level predictions.
#'
#' @param predictions Optional data.frame of fold-level predictions. If `NULL`,
#'   the function uses `attr(bestfit, "predictions")`.
#'
#' @param method Character. Ensemble strategy. Options are `"best"`,
#'   `"hard_voting"`, `"unweighted"`, `"weighted"`, and `"stacked"`.
#'
#' @param top_n Integer. Number of top-ranked models to include in the ensemble.
#'   Ignored when `model_ids` is supplied. Default is `10`.
#'
#' @param model_ids Optional vector of model IDs to use in the ensemble. If
#'   supplied, `top_n` is ignored.
#'
#' @param weight_metric Character. Column in `bestfit` used to compute weights
#'   for `method = "weighted"`. Default is `"CCC"`.
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
#' @param lambda Numeric. Ridge penalty used when `stacking_model = "ridge"`.
#'   Default is `1e-6`.
#'
#' @param stack_intercept Logical. If `TRUE`, includes an intercept in the
#'   stacking model. Default is `TRUE`.
#'
#' @param response_col Character. Name of the observed-response column in
#'   `predictions`. Default is `"observed"`.
#'
#' @param prediction_col Character. Name of the predicted-response column in
#'   `predictions`. Default is `"predicted"`.
#'
#' @param model_col Character. Name of the model ID column in both `bestfit`
#'   and `predictions`. Default is `"model_id"`.
#'
#' @param id_cols Character vector identifying each out-of-fold prediction row.
#'   By default, the function uses `c("group", "fold")`, which matches the
#'   output structure suggested for `find_bestfit()`.
#'
#' @param higher_is_better Logical. If `TRUE`, larger values of `weight_metric`
#'   indicate better models. Default is `TRUE`.
#'
#' @param verbose Logical. If `TRUE`, prints informative messages.
#'
#' @return A list with:
#' \itemize{
#'   \item `ensemble_summary`: performance metrics for the ensemble.
#'   \item `ensemble_predictions`: out-of-fold observed and ensemble-predicted values.
#'   \item `model_weights`: weights assigned to selected models.
#'   \item `selected_models`: subset of `bestfit` used in the ensemble.
#'   \item `method`: ensemble method used.
#' }
#'
#' @details
#' The function uses out-of-fold predictions from `find_bestfit()`. This is
#' important because ensemble performance should be evaluated using predictions
#' for observations that were not used to fit the corresponding base model.
#'
#' For `method = "stacked"`, the function fits a meta-model using the matrix of
#' out-of-fold predictions. To reduce optimistic bias, the reported stacked
#' predictions are generated using leave-one-out fitting at the meta-model level.
#'
#' For continuous responses, classical hard voting is not directly applicable.
#' Therefore, `method = "hard_voting"` is implemented as hard model selection,
#' equivalent to using the single best-ranked model.
#'
#' @export
ensemble_bestfit <- function(
    bestfit,
    predictions = NULL,
    method = c("unweighted", "weighted", "hard_voting", "best", "stacked"),
    top_n = 10,
    model_ids = NULL,
    weight_metric = "CCC",
    weight_transform = c("softmax", "positive", "rank_inverse", "uniform"),
    stacking_model = c("ridge", "lm"),
    lambda = 1e-6,
    stack_intercept = TRUE,
    response_col = "observed",
    prediction_col = "predicted",
    model_col = "model_id",
    id_cols = c("group", "fold"),
    higher_is_better = TRUE,
    verbose = TRUE
) {
  
  method <- match.arg(method)
  weight_transform <- match.arg(weight_transform)
  stacking_model <- match.arg(stacking_model)
  
  # ------------------------------------------------------------
  # Validations
  # ------------------------------------------------------------
  if (!is.data.frame(bestfit)) {
    stop("`bestfit` must be a data.frame returned by `find_bestfit()`.")
  }
  
  if (is.null(predictions)) {
    predictions <- attr(bestfit, "predictions")
  }
  
  if (is.null(predictions) || !is.data.frame(predictions)) {
    stop(
      "`predictions` must be supplied or available as attr(bestfit, 'predictions')."
    )
  }
  
  required_bestfit <- c(model_col, weight_metric)
  
  missing_bestfit <- setdiff(required_bestfit, names(bestfit))
  if (length(missing_bestfit) > 0) {
    stop(
      "The following required column(s) are missing from `bestfit`: ",
      paste(missing_bestfit, collapse = ", "),
      "."
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
  
  if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) || top_n < 1) {
    stop("`top_n` must be a positive integer.")
  }
  
  top_n <- as.integer(top_n)
  
  if (!is.null(model_ids)) {
    missing_models <- setdiff(model_ids, bestfit[[model_col]])
    if (length(missing_models) > 0) {
      stop(
        "The following `model_ids` were not found in `bestfit`: ",
        paste(missing_models, collapse = ", "),
        "."
      )
    }
  }
  
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be a non-negative numeric scalar.")
  }
  
  if (!is.logical(stack_intercept) || length(stack_intercept) != 1L) {
    stop("`stack_intercept` must be TRUE or FALSE.")
  }
  
  if (!is.logical(higher_is_better) || length(higher_is_better) != 1L) {
    stop("`higher_is_better` must be TRUE or FALSE.")
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
  # Select models
  # ------------------------------------------------------------
  bestfit_use <- bestfit
  
  if (!is.null(model_ids)) {
    
    selected_models <- bestfit_use[bestfit_use[[model_col]] %in% model_ids, , drop = FALSE]
    
  } else {
    
    ord_metric <- bestfit_use[[weight_metric]]
    
    if (!higher_is_better) {
      ord_metric <- -ord_metric
    }
    
    ord <- order(-ord_metric, na.last = NA)
    
    selected_models <- bestfit_use[ord, , drop = FALSE]
    selected_models <- utils::head(selected_models, top_n)
  }
  
  selected_ids <- selected_models[[model_col]]
  
  if (length(selected_ids) == 0) {
    stop("No models were selected for the ensemble.")
  }
  
  if (verbose) {
    message(
      "Selected ", length(selected_ids), " model(s): ",
      paste(selected_ids, collapse = ", ")
    )
  }
  
  pred_use <- predictions[predictions[[model_col]] %in% selected_ids, , drop = FALSE]
  
  if (nrow(pred_use) == 0) {
    stop("No predictions were found for the selected models.")
  }
  
  # ------------------------------------------------------------
  # Build observation key
  # ------------------------------------------------------------
  make_key <- function(dat, cols) {
    apply(dat[, cols, drop = FALSE], 1, paste, collapse = "||")
  }
  
  pred_use$.obs_key <- make_key(pred_use, id_cols)
  
  obs_df <- pred_use[, c(".obs_key", id_cols, response_col), drop = FALSE]
  obs_df <- obs_df[!duplicated(obs_df$.obs_key), , drop = FALSE]
  
  # ------------------------------------------------------------
  # Build prediction matrix manually
  # ------------------------------------------------------------
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
  
  # ------------------------------------------------------------
  # Compute base model weights
  # ------------------------------------------------------------
  compute_model_weights <- function(selected_models, selected_ids) {
    
    if (method %in% c("unweighted")) {
      w <- rep(1 / length(selected_ids), length(selected_ids))
      names(w) <- selected_ids
      return(w)
    }
    
    if (method %in% c("best", "hard_voting")) {
      w <- rep(0, length(selected_ids))
      w[1] <- 1
      names(w) <- selected_ids
      return(w)
    }
    
    metric <- selected_models[[weight_metric]]
    
    if (!higher_is_better) {
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
    
    names(w) <- selected_ids
    w
  }
  
  # ------------------------------------------------------------
  # Ensemble prediction
  # ------------------------------------------------------------
  model_weights <- NULL
  ensemble_pred <- NULL
  stack_weights <- NULL
  stack_intercepts <- NULL
  
  if (method %in% c("unweighted", "weighted", "best", "hard_voting")) {
    
    if (method == "weighted") {
      model_weights <- compute_model_weights(selected_models, selected_ids)
    } else {
      model_weights <- compute_model_weights(selected_models, selected_ids)
    }
    
    w <- as.numeric(model_weights)
    ensemble_pred <- as.numeric(pred_mat %*% w)
    
  } else if (method == "stacked") {
    
    X <- pred_mat
    y <- obs
    
    fit_stack <- function(X_train, y_train) {
      
      if (stack_intercept) {
        X_design <- cbind("(Intercept)" = 1, X_train)
      } else {
        X_design <- X_train
      }
      
      if (stacking_model == "lm") {
        
        fit <- stats::lm.fit(x = X_design, y = y_train)
        coef <- fit$coefficients
        coef[!is.finite(coef)] <- 0
        return(coef)
      }
      
      if (stacking_model == "ridge") {
        
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
      
      stop("Unsupported `stacking_model`.")
    }
    
    # LOOCV meta-predictions
    ensemble_pred <- numeric(nrow(X))
    coef_list <- vector("list", nrow(X))
    
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
      coef_list[[i]] <- coef_i
    }
    
    # final stack weights from all OOF rows
    coef_final <- fit_stack(X_train = X, y_train = y)
    
    if (stack_intercept) {
      stack_intercepts <- coef_final[1]
      stack_weights <- coef_final[-1]
    } else {
      stack_intercepts <- 0
      stack_weights <- coef_final
    }
    
    names(stack_weights) <- selected_ids
    
    model_weights <- stack_weights
  }
  
  # ------------------------------------------------------------
  # Output predictions
  # ------------------------------------------------------------
  ensemble_predictions <- obs_df[, id_cols, drop = FALSE]
  ensemble_predictions$observed <- obs
  ensemble_predictions$predicted_ensemble <- ensemble_pred
  
  # Also keep base model predictions
  pred_base_df <- as.data.frame(pred_mat, stringsAsFactors = FALSE)
  ensemble_predictions <- cbind(ensemble_predictions, pred_base_df)
  
  # ------------------------------------------------------------
  # Metrics
  # ------------------------------------------------------------
  ensemble_metrics <- ccc_lins(
    obs = ensemble_predictions$observed,
    pred = ensemble_predictions$predicted_ensemble
  )
  
  ensemble_summary <- data.frame(
    method = method,
    top_n = length(selected_ids),
    stacking_model = if (method == "stacked") stacking_model else NA_character_,
    weight_metric = if (method == "weighted") weight_metric else NA_character_,
    weight_transform = if (method == "weighted") weight_transform else NA_character_,
    ensemble_metrics,
    stringsAsFactors = FALSE
  )
  
  # ------------------------------------------------------------
  # Model weights table
  # ------------------------------------------------------------
  model_weights_df <- data.frame(
    model_id = selected_ids,
    weight = as.numeric(model_weights),
    stringsAsFactors = FALSE
  )
  
  if (method == "stacked") {
    attr(model_weights_df, "stack_intercept") <- stack_intercepts
  }
  
  selected_models_out <- selected_models
  
  # ------------------------------------------------------------
  # Return
  # ------------------------------------------------------------
  out <- list(
    ensemble_summary = ensemble_summary,
    ensemble_predictions = ensemble_predictions,
    model_weights = model_weights_df,
    selected_models = selected_models_out,
    method = method
  )
  
  out
}