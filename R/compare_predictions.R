#' Compare predicted outcomes between two exposure scenarios
#'
#' @param fit Fitted model (output of fit_epidlnm())
#' @param profiles1 Exposure profile(s): numeric vector or named list of vectors
#' @param profiles2 Exposure profile(s): numeric vector or named list of vectors
#' @param re "population" or "conditional" (passed to predict_outcome)
#' @param id Optional vector of IDs when re="conditional"
#' @param allow_new_levels Passed to predict_outcome
#' @param type Prediction scale passed to predict_outcome
#' @param eps Small constant to avoid division by zero in percent/ratio
#'
#' @return data.frame with pred1, pred2, diff, percent_change, ratio (and id if provided)
#' @export
compare_predictions <- function(
    fit,
    profiles1,
    profiles2,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link", "conditional"),
    eps = 1e-12
) {
  re <- match.arg(re)
  type <- match.arg(type)
  
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) stop("eps must be a positive numeric scalar.")
  
  p1 <- predict_outcome(
    fit = fit,
    profiles = profiles1,
    re = re,
    id = id,
    allow_new_levels = allow_new_levels,
    type = type
  )
  
  p2 <- predict_outcome(
    fit = fit,
    profiles = profiles2,
    re = re,
    id = id,
    allow_new_levels = allow_new_levels,
    type = type
  )
  
  # Align outputs (by id if present)
  common_cols <- intersect(names(p1), names(p2))
  id_col <- setdiff(common_cols, "prediction")
  
  if (length(id_col) == 1) {
    # join by id
    out <- merge(
      p1, p2,
      by = id_col,
      suffixes = c("_1", "_2"),
      all = FALSE
    )
    names(out)[names(out) == "prediction_1"] <- "pred1"
    names(out)[names(out) == "prediction_2"] <- "pred2"
  } else {
    # population output
    out <- data.frame(pred1 = p1$prediction, pred2 = p2$prediction)
  }
  
  out$diff <- out$pred2 - out$pred1
  out$percent_change <- (out$diff) / (abs(out$pred1) + eps) * 100
  out$ratio <- out$pred2 / (out$pred1 + eps)
  
  out
}