#' Compare predicted outcomes between multiple exposure scenarios
#'
#' @param fit Fitted model (output of fit_epidlnm())
#' @param profiles Can be:
#'   - list(profile1 = ..., profile2 = ...) (new recommended)
#'   - OR profiles1 / profiles2 (for backward compatibility)
#' @param profiles1 (legacy) first profile
#' @param profiles2 (legacy) second profile
#' @param re "population" or "conditional"
#' @param id Optional vector of IDs
#' @param allow_new_levels Passed to predict_outcome
#' @param type Prediction scale
#' @param uncertainty Logical; if TRUE, propagate uncertainty
#' @param output "summary" or "samples"
#' @param n_samples Number of samples
#' @param eps Small constant for stability
#'
#' @return data.frame
#' @export
compare_predictions <- function(
    fit,
    profiles = NULL,
    profiles1 = NULL,
    profiles2 = NULL,
    re = c("population", "conditional"),
    id = NULL,
    allow_new_levels = FALSE,
    type = c("response", "link", "conditional"),
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    eps = 1e-12
) {

  re <- match.arg(re)
  type <- match.arg(type)
  output <- match.arg(output)

  if (!is.numeric(eps) || eps <= 0) stop("eps must be positive.")

  # -----------------------
  # ✅ BACKWARD COMPATIBILITY
  # -----------------------
  if (is.null(profiles)) {
    if (is.null(profiles1) || is.null(profiles2)) {
      stop("Provide either `profiles` (list) OR both `profiles1` and `profiles2`.")
    }

    profiles <- list(
      scenario1 = profiles1,
      scenario2 = profiles2
    )
  }

  if (!is.list(profiles) || is.null(names(profiles))) {
    stop("`profiles` must be a named list.")
  }

  n_prof <- length(profiles)
  if (n_prof < 2) {
    stop("At least two profiles are required.")
  }

  # -----------------------
  # Run predictions
  # -----------------------
  preds_list <- lapply(names(profiles), function(nm) {

    p <- predict_outcome(
      fit = fit,
      profiles = profiles[[nm]],
      re = re,
      id = id,
      allow_new_levels = allow_new_levels,
      type = type,
      uncertainty = uncertainty,
      output = output,
      n_samples = n_samples
    )

    p$scenario <- nm
    p
  })

  preds <- do.call(rbind, preds_list)

  # -----------------------
  # ✅ NO UNCERTAINTY
  # -----------------------
  if (!uncertainty) {

    # reshape wide
    if (!is.null(id)) {

      wide <- reshape(
        preds,
        idvar = id,
        timevar = "scenario",
        direction = "wide"
      )

    } else {

      wide <- reshape(
        preds,
        idvar = NULL,
        timevar = "scenario",
        direction = "wide"
      )
    }

    # pairwise comparisons
    combs <- combn(names(profiles), 2, simplify = FALSE)

    out_list <- lapply(combs, function(cb) {

      s1 <- cb[1]
      s2 <- cb[2]

      pred1 <- wide[[paste0("prediction.", s1)]]
      pred2 <- wide[[paste0("prediction.", s2)]]

      df <- data.frame(
        scenario1 = s1,
        scenario2 = s2,
        pred1 = pred1,
        pred2 = pred2
      )

      df$diff <- df$pred2 - df$pred1
      df$percent_change <- df$diff / (abs(df$pred1) + eps) * 100
      df$ratio <- df$pred2 / (df$pred1 + eps)

      if (!is.null(id)) {
        df[[id]] <- wide[[id]]
        df <- df[, c(id, "scenario1", "scenario2", "pred1", "pred2",
                     "diff", "percent_change", "ratio")]
      }

      df
    })

    return(do.call(rbind, out_list))
  }

  # -----------------------
  # ✅ UNCERTAINTY: SAMPLES
  # -----------------------
  if (output == "samples") {

    # ensure join by sample (and id if exists)
    merge_cols <- c("scenario", "sample")
    if (!is.null(id)) merge_cols <- c(id, merge_cols)

    combs <- combn(names(profiles), 2, simplify = FALSE)

    out_list <- lapply(combs, function(cb) {

      s1 <- cb[1]
      s2 <- cb[2]

      p1 <- preds[preds$scenario == s1, ]
      p2 <- preds[preds$scenario == s2, ]

      by_cols <- "sample"
      if (!is.null(id)) by_cols <- c(id, "sample")

      m <- merge(p1, p2, by = by_cols, suffixes = c("_1", "_2"))

      names(m)[names(m) == "prediction_1"] <- "pred1"
      names(m)[names(m) == "prediction_2"] <- "pred2"

      m$diff <- m$pred2 - m$pred1
      m$percent_change <- m$diff / (abs(m$pred1) + eps) * 100
      m$ratio <- m$pred2 / (m$pred1 + eps)

      m$scenario1 <- s1
      m$scenario2 <- s2

      m
    })

    return(do.call(rbind, out_list))
  }

  # -----------------------
  # ✅ UNCERTAINTY: SUMMARY
  # -----------------------
  combs <- combn(names(profiles), 2, simplify = FALSE)

  out_list <- lapply(combs, function(cb) {

    s1 <- cb[1]
    s2 <- cb[2]

    p1 <- preds[preds$scenario == s1, ]
    p2 <- preds[preds$scenario == s2, ]

    df <- data.frame(
      scenario1 = s1,
      scenario2 = s2,
      pred1 = p1$prediction,
      pred2 = p2$prediction,
      sd1 = p1$sd,
      sd2 = p2$sd
    )

    df$diff <- df$pred2 - df$pred1
    df$sd_diff <- sqrt(df$sd1^2 + df$sd2^2)

    df$percent_change <- df$diff / (abs(df$pred1) + eps) * 100
    df$ratio <- df$pred2 / (df$pred1 + eps)

    if (!is.null(id)) {
      df[[id]] <- p1[[id]]
      df <- df[, c(id, "scenario1", "scenario2", "pred1", "pred2",
                   "sd1", "sd2", "sd_diff",
                   "diff", "percent_change", "ratio")]
    }

    df
  })

  return(do.call(rbind, out_list))
}
