#' Compare predicted outcomes between multiple exposure scenarios
#'
#' Compares predicted outcomes across two or more exposure scenarios,
#' using `predict_outcome()` as the computational backend.
#'
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
    reverse = FALSE,  # pass-through
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000,
    eps = 1e-12
) {

  re <- match.arg(re)
  type <- match.arg(type)
  output <- match.arg(output)

  if (!is.numeric(eps) || eps <= 0) {
    stop("eps must be positive.")
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  # -----------------------
  # BACKWARD COMPAT
  # -----------------------
  if (is.null(profiles)) {

    if (is.null(profiles1) || is.null(profiles2)) {
      stop("Provide either `profiles` OR both `profiles1` and `profiles2`.")
    }

    profiles <- list(
      scenario1 = profiles1,
      scenario2 = profiles2
    )
  }

  if (!is.list(profiles) || is.null(names(profiles))) {
    stop("`profiles` must be a named list.")
  }

  if (length(profiles) < 2) {
    stop("At least two profiles are required.")
  }

  # ------------------------------------------------------------
  # ALWAYS USE SAMPLES IF UNCERTAINTY
  # ------------------------------------------------------------
  use_output <- if (uncertainty) "samples" else "summary"

  preds_list <- lapply(names(profiles), function(nm) {

    p <- predict_outcome(
      fit = fit,
      profiles = profiles[[nm]],
      re = re,
      id = id,
      allow_new_levels = allow_new_levels,
      type = type,
      reverse = reverse,   # ✅ AQUI está o pass-through correto
      uncertainty = uncertainty,
      output = use_output,
      n_samples = n_samples
    )

    p$scenario <- nm
    p
  })

  preds <- do.call(rbind, preds_list)

  # ============================================================
  # ✅ NO UNCERTAINTY
  # ============================================================
  if (!uncertainty) {

    idvar_use <- if (!is.null(id)) id else NULL

    wide <- reshape(
      preds,
      idvar = idvar_use,
      timevar = "scenario",
      direction = "wide"
    )

    combs <- combn(names(profiles), 2, simplify = FALSE)

    out_list <- lapply(combs, function(cb) {

      s1 <- cb[1]
      s2 <- cb[2]

      p1 <- wide[[paste0("prediction.", s1)]]
      p2 <- wide[[paste0("prediction.", s2)]]

      df <- data.frame(
        scenario1 = s1,
        scenario2 = s2,
        pred1 = p1,
        pred2 = p2
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

  # ============================================================
  # ✅ UNCERTAINTY (SAMPLE-BASED)
  # ============================================================
  combs <- combn(names(profiles), 2, simplify = FALSE)

  out_samples <- lapply(combs, function(cb) {

    s1 <- cb[1]
    s2 <- cb[2]

    p1 <- preds[preds$scenario == s1, ]
    p2 <- preds[preds$scenario == s2, ]

    by_cols <- "sample"
    if (!is.null(id)) {
      by_cols <- c(id, "sample")
    }

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

  samples_df <- do.call(rbind, out_samples)

  if (output == "samples") {
    return(samples_df)
  }

  # ============================================================
  # ✅ SUMMARY (MEDIAN-BASED)
  # ============================================================
  group_vars <- intersect(c(id, "scenario1", "scenario2"), names(samples_df))

  summary_df <- samples_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
    dplyr::summarise(
      pred1 = stats::median(.data$pred1, na.rm = TRUE),
      pred2 = stats::median(.data$pred2, na.rm = TRUE),

      diff = stats::median(.data$diff, na.rm = TRUE),
      diff_sd = stats::sd(.data$diff, na.rm = TRUE),
      diff_lower = stats::quantile(.data$diff, 0.025, na.rm = TRUE, names = FALSE),
      diff_upper = stats::quantile(.data$diff, 0.975, na.rm = TRUE, names = FALSE),

      percent_change = stats::median(.data$percent_change, na.rm = TRUE),
      ratio = stats::median(.data$ratio, na.rm = TRUE),

      .groups = "drop"
    )

  as.data.frame(summary_df)
}
