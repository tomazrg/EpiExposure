#' Simulate epidemiological DLNM scenarios
#'
#' Builds exposure profiles over lag periods and predicts outcomes under
#' different environmental scenarios.
#'
#' @param fit Fitted model (from fit_epidlnm)
#' @param dat Original design matrix
#' @param cb_templates List of crossbasis templates (one per variable)
#' @param wx_ref Reference weather data
#' @param lag_windows Output from define_lag_windows()
#' @param scenarios Output from simulate_range() OR named list of scenarios
#' @param ref_vals Optional reference values (default = median)
#' @param pop_level Logical for population-level prediction
#'
#' @return data.frame (scenario or profile-ready)
#'
#' @export
simulate_scenarios <- function(fit,
                               dat,
                               cb_templates,
                               wx_ref,
                               lag_windows,
                               scenarios,
                               ref_vals = NULL,
                               pop_level = TRUE) {

  # ------------------------------------------------------------
  # ✅ Handle structured input from simulate_range()
  # ------------------------------------------------------------
  if (is.list(scenarios) && "scenarios" %in% names(scenarios)) {
    scenario_info <- scenarios$info
    scenarios     <- scenarios$scenarios
  } else {
    scenario_info <- NULL
  }

  # ------------------------------------------------------------
  # ✅ Extract coefficients (ALL engines)
  # ------------------------------------------------------------
  extract_beta <- function(fit) {

    if (inherits(fit, "glmmTMB")) return(fixef(fit)$cond)

    if (inherits(fit, "brmsfit")) {
      fe <- brms::fixef(fit)
      beta <- fe[, "Estimate"]
      names(beta) <- rownames(fe)
      return(beta)
    }

    if (inherits(fit, "inla")) {
      return(fit$summary.fixed$mean)
    }

    if (inherits(fit, "HLfit")) {
      return(spaMM::fixef(fit))
    }

    coef(fit)
  }

  beta <- extract_beta(fit)
  vars <- names(cb_templates)

  # ------------------------------------------------------------
  # ✅ Reference values
  # ------------------------------------------------------------
  if (is.null(ref_vals)) {
    ref_vals <- lapply(vars, function(v) {
      as.numeric(stats::median(wx_ref[[v]], na.rm = TRUE))
    })
    names(ref_vals) <- vars
  }

  # ------------------------------------------------------------
  # Lag settings
  # ------------------------------------------------------------
  lag_max <- attr(cb_templates[[1]], "lag")
  N <- lag_max + 1

  lag_to_idx <- function(lags, N) {
    idx <- N - lags
    idx[idx >= 1 & idx <= N]
  }

  # ------------------------------------------------------------
  # ✅ Build profile for ONE variable
  # ------------------------------------------------------------
  build_profile <- function(var, scen_values) {

    x <- rep(ref_vals[[var]], N)

    for (i in seq_len(nrow(lag_windows))) {

      w_id <- lag_windows$window_id[i]

      if (!is.null(scen_values[[w_id]][[var]])) {

        lags <- seq(lag_windows$lag_start[i],
                    lag_windows$lag_end[i])

        idx <- lag_to_idx(lags, N)

        x[idx] <- scen_values[[w_id]][[var]]
      }
    }

    x
  }

  # ------------------------------------------------------------
  # ✅ Extract CB row
  # ------------------------------------------------------------
  extract_cb <- function(x, cb_template) {

    cb <- dlnm::crossbasis(
      x,
      lag    = lag_max,
      argvar = attr(cb_template, "argvar"),
      arglag = attr(cb_template, "arglag")
    )

    as.numeric(cb[length(x), ])
  }

  get_cb_names <- function(dat, var) {
    grep(paste0("^cb_", var, "_"), names(dat), value = TRUE)
  }

  # ------------------------------------------------------------
  # ✅ SINGLE prediction (no recursion bug ✅)
  # ------------------------------------------------------------
  predict_single <- function(scen_values) {

    nd <- dat[1, , drop = FALSE]

    for (v in vars) {

      profile <- build_profile(v, scen_values)
      cb_row  <- extract_cb(profile, cb_templates[[v]])

      cols <- get_cb_names(dat, v)

      nd[cols] <- as.list(cb_row)
    }

    if (pop_level && inherits(fit, "glmmTMB")) {
      pred <- predict(fit, newdata = nd, type = "response", re.form = NA)
    } else {
      pred <- predict(fit, newdata = nd, type = "response")
    }

    as.numeric(pred)
  }

  # ------------------------------------------------------------
  # ✅ Predict scenario (handles vector/profile correctly)
  # ------------------------------------------------------------
  predict_scenario <- function(scen_values) {

    # detectar variável contínua (profile)
    lens <- sapply(scen_values[[1]], length)
    var_cont <- names(lens[lens > 1])

    # ------------------------------------------------------------
    # ✅ PROFILE MODE
    # ------------------------------------------------------------
    if (length(var_cont) > 0) {

      n_pts <- length(scen_values[[1]][[var_cont]])

      preds <- numeric(n_pts)

      for (i in seq_len(n_pts)) {

        scen_i <- lapply(scen_values, function(w) {
          lapply(w, function(val) {
            if (length(val) > 1) val[i] else val
          })
        })

        preds[i] <- predict_single(scen_i)
      }

      return(preds)
    }

    # ------------------------------------------------------------
    # ✅ STANDARD scenario
    # ------------------------------------------------------------
    predict_single(scen_values)
  }

  # ------------------------------------------------------------
  # ✅ Apply scenarios
  # ------------------------------------------------------------
  out <- purrr::imap_dfr(scenarios, function(scen, name) {

    preds <- predict_scenario(scen)

    # ----------------------------------------
    # ✅ scenario mode (no info_df)
    # ----------------------------------------
    if (is.null(scenario_info)) {

      return(data.frame(
        scenario   = name,
        prediction = preds
      ))
    }

    # ----------------------------------------
    # ✅ profile / structured output
    # ----------------------------------------
    df_info <- scenario_info[
      scenario_info$scenario == name,
      , drop = FALSE
    ]

    df_info$prediction <- preds

    df_info
  })

  out
}
