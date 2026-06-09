#' Summarise DLNM effects
#'
#' Summarises DLNM effects on either the daily or accumulated scale,
#' using exposure grids derived from the observed data and the DLNM
#' specification stored in the fitted model.
#'
#' This function supports both deterministic summaries and uncertainty
#' propagation. When `uncertainty = TRUE`, summaries are computed from
#' simulated or posterior draws of the model coefficients. If
#' `output = "summary"`, the central estimate is computed as the median
#' of the simulated effects, while interval limits are obtained from
#' empirical quantiles (default: 2.5% and 97.5%).
#'
#' **Important:** when `uncertainty = TRUE` and `output = "summary"`,
#' columns such as `eta`, `effect`, `delta`, `delta_pp`, `baseline`,
#' and `predicted` represent the *central estimate*, computed as the
#' median of the simulated distribution.
#'
#' @param fit Fitted model.
#' @param wx_long Long-format weather data.
#' @param var Exposure variable name (character scalar), character vector of variables,
#'   or `NULL`. If `NULL`, all variables stored in fit metadata are used.
#' @param lag_max Optional maximum lag. Ignored if `fit` contains `epiexposure_spec`.
#' @param df_var Optional degrees of freedom for exposure. Ignored if `fit`
#'   contains `epiexposure_spec`.
#' @param df_lag Optional degrees of freedom for lag. Ignored if `fit`
#'   contains `epiexposure_spec`.
#' @param fun_var Optional exposure basis function. Ignored if `fit`
#'   contains `epiexposure_spec`.
#' @param fun_lag Optional lag basis function. Ignored if `fit`
#'   contains `epiexposure_spec`.
#' @param scale Character. `"daily"` or `"accumulated"`.
#' @param lag_periods Optional lag-period table used when
#'   `scale = "accumulated"` and `incremental = FALSE`. The table must
#'   contain columns `period`, `lag_start`, and `lag_end`.
#' @param probs Quantiles used to define the exposure grid.
#' @param ref Reference definition list.
#' @param effect_measure Character. `"percent"`, `"ratio"`, or `"linear"`.
#' @param incremental Logical. If `TRUE` and `scale = "accumulated"`,
#'   returns cumulative effects from lag 0 up to each lag.
#' @param uncertainty Logical. If `TRUE`, quantify uncertainty using
#'   simulated or posterior coefficient draws.
#' @param output Character. `"summary"` or `"samples"`.
#' @param n_samples Integer. Number of samples used for uncertainty.
#'
#' @return A data.frame.
#'
#' If `uncertainty = FALSE`, returns deterministic summaries of DLNM effects.
#'
#' If `uncertainty = TRUE` and `output = "summary"`, returns one row per
#' grid value / lag / period combination (depending on `scale`) with
#' median-based central estimates and empirical interval limits.
#'
#' If `uncertainty = TRUE` and `output = "samples"`, returns one row per
#' simulated sample.
#'
#' @details
#' Uncertainty is propagated using model-consistent sampling:
#' - Bayesian models (e.g., `brms`, `INLA`, `bdlnm`) use posterior draws
#' - Frequentist models use simulation from the asymptotic coefficient distribution
#'
#' For summary outputs under uncertainty, the median is used instead of the
#' mean to provide a more robust central estimate under asymmetric effect
#' distributions, which are common in nonlinear DLNM settings.
#'
#' @export
summarise_effects <- function(
    fit,
    wx_long,
    var = NULL,
    lag_max = NULL,
    df_var = NULL,
    df_lag = NULL,
    fun_var = NULL,
    fun_lag = NULL,
    scale = c("daily", "accumulated"),
    lag_periods = NULL,
    probs = seq(0.05, 0.95, by = 0.01),
    ref = list(method = "median", value = NULL),
    effect_measure = c("percent", "ratio", "linear"),
    incremental = FALSE,
    uncertainty = FALSE,
    output = c("summary", "samples"),
    n_samples = 1000
) {

  scale <- match.arg(scale)
  effect_measure <- match.arg(effect_measure)
  output <- match.arg(output)

  if (!is.data.frame(wx_long)) stop("`wx_long` must be a data.frame.")
  if (!all(c("epi_id", "dpp") %in% names(wx_long))) {
    stop("`wx_long` must contain at least 'epi_id' and 'dpp'.")
  }
  if (!is.logical(incremental) || length(incremental) != 1L) {
    stop("`incremental` must be TRUE or FALSE.")
  }
  if (!is.logical(uncertainty) || length(uncertainty) != 1L) {
    stop("`uncertainty` must be TRUE or FALSE.")
  }
  if (!is.numeric(n_samples) || length(n_samples) != 1L || n_samples <= 0) {
    stop("`n_samples` must be a positive integer.")
  }
  n_samples <- as.integer(n_samples)

  if (scale == "accumulated" && !incremental && is.null(lag_periods)) {
    stop("`lag_periods` must be provided when scale = 'accumulated' and incremental = FALSE.")
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  safe_quantile <- function(x, probs = c(0.025, 0.975)) {
    stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  }

  fit_spec <- attr(fit, "epiexposure_spec")
  fit_vars <- attr(fit, "epiexposure_vars")
  family_fit <- attr(fit, "epiexposure_family")
  cb_cols_fit <- attr(fit, "epiexposure_cb_cols")
  dat_template <- attr(fit, "epiexposure_dat_template")

  # ----------------------------------------------------------
  # Resolve variables to summarize
  # ----------------------------------------------------------
  if (is.null(var)) {
    vars <- fit_vars
    if (is.null(vars) || length(vars) == 0) {
      stop("`var` is NULL and no variables were found in fit metadata.")
    }
  } else {
    if (!is.character(var) || length(var) < 1) {
      stop("`var` must be NULL, a character scalar, or a character vector.")
    }
    vars <- var
  }

  missing_vars <- setdiff(vars, names(wx_long))
  if (length(missing_vars) > 0) {
    stop("The following variables are not present in `wx_long`: ",
         paste(missing_vars, collapse = ", "))
  }

  # ----------------------------------------------------------
  # Internal single-variable worker
  # ----------------------------------------------------------
  .summarise_effect_one_var <- function(var_one) {

    # ----------------------------------------------------------
    # 1) pooled series
    # ----------------------------------------------------------
    build_pooled_series <- function(dat, var, sep_n) {
      ids <- unique(dat$epi_id)
      out <- vector("list", length(ids))

      for (i in seq_along(ids)) {
        v <- dat |>
          dplyr::filter(epi_id == ids[i]) |>
          dplyr::arrange(dpp) |>
          dplyr::pull(.data[[var]])

        out[[i]] <- c(v, rep(NA_real_, sep_n))
      }

      unlist(out)
    }

    # ----------------------------------------------------------
    # 2) Resolve basis spec (prefer fit metadata)
    # ----------------------------------------------------------
    if (!is.null(fit_spec) && !is.null(fit_spec[[var_one]])) {

      spec_v <- fit_spec[[var_one]]
      lag_max_use <- as.integer(max(spec_v$lag_max))
      argvar <- spec_v$argvar
      arglag <- spec_v$arglag

    } else {

      if (is.null(lag_max)) {
        stop("Model does not contain `epiexposure_spec` for variable '", var_one,
             "'. Please provide `lag_max`.")
      }

      lag_max_use <- as.integer(lag_max)

      fun_var_use <- fun_var %||% "ns"
      fun_lag_use <- fun_lag %||% "ns"
      df_var_use  <- df_var %||% 4
      df_lag_use  <- df_lag %||% 4

      argvar <- switch(
        fun_var_use,
        ns   = list(fun = "ns", df = df_var_use),
        bs   = list(fun = "bs", df = df_var_use),
        poly = list(fun = "poly", degree = df_var_use),
        lin  = list(fun = "lin"),
        stop("Unsupported fun_var: ", fun_var_use)
      )

      if (!is.null(argvar$fun) && argvar$fun != "lin") {
        argvar$intercept <- FALSE
      }

      arglag <- switch(
        fun_lag_use,
        ns  = list(fun = "ns", df = df_lag_use),
        ps  = list(fun = "ps", df = df_lag_use),
        lin = list(fun = "lin"),
        stop("Unsupported fun_lag: ", fun_lag_use)
      )
    }

    x_pool <- build_pooled_series(wx_long, var_one, lag_max_use)

    cb <- dlnm::crossbasis(
      x_pool,
      lag    = lag_max_use,
      argvar = argvar,
      arglag = arglag
    )

    cb_cols_var <- NULL
    if (!is.null(cb_cols_fit)) {
      cb_cols_var <- grep(paste0("^cb_", var_one, "_"), cb_cols_fit, value = TRUE)
    } else if (!is.null(dat_template)) {
      cb_cols_var <- grep(paste0("^cb_", var_one, "_"), names(dat_template), value = TRUE)
    }

    # ----------------------------------------------------------
    # 3) Extract beta / vcov / draws
    # ----------------------------------------------------------
    extract_coef_vcov <- function(model) {

      if (inherits(model, "glmmTMB")) {
        return(list(
          beta = glmmTMB::fixef(model)$cond,
          vcov = as.matrix(stats::vcov(model)$cond)
        ))
      }

      if (inherits(model, "merMod")) {
        return(list(
          beta = lme4::fixef(model),
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "lme")) {
        return(list(
          beta = nlme::fixef(model),
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "gls")) {
        return(list(
          beta = stats::coef(model),
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "gam")) {
        return(list(
          beta = stats::coef(model),
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "HLfit")) {
        return(list(
          beta = spaMM::fixef(model),
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "brmsfit")) {
        fe <- brms::fixef(model)
        beta <- fe[, "Estimate"]
        names(beta) <- rownames(fe)
        return(list(
          beta = beta,
          vcov = as.matrix(stats::vcov(model))
        ))
      }

      if (inherits(model, "inla")) {
        beta <- model$summary.fixed$mean
        V <- diag(model$summary.fixed$sd^2)
        rownames(V) <- names(beta)
        colnames(V) <- names(beta)
        return(list(beta = beta, vcov = V))
      }

      if (inherits(model, "bdlnm")) {
        beta <- model$coefficients.summary[, "mean"]
        V <- stats::cov(t(model$coefficients))
        return(list(beta = beta, vcov = V))
      }

      return(list(
        beta = stats::coef(model),
        vcov = as.matrix(stats::vcov(model))
      ))
    }

    extract_beta_draws <- function(model, n_samples, var, cb_ncol, cb_cols_var = NULL) {

      # ---------- brms: REAL posterior draws ----------
      if (inherits(model, "brmsfit")) {
        draws <- as.matrix(brms::as_draws_matrix(model))

        nm <- if (!is.null(cb_cols_var)) {
          paste0("b_", cb_cols_var)
        } else {
          grep(paste0("^b_cb_", var, "_"), colnames(draws), value = TRUE)
        }

        if (length(nm) == 0) {
          nm <- grep(paste0("^b_cb_", var, "_"), colnames(draws), value = TRUE)
        }
        if (length(nm) == 0) {
          nm <- grep(paste0("^b_.*", var), colnames(draws), value = TRUE)
        }

        if (length(nm) != cb_ncol) {
          stop("Could not match brms posterior draws to crossbasis columns.")
        }

        if (nrow(draws) > n_samples) {
          set.seed(1)
          keep <- sample(seq_len(nrow(draws)), n_samples)
          draws <- draws[keep, nm, drop = FALSE]
        } else {
          draws <- draws[, nm, drop = FALSE]
        }

        return(draws)
      }

      # ---------- INLA: REAL posterior draws ----------
      if (inherits(model, "inla")) {
        if (!requireNamespace("INLA", quietly = TRUE)) {
          stop("Package 'INLA' is required for INLA uncertainty.")
        }

        posterior <- INLA::inla.posterior.sample(n = n_samples, result = model)

        cb_names <- cb_cols_var %||% paste0("cb_", var, "_", seq_len(cb_ncol))

        draws <- do.call(rbind, lapply(posterior, function(s) {
          latent <- s$latent
          names(latent) <- gsub(":1$", "", names(latent))
          as.numeric(latent[cb_names])
        }))

        colnames(draws) <- cb_names
        return(draws)
      }

      # ---------- bdlnm: REAL posterior draws ----------
      if (inherits(model, "bdlnm")) {
        beta_draws <- model$coefficients
        if (is.null(dim(beta_draws))) {
          beta_draws <- matrix(beta_draws, ncol = 1)
        }

        cb_names <- intersect(rownames(beta_draws), cb_cols_var %||% character(0))
        if (length(cb_names) == 0) {
          cb_names <- intersect(rownames(beta_draws), paste0("cb_", var, "_", seq_len(cb_ncol)))
        }
        if (length(cb_names) == 0) {
          cb_names <- grep(paste0("^cb_", var, "_"), rownames(beta_draws), value = TRUE)
        }

        if (length(cb_names) != cb_ncol) {
          stop("Could not match bdlnm posterior draws to crossbasis columns.")
        }

        if (ncol(beta_draws) > n_samples) {
          set.seed(1)
          keep <- sample(seq_len(ncol(beta_draws)), n_samples)
          beta_draws <- beta_draws[cb_names, keep, drop = FALSE]
        } else {
          beta_draws <- beta_draws[cb_names, , drop = FALSE]
        }

        return(t(beta_draws))
      }

      # ---------- Frequentist: normal approximation ----------
      cv <- extract_coef_vcov(model)
      beta_hat <- cv$beta
      V_hat <- cv$vcov

      cb_names <- if (!is.null(cb_cols_var)) {
        intersect(cb_cols_var, names(beta_hat))
      } else {
        grep(paste0("^cb_", var, "_"), names(beta_hat), value = TRUE)
      }

      if (length(cb_names) != cb_ncol) {
        stop("Could not match coefficient names to crossbasis structure.")
      }

      beta_hat <- beta_hat[cb_names]
      V_hat <- V_hat[cb_names, cb_names, drop = FALSE]

      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for frequentist uncertainty approximation.")
      }

      draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)
      if (is.null(dim(draws))) {
        draws <- matrix(draws, nrow = 1)
      }
      colnames(draws) <- cb_names
      draws
    }

    bv <- extract_coef_vcov(fit)
    cf <- bv$beta
    vc <- bv$vcov

    idx <- grepl(paste0("^cb_", var_one, "_"), names(cf))
    beta <- cf[idx]
    vc_sub <- vc[idx, idx, drop = FALSE]

    if (length(beta) != ncol(cb)) stop("length(beta) != ncol(cb) for variable: ", var_one)
    if (!(nrow(vc_sub) == ncol(cb) && ncol(vc_sub) == ncol(cb))) {
      stop("Variance-covariance matrix does not match crossbasis dimensions for variable: ", var_one)
    }

    # ----------------------------------------------------------
    # 4) exposure grid
    # ----------------------------------------------------------
    x_all <- wx_long[[var_one]]

    at_vals <- sort(unique(
      as.numeric(stats::quantile(x_all, probs = probs, na.rm = TRUE))
    ))

    cen <- switch(
      ref$method,
      median     = stats::median(x_all, na.rm = TRUE),
      percentile = stats::quantile(x_all, ref$value, na.rm = TRUE),
      fixed      = ref$value,
      stop("Invalid ref$method")
    )
    cen <- as.numeric(cen)

    # ----------------------------------------------------------
    # 5) deterministic crosspred
    # ----------------------------------------------------------
    cp <- dlnm::crosspred(
      cb,
      coef  = beta,
      vcov  = vc_sub,
      at    = at_vals,
      cen   = cen,
      bylag = 1
    )

    mf <- cp$matfit

    lag_idx <- suppressWarnings(as.integer(gsub("lag", "", colnames(mf))))
    if (anyNA(lag_idx)) lag_idx <- 0:(ncol(mf) - 1)

    ord <- order(lag_idx)
    lag_idx <- lag_idx[ord]
    mf <- mf[, ord, drop = FALSE]

    # ----------------------------------------------------------
    # 6) Baseline severity for DELTA
    # ----------------------------------------------------------
    get_linkfun <- function(family_fit) {
      if (is.character(family_fit)) {
        if (family_fit %in% c("beta", "binomial")) return(qlogis)
        if (family_fit %in% c("poisson", "gamma")) return(log)
        if (family_fit %in% c("gaussian")) return(identity)
      }
      if (inherits(family_fit, "family") && !is.null(family_fit$link)) {
        if (family_fit$link == "logit") return(qlogis)
        if (family_fit$link == "log") return(log)
        if (family_fit$link == "identity") return(identity)
      }
      qlogis
    }

    get_linkinv <- function(family_fit) {
      if (is.character(family_fit)) {
        if (family_fit %in% c("beta", "binomial")) return(plogis)
        if (family_fit %in% c("poisson", "gamma")) return(exp)
        if (family_fit %in% c("gaussian")) return(identity)
      }
      if (inherits(family_fit, "family") && !is.null(family_fit$linkinv)) {
        return(family_fit$linkinv)
      }
      plogis
    }

    linkfun <- get_linkfun(family_fit)
    linkinv <- get_linkinv(family_fit)

    baseline_response <- NA_real_

    if (!is.null(fit_vars) && all(fit_vars %in% names(wx_long))) {

      ref_profiles <- lapply(fit_vars, function(v) {
        spec_v2 <- fit_spec[[v]]
        lag_v <- if (!is.null(spec_v2$lag_max)) max(spec_v2$lag_max) else lag_max_use
        ref_v <- as.numeric(stats::median(wx_long[[v]], na.rm = TRUE))
        rep(ref_v, lag_v + 1L)
      })
      names(ref_profiles) <- fit_vars

      baseline_response <- tryCatch(
        {
          as.numeric(predict_outcome(
            fit = fit,
            profiles = ref_profiles,
            re = "population",
            type = "response",
            uncertainty = FALSE
          )$prediction[1])
        },
        error = function(e) NA_real_
      )
    }

    transform_effect <- function(eta) {
      if (effect_measure == "linear") return(eta)
      if (effect_measure == "ratio")  return(exp(eta))
      (exp(eta) - 1) * 100
    }

    add_delta <- function(eta_vec) {
      if (is.na(baseline_response)) {
        return(list(
          delta = rep(NA_real_, length(eta_vec)),
          delta_pp = rep(NA_real_, length(eta_vec)),
          baseline = rep(NA_real_, length(eta_vec)),
          predicted = rep(NA_real_, length(eta_vec))
        ))
      }

      pred_resp <- linkinv(linkfun(baseline_response) + eta_vec)
      delta <- pred_resp - baseline_response
      delta_pp <- 100 * delta

      list(
        delta = delta,
        delta_pp = delta_pp,
        baseline = rep(baseline_response, length(eta_vec)),
        predicted = pred_resp
      )
    }

    build_daily_df <- function(eta_mat) {

      df <- as.data.frame(eta_mat)
      colnames(df) <- paste0("lag_", lag_idx)

      out <- df |>
        dplyr::mutate(value = at_vals) |>
        tidyr::pivot_longer(
          cols = tidyselect::starts_with("lag_"),
          names_to = "lag",
          values_to = "eta"
        ) |>
        dplyr::mutate(
          lag = as.integer(gsub("lag_", "", lag)),
          effect = transform_effect(eta),
          scale = "daily",
          var = var_one
        )

      dd <- add_delta(out$eta)
      out$delta <- dd$delta
      out$delta_pp <- dd$delta_pp
      out$baseline <- dd$baseline
      out$predicted <- dd$predicted

      out
    }

    build_accumulated_df <- function(eta_mat) {

      if (incremental) {
        eta_cum <- t(apply(eta_mat, 1, cumsum))

        df <- as.data.frame(eta_cum)
        colnames(df) <- paste0("lag_", lag_idx)

        out <- df |>
          dplyr::mutate(value = at_vals) |>
          tidyr::pivot_longer(
            cols = tidyselect::starts_with("lag_"),
            names_to = "lag",
            values_to = "eta"
          ) |>
          dplyr::mutate(
            lag = as.integer(gsub("lag_", "", lag)),
            period = paste0("0-", lag),
            effect = transform_effect(eta),
            scale = "accumulated",
            var = var_one
          )

        dd <- add_delta(out$eta)
        out$delta <- dd$delta
        out$delta_pp <- dd$delta_pp
        out$baseline <- dd$baseline
        out$predicted <- dd$predicted

        return(out)
      }

      if (!all(c("period", "lag_start", "lag_end") %in% names(lag_periods))) {
        stop("`lag_periods` must contain columns: period, lag_start, lag_end.")
      }

      purrr::map_dfr(seq_len(nrow(lag_periods)), function(i) {

        sel <- which(
          lag_idx >= lag_periods$lag_start[i] &
            lag_idx <= lag_periods$lag_end[i]
        )

        eta_sum <- rowSums(eta_mat[, sel, drop = FALSE])

        out <- data.frame(
          value  = at_vals,
          period = lag_periods$period[i],
          eta    = eta_sum,
          effect = transform_effect(eta_sum),
          scale  = "accumulated",
          var    = var_one
        )

        dd <- add_delta(out$eta)
        out$delta <- dd$delta
        out$delta_pp <- dd$delta_pp
        out$baseline <- dd$baseline
        out$predicted <- dd$predicted

        out
      })
    }

    # ----------------------------------------------------------
    # 7) no uncertainty
    # ----------------------------------------------------------
    if (!uncertainty) {
      if (scale == "daily") {
        return(build_daily_df(mf))
      } else {
        return(build_accumulated_df(mf))
      }
    }

    # ----------------------------------------------------------
    # 8) uncertainty
    # ----------------------------------------------------------
    beta_draws <- extract_beta_draws(
      model = fit,
      n_samples = n_samples,
      var = var_one,
      cb_ncol = ncol(cb),
      cb_cols_var = cb_cols_var
    )

    if (output == "samples") {

      out_list <- vector("list", nrow(beta_draws))

      for (i in seq_len(nrow(beta_draws))) {

        cp_i <- dlnm::crosspred(
          cb,
          coef  = beta_draws[i, ],
          vcov  = diag(0, length(beta_draws[i, ])),
          at    = at_vals,
          cen   = cen,
          bylag = 1
        )

        mf_i <- cp_i$matfit
        mf_i <- mf_i[, ord, drop = FALSE]

        tmp <- if (scale == "daily") {
          build_daily_df(mf_i)
        } else {
          build_accumulated_df(mf_i)
        }

        tmp$sample <- i
        out_list[[i]] <- tmp
      }

      out <- do.call(rbind, out_list)
      nm <- names(out)
      out <- out[, c("sample", nm[nm != "sample"]), drop = FALSE]

      return(out)
    }

    # summary
    draw_list <- vector("list", nrow(beta_draws))

    for (i in seq_len(nrow(beta_draws))) {

      cp_i <- dlnm::crosspred(
        cb,
        coef  = beta_draws[i, ],
        vcov  = NULL,
        at    = at_vals,
        cen   = cen,
        bylag = 1
      )

      mf_i <- cp_i$matfit
      mf_i <- mf_i[, ord, drop = FALSE]

      draw_list[[i]] <- if (scale == "daily") {
        build_daily_df(mf_i)
      } else {
        build_accumulated_df(mf_i)
      }
    }

    all_draws <- dplyr::bind_rows(draw_list, .id = "sample")
    all_draws$sample <- as.integer(all_draws$sample)

    group_vars <- intersect(c("var", "value", "lag", "period", "scale"), names(all_draws))

    summarised <- all_draws |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::summarise(
        eta = stats::median(.data$eta, na.rm = TRUE),
        eta_sd = stats::sd(.data$eta, na.rm = TRUE),
        eta_lower = safe_quantile(.data$eta)[1],
        eta_upper = safe_quantile(.data$eta)[2],

        effect = stats::median(.data$effect, na.rm = TRUE),
        effect_sd = stats::sd(.data$effect, na.rm = TRUE),
        effect_lower = safe_quantile(.data$effect)[1],
        effect_upper = safe_quantile(.data$effect)[2],

        delta = stats::median(.data$delta, na.rm = TRUE),
        delta_sd = stats::sd(.data$delta, na.rm = TRUE),
        delta_lower = safe_quantile(.data$delta)[1],
        delta_upper = safe_quantile(.data$delta)[2],

        delta_pp = stats::median(.data$delta_pp, na.rm = TRUE),
        delta_pp_sd = stats::sd(.data$delta_pp, na.rm = TRUE),
        delta_pp_lower = safe_quantile(.data$delta_pp)[1],
        delta_pp_upper = safe_quantile(.data$delta_pp)[2],

        baseline = stats::median(.data$baseline, na.rm = TRUE),
        predicted = stats::median(.data$predicted, na.rm = TRUE),
        predicted_sd = stats::sd(.data$predicted, na.rm = TRUE),
        predicted_lower = safe_quantile(.data$predicted)[1],
        predicted_upper = safe_quantile(.data$predicted)[2],
        .groups = "drop"
      )

    as.data.frame(summarised)
  }

  # ----------------------------------------------------------
  # Apply over one or many variables
  # ----------------------------------------------------------
  out_list <- lapply(vars, .summarise_effect_one_var)
  out <- do.call(rbind, out_list)

  rownames(out) <- NULL
  out
}
