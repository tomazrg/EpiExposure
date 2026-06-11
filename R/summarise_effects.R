#' Summarise DLNM effects
#'
#' Summarises DLNM effects on either the daily or accumulated scale,
#' using exposure grids derived from the observed data and the DLNM
#' specification stored in the fitted model.
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

  safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1) return(0)
    stats::sd(x)
  }

  # ----------------------------------------------------------
  # ✅ helper: matching robusto de nomes para draws do brms
  # ----------------------------------------------------------
  match_brms_draw_names <- function(cb_names_ref, draw_colnames) {

    bnames <- paste0("b_", cb_names_ref)

    if (all(bnames %in% draw_colnames)) {
      return(bnames)
    }

    # fallback defensivo
    raw_match <- cb_names_ref[cb_names_ref %in% draw_colnames]
    if (length(raw_match) == length(cb_names_ref)) {
      return(raw_match)
    }

    stop("Could not match brms posterior draw names to crossbasis columns.")
  }

  fit_spec     <- attr(fit, "epiexposure_spec")
  fit_vars     <- attr(fit, "epiexposure_vars")
  family_fit   <- attr(fit, "epiexposure_family")
  cb_cols_fit  <- attr(fit, "epiexposure_cb_cols")
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

  # ==========================================================
  # ✅ AJUSTE 4 — helper para ordenar coeficientes cb_* por índice
  # ==========================================================
  sort_cb_names <- function(x) {
    if (length(x) == 0) return(x)
    idx <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))
    idx[is.na(idx)] <- seq_along(x)
    x[order(idx)]
  }

  # ----------------------------------------------------------
  # Worker por variável
  # ----------------------------------------------------------
  .summarise_effect_one_var <- function(var_one) {

    # ----------------------------------------------------------
    # pooled series
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
    # Resolve basis spec
    # ----------------------------------------------------------
    if (!is.null(fit_spec) && !is.null(fit_spec[[var_one]])) {

      spec_v <- fit_spec[[var_one]]
      # ======================================================
      # ✅ AJUSTE 1 — DIMENSÃO TEMPORAL
      # ======================================================
      lag_max_use <- as.integer(max(spec_v$lag_max))
      argvar <- spec_v$argvar
      arglag <- spec_v$arglag

    } else {

      if (is.null(lag_max)) {
        stop("Model does not contain `epiexposure_spec` for variable '", var_one,
             "'. Please provide `lag_max`.")
      }

      # ======================================================
      # ✅ AJUSTE 1 — DIMENSÃO TEMPORAL (fallback)
      # ======================================================
      lag_max_use <- as.integer(max(lag_max))

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

    # ======================================================
    # ✅ AJUSTE EXTRA — CHECK AUTOMÁTICO DE COBERTURA TEMPORAL
    # ======================================================
    .check_lag_coverage <- function(dat, lag_max) {
      n_required <- lag_max + 1

      bad_ids <- dat |>
        dplyr::group_by(epi_id) |>
        dplyr::summarise(n_days = dplyr::n_distinct(dpp), .groups = "drop") |>
        dplyr::filter(n_days < n_required)

      if (nrow(bad_ids) > 0) {
        stop(
          paste0(
            "Some epidemics do not have enough temporal coverage for lag_max.\n",
            "Required days per epi_id: ", n_required, "\n",
            "Example problematic epi_id: ",
            paste(head(bad_ids$epi_id, 5), collapse = ", ")
          )
        )
      }
    }

    .check_lag_coverage(wx_long, lag_max_use)

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
      cb_cols_var <- sort_cb_names(cb_cols_var)
    } else if (!is.null(dat_template)) {
      cb_cols_var <- grep(paste0("^cb_", var_one, "_"), names(dat_template), value = TRUE)
      cb_cols_var <- sort_cb_names(cb_cols_var)
    }

    # ----------------------------------------------------------
    # Extract coef / vcov
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
        V <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
        if (is.null(V)) stop("Could not extract vcov from spaMM model.")
        return(list(
          beta = spaMM::fixef(model),
          vcov = V
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

    bv <- extract_coef_vcov(fit)
    cf <- bv$beta
    vc <- bv$vcov

    # ==========================================================
    # ✅ AJUSTE 4 — ORDEM DOS COEFICIENTES: REFERÊNCIA ÚNICA
    # ==========================================================
    cb_names_ref <- if (!is.null(cb_cols_var) && length(cb_cols_var) > 0) {
      cb_cols_var[cb_cols_var %in% names(cf)]
    } else {
      sort_cb_names(grep(paste0("^cb_", var_one, "_"), names(cf), value = TRUE))
    }

    if (length(cb_names_ref) != ncol(cb)) {
      stop(
        "Could not match coefficient names to crossbasis columns for variable '", var_one,
        "'. Expected ", ncol(cb), " and found ", length(cb_names_ref), "."
      )
    }

    beta <- cf[cb_names_ref]
    vc_sub <- vc[cb_names_ref, cb_names_ref, drop = FALSE]

    # ----------------------------------------------------------
    # draws por engine
    # ----------------------------------------------------------
    extract_beta_draws <- function(model, n_samples, cb_names_ref) {

      # ---------- brms ----------
      if (inherits(model, "brmsfit")) {
        draws <- posterior::as_draws_df(model)
        bnames <- match_brms_draw_names(cb_names_ref, names(draws))

        draws <- as.matrix(draws[, bnames, drop = FALSE])

        if (nrow(draws) > n_samples) {
          set.seed(1)
          keep <- sample(seq_len(nrow(draws)), n_samples)
          draws <- draws[keep, , drop = FALSE]
        }

        colnames(draws) <- cb_names_ref
        return(draws)
      }

      # ---------- INLA ----------
      if (inherits(model, "inla")) {
        if (!requireNamespace("INLA", quietly = TRUE)) {
          stop("Package 'INLA' is required for INLA uncertainty.")
        }

        posterior <- INLA::inla.posterior.sample(n = n_samples, result = model)

        draws <- do.call(rbind, lapply(posterior, function(s) {
          latent <- s$latent
          names(latent) <- gsub(":1$", "", names(latent))
          as.numeric(latent[cb_names_ref])
        }))

        colnames(draws) <- cb_names_ref
        return(draws)
      }

      # ---------- bdlnm ----------
      if (inherits(model, "bdlnm")) {
        beta_draws <- model$coefficients
        if (is.null(dim(beta_draws))) {
          beta_draws <- matrix(beta_draws, ncol = 1)
        }

        if (!all(cb_names_ref %in% rownames(beta_draws))) {
          stop("Could not match bdlnm posterior draws to crossbasis columns.")
        }

        beta_draws <- beta_draws[cb_names_ref, , drop = FALSE]

        if (ncol(beta_draws) > n_samples) {
          set.seed(1)
          keep <- sample(seq_len(ncol(beta_draws)), n_samples)
          beta_draws <- beta_draws[, keep, drop = FALSE]
        }

        out <- t(beta_draws)
        colnames(out) <- cb_names_ref
        return(out)
      }

      # ---------- frequentista ----------
      beta_hat <- cf[cb_names_ref]
      V_hat <- vc[cb_names_ref, cb_names_ref, drop = FALSE]

      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop("Package 'MASS' is required for frequentist uncertainty approximation.")
      }

      draws <- MASS::mvrnorm(n = n_samples, mu = beta_hat, Sigma = V_hat)

      if (is.null(dim(draws))) {
        draws <- matrix(draws, nrow = 1)
      }

      colnames(draws) <- cb_names_ref
      draws
    }

    # ----------------------------------------------------------
    # exposure grid
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
    # deterministic crosspred
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
    # baseline
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
        lag_v <- if (!is.null(spec_v2$lag_max)) as.integer(max(spec_v2$lag_max)) else lag_max_use
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
        error = function(e) {
          warning("Baseline prediction failed in summarise_effects(): ", conditionMessage(e))
          NA_real_
        }
      )
    }

    transform_effect <- function(eta_mat) {
      if (effect_measure == "linear") return(eta_mat)
      if (effect_measure == "ratio")  return(exp(eta_mat))
      (exp(eta_mat) - 1) * 100
    }

    # ----------------------------------------------------------
    # deterministic output helper
    # ----------------------------------------------------------
    build_daily_df <- function(eta_mat) {
      grid <- expand.grid(
        value = at_vals,
        lag   = lag_idx,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      df <- data.frame(
        value = grid$value,
        lag = grid$lag,
        eta = as.vector(eta_mat),
        stringsAsFactors = FALSE
      )

      df$effect <- as.vector(transform_effect(eta_mat))
      df$scale <- "daily"
      df$var <- var_one

      if (is.na(baseline_response)) {
        df$delta <- NA_real_
        df$delta_pp <- NA_real_
        df$baseline <- NA_real_
        df$predicted <- NA_real_
      } else {
        pred_mat <- linkinv(linkfun(baseline_response) + eta_mat)
        delta_mat <- pred_mat - baseline_response
        df$delta <- as.vector(delta_mat)
        df$delta_pp <- as.vector(100 * delta_mat)
        df$baseline <- baseline_response
        df$predicted <- as.vector(pred_mat)
      }

      df
    }

    build_accumulated_df <- function(eta_mat) {

      if (incremental) {
        eta_cum <- t(apply(eta_mat, 1, cumsum))

        grid <- expand.grid(
          value = at_vals,
          lag   = lag_idx,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )

        df <- data.frame(
          value = grid$value,
          lag = grid$lag,
          period = paste0("0-", grid$lag),
          eta = as.vector(eta_cum),
          stringsAsFactors = FALSE
        )

        df$effect <- as.vector(transform_effect(eta_cum))
        df$scale <- "accumulated"
        df$var <- var_one

        if (is.na(baseline_response)) {
          df$delta <- NA_real_
          df$delta_pp <- NA_real_
          df$baseline <- NA_real_
          df$predicted <- NA_real_
        } else {
          pred_mat <- linkinv(linkfun(baseline_response) + eta_cum)
          delta_mat <- pred_mat - baseline_response
          df$delta <- as.vector(delta_mat)
          df$delta_pp <- as.vector(100 * delta_mat)
          df$baseline <- baseline_response
          df$predicted <- as.vector(pred_mat)
        }

        return(df)
      }

      if (!all(c("period", "lag_start", "lag_end") %in% names(lag_periods))) {
        stop("`lag_periods` must contain columns: period, lag_start, lag_end.")
      }

      out <- purrr::map_dfr(seq_len(nrow(lag_periods)), function(i) {
        sel <- which(lag_idx >= lag_periods$lag_start[i] &
                       lag_idx <= lag_periods$lag_end[i])

        eta_sum <- rowSums(eta_mat[, sel, drop = FALSE])
        eff_sum <- transform_effect(eta_sum)

        tmp <- data.frame(
          value = at_vals,
          period = lag_periods$period[i],
          eta = eta_sum,
          effect = eff_sum,
          scale = "accumulated",
          var = var_one,
          stringsAsFactors = FALSE
        )

        if (is.na(baseline_response)) {
          tmp$delta <- NA_real_
          tmp$delta_pp <- NA_real_
          tmp$baseline <- NA_real_
          tmp$predicted <- NA_real_
        } else {
          pred <- linkinv(linkfun(baseline_response) + eta_sum)
          delta <- pred - baseline_response
          tmp$delta <- delta
          tmp$delta_pp <- 100 * delta
          tmp$baseline <- baseline_response
          tmp$predicted <- pred
        }

        tmp
      })

      out
    }

    # ----------------------------------------------------------
    # sem incerteza
    # ----------------------------------------------------------
    if (!uncertainty) {
      if (scale == "daily") {
        return(build_daily_df(mf))
      } else {
        return(build_accumulated_df(mf))
      }
    }

    # ----------------------------------------------------------
    # com incerteza: estratégia por array
    # ----------------------------------------------------------
    beta_draws <- extract_beta_draws(
      model = fit,
      n_samples = n_samples,
      cb_names_ref = cb_names_ref
    )

    # ==========================================================
    # ✅ AJUSTE 3 — DIAGNÓSTICO DE VARIABILIDADE DOS DRAWS
    # ==========================================================
    draw_sd <- apply(beta_draws, 2, stats::sd)
    if (all(!is.finite(draw_sd)) || all(draw_sd < 1e-12, na.rm = TRUE)) {
      warning(
        "Near-zero variability detected in coefficient draws for variable '", var_one,
        "'. Summary intervals may collapse to a single value."
      )
    }

    # gerar uma matriz de eta por draw
    eta_list <- vector("list", nrow(beta_draws))

    for (i in seq_len(nrow(beta_draws))) {
      beta_i <- beta_draws[i, ]

      cp_i <- dlnm::crosspred(
        cb,
        coef  = beta_i,
        # ======================================================
        # ✅ AJUSTE 2 — API do crosspred()
        # ======================================================
        vcov  = diag(0, length(beta_i)),
        at    = at_vals,
        cen   = cen,
        bylag = 1
      )

      mf_i <- cp_i$matfit
      mf_i <- mf_i[, ord, drop = FALSE]
      eta_list[[i]] <- mf_i
    }

    if (output == "samples") {
      out_list <- vector("list", length(eta_list))

      for (i in seq_along(eta_list)) {
        tmp <- if (scale == "daily") {
          build_daily_df(eta_list[[i]])
        } else {
          build_accumulated_df(eta_list[[i]])
        }
        tmp$sample <- i
        out_list[[i]] <- tmp
      }

      out <- do.call(rbind, out_list)
      nm <- names(out)
      out <- out[, c("sample", nm[nm != "sample"]), drop = FALSE]
      return(out)
    }

    # ----------------------------------------------------------
    # SUMMARY por array
    # ----------------------------------------------------------
    if (scale == "daily") {

      n_at <- length(at_vals)
      n_lag <- length(lag_idx)
      n_draw <- length(eta_list)

      eta_arr <- array(NA_real_, dim = c(n_at, n_lag, n_draw))
      for (i in seq_len(n_draw)) {
        eta_arr[, , i] <- eta_list[[i]]
      }

      summarize_arr <- function(arr) {
        list(
          est   = apply(arr, c(1, 2), stats::median, na.rm = TRUE),
          sd    = apply(arr, c(1, 2), safe_sd),
          lower = apply(arr, c(1, 2), function(z) safe_quantile(z)[1]),
          upper = apply(arr, c(1, 2), function(z) safe_quantile(z)[2])
        )
      }

      eta_sum <- summarize_arr(eta_arr)
      eff_arr <- transform_effect(eta_arr)
      eff_sum <- summarize_arr(eff_arr)

      if (is.na(baseline_response)) {
        delta_sum <- list(
          est   = matrix(NA_real_, n_at, n_lag),
          sd    = matrix(NA_real_, n_at, n_lag),
          lower = matrix(NA_real_, n_at, n_lag),
          upper = matrix(NA_real_, n_at, n_lag)
        )
        delta_pp_sum <- delta_sum
        pred_sum <- delta_sum
        baseline_mat <- matrix(NA_real_, n_at, n_lag)
      } else {
        pred_arr <- linkinv(linkfun(baseline_response) + eta_arr)
        delta_arr <- pred_arr - baseline_response
        delta_pp_arr <- 100 * delta_arr

        pred_sum <- summarize_arr(pred_arr)
        delta_sum <- summarize_arr(delta_arr)
        delta_pp_sum <- summarize_arr(delta_pp_arr)
        baseline_mat <- matrix(baseline_response, nrow = n_at, ncol = n_lag)
      }

      grid <- expand.grid(
        value = at_vals,
        lag   = lag_idx,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      out <- data.frame(
        var = var_one,
        lag = grid$lag,
        scale = "daily",
        value = grid$value,

        eta = as.vector(eta_sum$est),
        eta_sd = as.vector(eta_sum$sd),
        eta_lower = as.vector(eta_sum$lower),
        eta_upper = as.vector(eta_sum$upper),

        effect = as.vector(eff_sum$est),
        effect_sd = as.vector(eff_sum$sd),
        effect_lower = as.vector(eff_sum$lower),
        effect_upper = as.vector(eff_sum$upper),

        delta = as.vector(delta_sum$est),
        delta_sd = as.vector(delta_sum$sd),
        delta_lower = as.vector(delta_sum$lower),
        delta_upper = as.vector(delta_sum$upper),

        delta_pp = as.vector(delta_pp_sum$est),
        delta_pp_sd = as.vector(delta_pp_sum$sd),
        delta_pp_lower = as.vector(delta_pp_sum$lower),
        delta_pp_upper = as.vector(delta_pp_sum$upper),

        baseline = as.vector(baseline_mat),
        predicted = as.vector(pred_sum$est),
        predicted_sd = as.vector(pred_sum$sd),
        predicted_lower = as.vector(pred_sum$lower),
        predicted_upper = as.vector(pred_sum$upper),

        stringsAsFactors = FALSE
      )

      return(out)
    }

    # ----------------------------------------------------------
    # SUMMARY acumulado incremental
    # ----------------------------------------------------------
    if (incremental) {

      n_at <- length(at_vals)
      n_lag <- length(lag_idx)
      n_draw <- length(eta_list)

      eta_arr <- array(NA_real_, dim = c(n_at, n_lag, n_draw))
      for (i in seq_len(n_draw)) {
        eta_arr[, , i] <- t(apply(eta_list[[i]], 1, cumsum))
      }

      summarize_arr <- function(arr) {
        list(
          est   = apply(arr, c(1, 2), stats::median, na.rm = TRUE),
          sd    = apply(arr, c(1, 2), safe_sd),
          lower = apply(arr, c(1, 2), function(z) safe_quantile(z)[1]),
          upper = apply(arr, c(1, 2), function(z) safe_quantile(z)[2])
        )
      }

      eta_sum <- summarize_arr(eta_arr)
      eff_arr <- transform_effect(eta_arr)
      eff_sum <- summarize_arr(eff_arr)

      if (is.na(baseline_response)) {
        delta_sum <- list(
          est   = matrix(NA_real_, n_at, n_lag),
          sd    = matrix(NA_real_, n_at, n_lag),
          lower = matrix(NA_real_, n_at, n_lag),
          upper = matrix(NA_real_, n_at, n_lag)
        )
        delta_pp_sum <- delta_sum
        pred_sum <- delta_sum
        baseline_mat <- matrix(NA_real_, n_at, n_lag)
      } else {
        pred_arr <- linkinv(linkfun(baseline_response) + eta_arr)
        delta_arr <- pred_arr - baseline_response
        delta_pp_arr <- 100 * delta_arr

        pred_sum <- summarize_arr(pred_arr)
        delta_sum <- summarize_arr(delta_arr)
        delta_pp_sum <- summarize_arr(delta_pp_arr)
        baseline_mat <- matrix(baseline_response, nrow = n_at, ncol = n_lag)
      }

      grid <- expand.grid(
        value = at_vals,
        lag   = lag_idx,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      out <- data.frame(
        var = var_one,
        lag = grid$lag,
        period = paste0("0-", grid$lag),
        scale = "accumulated",
        value = grid$value,

        eta = as.vector(eta_sum$est),
        eta_sd = as.vector(eta_sum$sd),
        eta_lower = as.vector(eta_sum$lower),
        eta_upper = as.vector(eta_sum$upper),

        effect = as.vector(eff_sum$est),
        effect_sd = as.vector(eff_sum$sd),
        effect_lower = as.vector(eff_sum$lower),
        effect_upper = as.vector(eff_sum$upper),

        delta = as.vector(delta_sum$est),
        delta_sd = as.vector(delta_sum$sd),
        delta_lower = as.vector(delta_sum$lower),
        delta_upper = as.vector(delta_sum$upper),

        delta_pp = as.vector(delta_pp_sum$est),
        delta_pp_sd = as.vector(delta_pp_sum$sd),
        delta_pp_lower = as.vector(delta_pp_sum$lower),
        delta_pp_upper = as.vector(delta_pp_sum$upper),

        baseline = as.vector(baseline_mat),
        predicted = as.vector(pred_sum$est),
        predicted_sd = as.vector(pred_sum$sd),
        predicted_lower = as.vector(pred_sum$lower),
        predicted_upper = as.vector(pred_sum$upper),

        stringsAsFactors = FALSE
      )

      return(out)
    }

    # ----------------------------------------------------------
    # SUMMARY acumulado por lag_periods
    # ----------------------------------------------------------
    if (!all(c("period", "lag_start", "lag_end") %in% names(lag_periods))) {
      stop("`lag_periods` must contain columns: period, lag_start, lag_end.")
    }

    # ==========================================================
    # ✅ AJUSTE CRÍTICO — BUG FIX
    # Antes, `n_at` podia não existir aqui
    # ==========================================================
    n_at <- length(at_vals)

    period_list <- vector("list", nrow(lag_periods))

    summarize_mat_draws <- function(mat_draws) {
      list(
        est   = apply(mat_draws, 1, stats::median, na.rm = TRUE),
        sd    = apply(mat_draws, 1, safe_sd),
        lower = apply(mat_draws, 1, function(z) safe_quantile(z)[1]),
        upper = apply(mat_draws, 1, function(z) safe_quantile(z)[2])
      )
    }

    for (pp in seq_len(nrow(lag_periods))) {

      sel <- which(
        lag_idx >= lag_periods$lag_start[pp] &
          lag_idx <= lag_periods$lag_end[pp]
      )

      eta_draws_mat <- sapply(eta_list, function(m) rowSums(m[, sel, drop = FALSE], na.rm = TRUE))
      if (is.null(dim(eta_draws_mat))) {
        eta_draws_mat <- matrix(eta_draws_mat, ncol = 1)
      }

      eta_sum <- summarize_mat_draws(eta_draws_mat)

      eff_draws_mat <- apply(eta_draws_mat, 2, transform_effect)
      if (is.null(dim(eff_draws_mat))) eff_draws_mat <- matrix(eff_draws_mat, ncol = 1)
      eff_sum <- summarize_mat_draws(eff_draws_mat)

      if (is.na(baseline_response)) {
        delta_sum <- list(
          est   = rep(NA_real_, n_at),
          sd    = rep(NA_real_, n_at),
          lower = rep(NA_real_, n_at),
          upper = rep(NA_real_, n_at)
        )
        delta_pp_sum <- delta_sum
        pred_sum <- delta_sum
        baseline_vec <- rep(NA_real_, n_at)
      } else {
        pred_draws_mat <- linkinv(linkfun(baseline_response) + eta_draws_mat)
        delta_draws_mat <- pred_draws_mat - baseline_response
        delta_pp_draws_mat <- 100 * delta_draws_mat

        pred_sum <- summarize_mat_draws(pred_draws_mat)
        delta_sum <- summarize_mat_draws(delta_draws_mat)
        delta_pp_sum <- summarize_mat_draws(delta_pp_draws_mat)
        baseline_vec <- rep(baseline_response, length(at_vals))
      }

      period_list[[pp]] <- data.frame(
        var = var_one,
        period = lag_periods$period[pp],
        scale = "accumulated",
        value = at_vals,

        eta = eta_sum$est,
        eta_sd = eta_sum$sd,
        eta_lower = eta_sum$lower,
        eta_upper = eta_sum$upper,

        effect = eff_sum$est,
        effect_sd = eff_sum$sd,
        effect_lower = eff_sum$lower,
        effect_upper = eff_sum$upper,

        delta = delta_sum$est,
        delta_sd = delta_sum$sd,
        delta_lower = delta_sum$lower,
        delta_upper = delta_sum$upper,

        delta_pp = delta_pp_sum$est,
        delta_pp_sd = delta_pp_sum$sd,
        delta_pp_lower = delta_pp_sum$lower,
        delta_pp_upper = delta_pp_sum$upper,

        baseline = baseline_vec,
        predicted = pred_sum$est,
        predicted_sd = pred_sum$sd,
        predicted_lower = pred_sum$lower,
        predicted_upper = pred_sum$upper,

        stringsAsFactors = FALSE
      )
    }

    return(dplyr::bind_rows(period_list))
  }

  # ----------------------------------------------------------
  # Apply over one or many variables
  # ----------------------------------------------------------
  out_list <- lapply(vars, .summarise_effect_one_var)
  out <- do.call(rbind, out_list)

  rownames(out) <- NULL
  out
}
