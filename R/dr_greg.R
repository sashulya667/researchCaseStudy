# ============================================================
# R/dr_greg.R
# Doubly robust (DR-GREG / AIPW) estimators for nonprob + probability reference
#
# Requires your existing R/ipw.R to be loaded (for):
#   make_ipw_training_data(), compute_overlap(), drop_bad_cells(),
#   make_cell_counts(), fit_selection_logit_grouped(), predict_p_hat(),
#   attach_ipw(), effective_sample_size(), msg()
#
# Core setup:
# - Nonprob sample: S=1, has outcome Y (e.g., EDI)
# - Reference sample: S=0, has design weights d_i = 1/pi_i (design_weight)
# - Shared covariates X (e.g., REG, UEP, DOU)
#
# Selection model:
#   p(x) = P(S=1 | X=x)
# Density ratio up to constant:
#   r(x) ∝ (1 - p(x)) / p(x)
#
# DR mean:
#   mu_DR = E_ref[ mhat(X) ] + E_np,r[ Y - mhat(X) ]
#
# DR CDF at threshold t:
#   F_DR(t) = E_ref[ mhat_t(X) ] + E_np,r[ 1{Y<=t} - mhat_t(X) ]
#
# Then DR median solves F_DR(med)=0.5, and DR AROP = F_DR(0.6*med).
# ============================================================

# -------------------------
# Small utilities
# -------------------------

assert_required_cols <- function(df, cols, fn_name = "function") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing required columns: %s",
                 fn_name, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# Normalize weights to sum to 1
normalize_weights <- function(w) {
  w <- as.numeric(w)
  ok <- is.finite(w) & !is.na(w) & w > 0
  if (!any(ok)) stop("normalize_weights: no positive finite weights", call. = FALSE)
  w2 <- w
  w2[!ok] <- NA_real_
  s <- sum(w2, na.rm = TRUE)
  w2 / s
}

# Weighted mean with sanity checks
wmean_safe <- function(x, w) {
  x <- as.numeric(x); w <- as.numeric(w)
  ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w) & w > 0
  if (!any(ok)) stop("wmean_safe: no valid (x,w)", call. = FALSE)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

# Weighted quantile (used for IPW median)
weighted_quantile <- function(x, w, p = 0.5) {
  x <- as.numeric(x); w <- as.numeric(w)
  ok <- is.finite(x) & !is.na(x) & is.finite(w) & !is.na(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) stop("weighted_quantile: no valid data", call. = FALSE)
  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w)
  x[which(cw >= p * sum(w))[1]]
}

# -------------------------
# Density-ratio weights from p_hat
# -------------------------

# r(x) ∝ (1-p)/p (clipped)
density_ratio_from_p <- function(p_hat, eps = 1e-6) {
  p <- pmin(pmax(as.numeric(p_hat), eps), 1 - eps)
  (1 - p) / p
}

# Get normalized r-weights on the nonprob sample
# These approximate the population distribution over X when applied to nonprob.
dr_weights_np <- function(p_hat_np, eps = 1e-6) {
  r <- density_ratio_from_p(p_hat_np, eps = eps)
  normalize_weights(r)
}

# -------------------------
# Simple IPW estimators (baseline)
# -------------------------

ipw_mean <- function(y, w_norm) {
  wmean_safe(y, w_norm)
}

ipw_cdf_at <- function(y, w_norm, t) {
  y <- as.numeric(y)
  ind <- as.numeric(y <= t)
  wmean_safe(ind, w_norm)
}

ipw_median <- function(y, w_norm) {
  weighted_quantile(y, w_norm, p = 0.5)
}

ipw_arop <- function(y, w_norm, poverty_frac = 0.6) {
  med <- ipw_median(y, w_norm)
  line <- poverty_frac * med
  arop <- ipw_cdf_at(y, w_norm, line)
  tibble::tibble(
    method = "IPW",
    median = med,
    poverty_line = line,
    arop = arop
  )
}

# -------------------------
# Outcome models for DR
# -------------------------

# Fit outcome regression on nonprob for continuous Y
# Keep it simple/defensible: main effects by default.
fit_outcome_lm_np <- function(nonprob_df, y_var, rhs_x) {
  assert_required_cols(nonprob_df, c(y_var), "fit_outcome_lm_np")
  rhs <- paste(deparse(rhs_x[[2]]), collapse = "")
  f <- stats::as.formula(paste0(y_var, " ~ ", rhs))
  stats::lm(f, data = nonprob_df)
}

# Fit distribution regression (logistic) on nonprob for I(Y<=t)
fit_outcome_logit_cdf_np <- function(nonprob_df, y_var, rhs_x, t) {
  assert_required_cols(nonprob_df, c(y_var), "fit_outcome_logit_cdf_np")
  rhs <- paste(deparse(rhs_x[[2]]), collapse = "")
  nonprob_df2 <- nonprob_df
  nonprob_df2[[".__I__"]] <- as.integer(as.numeric(nonprob_df2[[y_var]]) <= t)
  f <- stats::as.formula(paste0(".__I__ ~ ", rhs))
  stats::glm(f, data = nonprob_df2, family = stats::binomial())
}

# Predict helper (respects xlevels like your predict_p_hat)
predict_response <- function(fit, df_new) {
  # align factor levels if present
  if (!is.null(fit$xlevels)) {
    for (v in names(fit$xlevels)) {
      if (v %in% names(df_new)) {
        df_new[[v]] <- factor(df_new[[v]], levels = fit$xlevels[[v]])
      }
    }
  }
  as.numeric(stats::predict(fit, newdata = df_new, type = "response"))
}

# -------------------------
# DR mean (continuous Y)
# -------------------------

dr_mean <- function(nonprob_df,
                    ref_df,
                    y_var,
                    ref_weight_var = "design_weight",
                    w_np_norm,
                    rhs_x) {

  assert_required_cols(nonprob_df, c(y_var), "dr_mean")
  assert_required_cols(ref_df, c(ref_weight_var), "dr_mean")

  # Outcome model on nonprob
  fit_m <- fit_outcome_lm_np(nonprob_df, y_var, rhs_x)

  m_np  <- predict_response(fit_m, nonprob_df)
  m_ref <- predict_response(fit_m, ref_df)

  # E_ref[m(X)]
  mu_or <- wmean_safe(m_ref, ref_df[[ref_weight_var]])

  # E_np,r[ Y - m(X) ]
  y_np <- as.numeric(nonprob_df[[y_var]])
  mu_corr <- wmean_safe(y_np - m_np, w_np_norm)

  mu_dr <- mu_or + mu_corr

  tibble::tibble(
    method = "DR_MEAN",
    mu_or = mu_or,
    mu_corr = mu_corr,
    mu_dr = mu_dr
  )
}

# -------------------------
# DR CDF grid + median + AROP
# -------------------------

# Compute DR CDF on a grid of thresholds t_grid.
# Also returns IPW CDF on same grid (baseline).
dr_cdf_grid <- function(nonprob_df,
                        ref_df,
                        y_var,
                        ref_weight_var = "design_weight",
                        w_np_norm,
                        rhs_x,
                        t_grid) {

  assert_required_cols(nonprob_df, c(y_var), "dr_cdf_grid")
  assert_required_cols(ref_df, c(ref_weight_var), "dr_cdf_grid")

  y_np <- as.numeric(nonprob_df[[y_var]])

  out <- vector("list", length(t_grid))
  for (k in seq_along(t_grid)) {
    t <- t_grid[[k]]

    # Fit m_t(x) = P(Y<=t | X) on nonprob
    fit_mt <- fit_outcome_logit_cdf_np(nonprob_df, y_var, rhs_x, t)

    mt_np  <- predict_response(fit_mt, nonprob_df)
    mt_ref <- predict_response(fit_mt, ref_df)

    # DR CDF
    F_or   <- wmean_safe(mt_ref, ref_df[[ref_weight_var]])
    ind_np <- as.numeric(y_np <= t)
    F_corr <- wmean_safe(ind_np - mt_np, w_np_norm)
    F_dr   <- F_or + F_corr

    # IPW baseline CDF
    F_ipw <- ipw_cdf_at(y_np, w_np_norm, t)

    out[[k]] <- tibble::tibble(
      t = t,
      F_ipw = F_ipw,
      F_or = F_or,
      F_corr = F_corr,
      F_dr = F_dr
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::arrange(.data$t)
}

# Interpolate a monotone-ish CDF to get quantile
cdf_quantile_from_grid <- function(cdf_tbl, p = 0.5, method = c("dr", "ipw")) {
  method <- match.arg(method)
  stopifnot(all(c("t", "F_dr", "F_ipw") %in% names(cdf_tbl)))

  F <- if (method == "dr") cdf_tbl$F_dr else cdf_tbl$F_ipw
  t <- cdf_tbl$t

  # Clip to [0,1] for safety
  F <- pmin(pmax(as.numeric(F), 0), 1)

  # Ensure increasing in a weak sense (monotone smoothing by cummax)
  F_mon <- cummax(F)

  # If p outside grid range, extrapolate by endpoints
  if (p <= min(F_mon, na.rm = TRUE)) return(min(t, na.rm = TRUE))
  if (p >= max(F_mon, na.rm = TRUE)) return(max(t, na.rm = TRUE))

  # Linear interpolation t(F) at level p
  stats::approx(x = F_mon, y = t, xout = p, ties = "ordered")$y
}

# Interpolate CDF value at arbitrary t0 from grid
cdf_value_from_grid <- function(cdf_tbl, t0, method = c("dr", "ipw")) {
  method <- match.arg(method)
  F <- if (method == "dr") cdf_tbl$F_dr else cdf_tbl$F_ipw
  t <- cdf_tbl$t

  F <- pmin(pmax(as.numeric(F), 0), 1)
  F_mon <- cummax(F)

  # interpolate F(t0)
  stats::approx(x = t, y = F_mon, xout = t0, rule = 2)$y
}

# DR AROP via DR median and DR CDF at poverty line
dr_arop_via_cdf <- function(cdf_tbl, poverty_frac = 0.6) {
  med_dr <- cdf_quantile_from_grid(cdf_tbl, p = 0.5, method = "dr")
  line_dr <- poverty_frac * med_dr
  arop_dr <- cdf_value_from_grid(cdf_tbl, line_dr, method = "dr")

  med_ipw <- cdf_quantile_from_grid(cdf_tbl, p = 0.5, method = "ipw")
  line_ipw <- poverty_frac * med_ipw
  arop_ipw <- cdf_value_from_grid(cdf_tbl, line_ipw, method = "ipw")

  tibble::tibble(
    method = c("IPW_CDF", "DR_CDF"),
    median = c(med_ipw, med_dr),
    poverty_line = c(line_ipw, line_dr),
    arop = c(arop_ipw, arop_dr)
  )
}
