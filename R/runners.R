# ============================================================
# R/runners.R
# Method runners that REUSE your existing functions from:
# - R/ipw.R
# - R/dr_greg.R
#
# Each runner takes (nonprob_df, ref_df, ...) and returns a tibble
# with: method, median, poverty_line, arop  (and optionally others).
# ============================================================

# small local assertion (keeps errors clean)
assert_required_cols2 <- function(df, cols, fn_name = "function") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing required columns: %s",
                 fn_name, paste(missing, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# Harmonize factor levels of X between nonprob and ref (safe + reproducible)
harmonize_x_levels <- function(nonprob_df, ref_df, x_vars = c("REG","UEP","DOU")) {
  out_np <- nonprob_df
  out_rf <- ref_df

  for (v in x_vars) {
    lev <- sort(unique(c(as.character(out_np[[v]]), as.character(out_rf[[v]]))))
    out_np[[v]] <- factor(out_np[[v]], levels = lev)
    out_rf[[v]] <- factor(out_rf[[v]], levels = lev)
  }
  list(nonprob = out_np, ref = out_rf)
}

# -------------------------
# Runner 1: IPW (GROUPED selection) + AROP via weighted median + binary
# This matches your "simple IPW then AROP (binary)"
# -------------------------
run_ipw_grouped_arop <- function(nonprob_df,
                                 ref_df,
                                 x_vars = c("REG","UEP","DOU"),
                                 income_var = "EDI",
                                 poverty_frac = 0.6,
                                 rhs_x = ~ REG + UEP + DOU,
                                 ref_weight_var = "design_weight") {

  assert_required_cols2(nonprob_df, c(x_vars, income_var), "run_ipw_grouped_arop(nonprob_df)")
  assert_required_cols2(ref_df, c(x_vars, ref_weight_var), "run_ipw_grouped_arop(ref_df)")

  # harmonize levels to avoid predict() NA
  h <- harmonize_x_levels(nonprob_df, ref_df, x_vars)
  nonprob_df <- h$nonprob
  ref_df <- h$ref

  # Build combined training
  train <- make_ipw_training_data(nonprob_df, ref_df, ref_weight_var = ref_weight_var)
  combined <- train$combined

  # overlap + drop unsupported
  ov <- compute_overlap(combined, x_vars = x_vars, s_var = "S", w_var = "w_fit")
  bad_cells <- ov$bad_no_ref |> dplyr::select(dplyr::all_of(x_vars)) |> dplyr::distinct()

  drop_obj <- drop_bad_cells(nonprob_df, bad_cells, x_vars = x_vars)
  nonprob_support <- drop_obj$df
  combined_fit <- combined |> dplyr::anti_join(bad_cells, by = x_vars)

  # fit grouped selection
  cell <- make_cell_counts(combined_fit, x_vars = x_vars, s_var = "S", weight_var = "w_fit")
  fit_sel <- fit_selection_logit_grouped(cell, rhs_x)

  # predict + compute ipw
  p_hat <- predict_p_hat(fit_sel, nonprob_support)
  w_ipw <- ipw_from_p(p_hat)

  df_ipw <- attach_ipw(nonprob_support, p_hat, w_ipw)
  out <- arop_weighted(df_ipw, income_var = income_var, weight_var = "w_ipw", poverty_frac = poverty_frac)

  out |>
    dplyr::mutate(method = "IPW_GROUPED_BOOL") |>
    dplyr::select(method, median, poverty_line, arop)
}

# -------------------------
# Runner 2: DR-CDF (plus IPW_CDF baseline on same grid)
# This is your "more correct with CDF"
# Returns 2 rows: IPW_CDF and DR_CDF
# -------------------------
run_dr_cdf_arop <- function(nonprob_df,
                            ref_df,
                            x_vars = c("REG","UEP","DOU"),
                            income_var = "EDI",
                            poverty_frac = 0.6,
                            rhs_x_sel = ~ REG + UEP + DOU,
                            rhs_y = ~ REG + UEP + DOU,
                            ref_weight_var = "design_weight",
                            n_grid = 121,
                            grid_trim_q = c(0.01, 0.99)) {

  assert_required_cols2(nonprob_df, c(x_vars, income_var), "run_dr_cdf_arop(nonprob_df)")
  assert_required_cols2(ref_df, c(x_vars, ref_weight_var), "run_dr_cdf_arop(ref_df)")

  # harmonize levels
  h <- harmonize_x_levels(nonprob_df, ref_df, x_vars)
  nonprob_df <- h$nonprob
  ref_df <- h$ref

  # Selection training
  train <- make_ipw_training_data(nonprob_df, ref_df, ref_weight_var = ref_weight_var)
  combined <- train$combined

  ov <- compute_overlap(combined, x_vars = x_vars, s_var = "S", w_var = "w_fit")
  bad_cells <- ov$bad_no_ref |> dplyr::select(dplyr::all_of(x_vars)) |> dplyr::distinct()

  drop_obj <- drop_bad_cells(nonprob_df, bad_cells, x_vars = x_vars)
  nonprob_support <- drop_obj$df
  combined_fit <- combined |> dplyr::anti_join(bad_cells, by = x_vars)

  # Fit grouped selection p(S=1|X)
  cell <- make_cell_counts(combined_fit, x_vars = x_vars, s_var = "S", weight_var = "w_fit")
  fit_sel <- fit_selection_logit_grouped(cell, rhs_x_sel)

  p_hat <- predict_p_hat(fit_sel, nonprob_support)
  w_norm <- dr_weights_np(p_hat, eps = 1e-6)

  nonprob_out <- nonprob_support |>
    dplyr::mutate(p_hat = as.numeric(p_hat), w_norm = as.numeric(w_norm))

  y_np <- as.numeric(nonprob_out[[income_var]])
  qlo <- unname(stats::quantile(y_np, grid_trim_q[1], na.rm = TRUE))
  qhi <- unname(stats::quantile(y_np, grid_trim_q[2], na.rm = TRUE))
  t_grid <- seq(qlo, qhi, length.out = n_grid)

  cdf_tbl <- dr_cdf_grid(
    nonprob_df = nonprob_out,
    ref_df = ref_df,
    y_var = income_var,
    ref_weight_var = ref_weight_var,
    w_np_norm = nonprob_out$w_norm,
    rhs_x = rhs_y,
    t_grid = t_grid
  )

  dr_arop_via_cdf(cdf_tbl, poverty_frac = poverty_frac) |>
    dplyr::select(method, median, poverty_line, arop)
}

# -------------------------
# Convenience: run all methods you want to compare
# -------------------------
run_all_methods <- function(nonprob_df,
                            ref_df,
                            x_vars = c("REG","UEP","DOU"),
                            income_var = "EDI",
                            poverty_frac = 0.6,
                            rhs_x_sel = ~ REG + UEP + DOU,
                            rhs_y = ~ REG + UEP + DOU,
                            ref_weight_var = "design_weight",
                            n_grid = 121,
                            grid_trim_q = c(0.01, 0.99)) {

  ipw_tbl <- run_ipw_grouped_arop(
    nonprob_df = nonprob_df,
    ref_df = ref_df,
    x_vars = x_vars,
    income_var = income_var,
    poverty_frac = poverty_frac,
    rhs_x = rhs_x_sel,
    ref_weight_var = ref_weight_var
  )

  dr_tbl <- run_dr_cdf_arop(
    nonprob_df = nonprob_df,
    ref_df = ref_df,
    x_vars = x_vars,
    income_var = income_var,
    poverty_frac = poverty_frac,
    rhs_x_sel = rhs_x_sel,
    rhs_y = rhs_y,
    ref_weight_var = ref_weight_var,
    n_grid = n_grid,
    grid_trim_q = grid_trim_q
  )

  dplyr::bind_rows(ipw_tbl, dr_tbl)
}
