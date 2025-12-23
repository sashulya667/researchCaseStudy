# ============================================================
# analysis/03_dr_greg.R
# Sequential run: IPW baseline + DR-GREG (mean) + DR-CDF (median/AROP)
# for Dataset_B and Dataset_E.
#
# Uses:
# - R/ipw.R (selection model + overlap + helpers)
# - R/dr_greg.R (this file’s functions)
# ============================================================

devtools::load_all()
rm(list = ls())

library(tidyverse)
library(here)

# -------------------------
# Settings
# -------------------------
X_VARS <- c("REG", "UEP", "DOU")
INCOME_VAR <- "EDI"
POVERTY_FRAC <- 0.6

# selection model RHS (keep simple + stable)
RHS_X <- ~ REG + UEP + DOU

# outcome model RHS for DR (keep same for clarity/defensibility)
RHS_Y <- ~ REG + UEP + DOU

# CDF grid controls
N_GRID <- 121  # number of thresholds for DR CDF grid
GRID_TRIM_Q <- c(0.01, 0.99)  # build grid from these quantiles of Y

SAVE_PLOTS <- TRUE
OUT_DIR <- here("outputs", "dr_greg_stepwise")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Load
# -------------------------
load_data <- function(filename) read.csv(file.path(here("data"), filename))

df_B_raw <- load_data("Dataset_B.csv")
df_E_raw <- load_data("Dataset_E.csv")
ref_raw  <- load_data("reference_sample_probabilistic.csv")

msg("Loaded rows: B=", nrow(df_B_raw), " | E=", nrow(df_E_raw), " | Ref=", nrow(ref_raw))

# -------------------------
# Harmonize factor levels
# -------------------------
lev_REG <- sort(unique(c(df_B_raw$REG, df_E_raw$REG, ref_raw$REG)))
lev_UEP <- sort(unique(c(df_B_raw$UEP, df_E_raw$UEP, ref_raw$UEP)))
lev_DOU <- sort(unique(c(df_B_raw$DOU, df_E_raw$DOU, ref_raw$DOU)))

df_B <- df_B_raw |>
  mutate(REG = factor(REG, levels = lev_REG),
         UEP = factor(UEP, levels = lev_UEP),
         DOU = factor(DOU, levels = lev_DOU))

df_E <- df_E_raw |>
  mutate(REG = factor(REG, levels = lev_REG),
         UEP = factor(UEP, levels = lev_UEP),
         DOU = factor(DOU, levels = lev_DOU))

ref <- ref_raw |>
  mutate(REG = factor(REG, levels = lev_REG),
         UEP = factor(UEP, levels = lev_UEP),
         DOU = factor(DOU, levels = lev_DOU))

# -------------------------
# Drop missing X
# -------------------------
before_B <- nrow(df_B); df_B <- df_B |> filter(!if_any(all_of(X_VARS), is.na)); msg("B drop_missing_x: dropped ", before_B - nrow(df_B))
before_E <- nrow(df_E); df_E <- df_E |> filter(!if_any(all_of(X_VARS), is.na)); msg("E drop_missing_x: dropped ", before_E - nrow(df_E))
before_R <- nrow(ref ); ref  <- ref  |> filter(!if_any(all_of(X_VARS), is.na)); msg("Ref drop_missing_x: dropped ", before_R - nrow(ref))

# ============================================================
# Helper: run one dataset end-to-end
# ============================================================
run_one <- function(nonprob_df, ref_df, dataset_name) {
  msg("\n==============================")
  msg("RUN: ", dataset_name)
  msg("==============================")

  # ---- Build training data
  train <- make_ipw_training_data(nonprob_df, ref_df, ref_weight_var = "design_weight")
  combined <- train$combined
  msg("Training rows: nonprob=", train$n_nonprob, " | ref=", train$n_ref)

  # ---- Overlap + drop unsupported
  ov <- compute_overlap(combined, x_vars = X_VARS, s_var = "S", w_var = "w_fit")
  msg("Overlap summary:"); print(ov$summary)

  bad_cells <- ov$bad_no_ref |> select(all_of(X_VARS)) |> distinct()
  drop_obj <- drop_bad_cells(nonprob_df, bad_cells, x_vars = X_VARS)
  nonprob_support <- drop_obj$df
  combined_fit <- combined |> anti_join(bad_cells, by = X_VARS)
  msg("Dropped nonprob rows (no ref support): ", drop_obj$n_dropped)

  # ---- Fit GROUPED selection model
  cell <- make_cell_counts(combined_fit, x_vars = X_VARS, s_var = "S", weight_var = "w_fit")
  fit_sel <- fit_selection_logit_grouped(cell, RHS_X)
  msg("Selection GROUPED glm: converged=", fit_sel$converged, " iter=", fit_sel$iter)

  # ---- Predict p_hat on supported nonprob, compute density-ratio weights
  p_np <- predict_p_hat(fit_sel, nonprob_support)
  w_norm <- dr_weights_np(p_np, eps = 1e-6)  # sums to 1

  # Attach for convenience
  r_raw <- density_ratio_from_p(p_np, eps = 1e-6)
  nonprob_out <- nonprob_support |>
    mutate(
      p_hat = as.numeric(p_np),
      r = as.numeric(r_raw),
      w_norm = as.numeric(w_norm)
    )

  # Diagnostics
  p_q <- quantile(nonprob_out$p_hat, c(0, .001, .01, .5, .99, .999, 1), na.rm = TRUE)
  msg("p_hat quantiles:"); print(p_q)

  w_summ <- weight_summary(nonprob_out$w_norm)
  msg("w_norm summary (should sum to 1, small sd):"); print(w_summ)
  msg("ESS using raw r (informal): ", effective_sample_size(nonprob_out$r))

  # ---- IPW baseline AROP (using normalized r-weights)
  y_np <- nonprob_out[[INCOME_VAR]]
  ipw_arop_tbl <- ipw_arop(y_np, nonprob_out$w_norm, poverty_frac = POVERTY_FRAC) |>
    mutate(dataset = dataset_name) |>
    select(dataset, everything())

  # ---- DR mean (optional, for report)
  dr_mean_tbl <- dr_mean(
    nonprob_df = nonprob_out,
    ref_df = ref_df,
    y_var = INCOME_VAR,
    ref_weight_var = "design_weight",
    w_np_norm = nonprob_out$w_norm,
    rhs_x = RHS_Y
  ) |>
    mutate(dataset = dataset_name) |>
    select(dataset, everything())

  # ---- DR CDF grid + DR AROP via CDF
  # Build grid from nonprob Y quantiles
  qlo <- unname(quantile(y_np, GRID_TRIM_Q[1], na.rm = TRUE))
  qhi <- unname(quantile(y_np, GRID_TRIM_Q[2], na.rm = TRUE))
  t_grid <- seq(qlo, qhi, length.out = N_GRID)

  cdf_tbl <- dr_cdf_grid(
    nonprob_df = nonprob_out,
    ref_df = ref_df,
    y_var = INCOME_VAR,
    ref_weight_var = "design_weight",
    w_np_norm = nonprob_out$w_norm,
    rhs_x = RHS_Y,
    t_grid = t_grid
  ) |>
    mutate(dataset = dataset_name)

  dr_arop_tbl <- dr_arop_via_cdf(cdf_tbl, poverty_frac = POVERTY_FRAC) |>
    mutate(dataset = dataset_name) |>
    select(dataset, everything())

  # Optional plot: compare IPW vs DR CDF
  p_cdf <- NULL
  if (isTRUE(SAVE_PLOTS)) {
    p_cdf <- ggplot(cdf_tbl, aes(x = t)) +
      geom_line(aes(y = F_ipw), linewidth = 0.8) +
      geom_line(aes(y = F_dr), linewidth = 0.8) +
      labs(
        title = paste0("CDF: IPW vs DR — ", dataset_name),
        x = "t",
        y = "F(t)"
      ) +
      theme_minimal()
  }

  list(
    dataset = dataset_name,
    nonprob_out = nonprob_out,
    fit_sel = fit_sel,
    overlap_summary = ov$summary,
    bad_cells = bad_cells,
    ipw_arop = ipw_arop_tbl,
    dr_mean = dr_mean_tbl,
    cdf = cdf_tbl,
    dr_arop = dr_arop_tbl,
    plot_cdf = p_cdf
  )
}

# ============================================================
# Run both datasets
# ============================================================
res_B <- run_one(df_B, ref, "Dataset_B")
res_E <- run_one(df_E, ref, "Dataset_E")

# ============================================================
# Collect tables
# ============================================================
tbl_ipw_arop <- bind_rows(res_B$ipw_arop, res_E$ipw_arop)
tbl_dr_mean  <- bind_rows(res_B$dr_mean,  res_E$dr_mean)
tbl_dr_arop  <- bind_rows(res_B$dr_arop,  res_E$dr_arop)

msg("\nIPW AROP (baseline):"); print(tbl_ipw_arop)
msg("\nDR mean (reporting):"); print(tbl_dr_mean)
msg("\nAROP via CDF (IPW_CDF vs DR_CDF):"); print(tbl_dr_arop)

# ============================================================
# Save outputs
# ============================================================
saveRDS(list(
  res_B = res_B,
  res_E = res_E,
  ipw_arop = tbl_ipw_arop,
  dr_mean = tbl_dr_mean,
  dr_arop = tbl_dr_arop
), file.path(OUT_DIR, "dr_greg_results.rds"))

readr::write_csv(tbl_ipw_arop, file.path(OUT_DIR, "ipw_arop_baseline.csv"))
readr::write_csv(tbl_dr_mean,  file.path(OUT_DIR, "dr_mean.csv"))
readr::write_csv(tbl_dr_arop,  file.path(OUT_DIR, "arop_via_cdf_ipw_vs_dr.csv"))

if (isTRUE(SAVE_PLOTS)) {
  ggsave(file.path(OUT_DIR, "cdf_ipw_vs_dr_dataset_B.png"), res_B$plot_cdf, width = 8, height = 5, dpi = 160)
  ggsave(file.path(OUT_DIR, "cdf_ipw_vs_dr_dataset_E.png"), res_E$plot_cdf, width = 8, height = 5, dpi = 160)
}

msg("\nSaved outputs to: ", OUT_DIR)

