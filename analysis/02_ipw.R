# ============================================================
# analysis/02_ipw.R
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

SAVE_PLOTS <- FALSE
OUT_DIR <- here("outputs", "ipw_stepwise")
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
# Harmonize factor levels across all sources
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
# Explicit missing-X removal with prints (B, E, Ref)
# -------------------------
before_B <- nrow(df_B); df_B <- df_B |> filter(!if_any(all_of(X_VARS), is.na)); msg("B drop_missing_x: dropped ", before_B - nrow(df_B))
before_E <- nrow(df_E); df_E <- df_E |> filter(!if_any(all_of(X_VARS), is.na)); msg("E drop_missing_x: dropped ", before_E - nrow(df_E))
before_R <- nrow(ref ); ref  <- ref  |> filter(!if_any(all_of(X_VARS), is.na)); msg("Ref drop_missing_x: dropped ", before_R - nrow(ref))

# Selection model formulas
FORMULA_S <- S ~ REG + UEP + DOU
RHS_X <- ~ REG + UEP + DOU

# ============================================================
# STEP 1) Build training data for BOTH datasets
# ============================================================
msg("\nSTEP 1 | Build combined training tables")
train_B <- make_ipw_training_data(df_B, ref, ref_weight_var = "design_weight")
train_E <- make_ipw_training_data(df_E, ref, ref_weight_var = "design_weight")

combined_B <- train_B$combined
combined_E <- train_E$combined

msg("Training rows: B nonprob=", train_B$n_nonprob, " | B ref=", train_B$n_ref)
msg("Training rows: E nonprob=", train_E$n_nonprob, " | E ref=", train_E$n_ref)

# ============================================================
# STEP 2) Overlap diagnostics for BOTH datasets
# ============================================================
msg("\nSTEP 2 | Overlap diagnostics (REG×UEP×DOU)")
ov_B <- compute_overlap(combined_B, x_vars = X_VARS, s_var = "S", w_var = "w_fit")
ov_E <- compute_overlap(combined_E, x_vars = X_VARS, s_var = "S", w_var = "w_fit")

msg("Overlap summary B: "); print(ov_B$summary)
msg("Overlap summary E: "); print(ov_E$summary)

bad_cells_B <- ov_B$bad_no_ref |> select(all_of(X_VARS)) |> distinct()
bad_cells_E <- ov_E$bad_no_ref |> select(all_of(X_VARS)) |> distinct()

# ============================================================
# STEP 3) Drop unsupported cells in nonprob and filter combined_*_fit
# ============================================================
msg("\nSTEP 3 | Drop nonprob rows in no-ref-support cells (B and E)")

drop_B <- drop_bad_cells(df_B, bad_cells_B, x_vars = X_VARS)
drop_E <- drop_bad_cells(df_E, bad_cells_E, x_vars = X_VARS)

df_B_support <- drop_B$df
df_E_support <- drop_E$df

combined_B_fit <- combined_B |> anti_join(bad_cells_B, by = X_VARS)
combined_E_fit <- combined_E |> anti_join(bad_cells_E, by = X_VARS)

msg("Dropped B nonprob rows (no ref support): ", drop_B$n_dropped)
msg("Dropped E nonprob rows (no ref support): ", drop_E$n_dropped)

# ============================================================
# STEP 4) Fit TWO selection models for BOTH datasets
#   (A) row-level weighted logistic
#   (B) grouped-cell logistic
# ============================================================
msg("\nSTEP 4 | Fit selection models (ROW and GROUPED)")

# ---- Row level ----
fit_B_row <- fit_selection_logit_row(combined_B_fit, FORMULA_S, weight_var = "w_fit")
fit_E_row <- fit_selection_logit_row(combined_E_fit, FORMULA_S, weight_var = "w_fit")
msg("B row glm: converged=", fit_B_row$converged, " iter=", fit_B_row$iter)
msg("E row glm: converged=", fit_E_row$converged, " iter=", fit_E_row$iter)

diag_B_row <- diagnose_row_glm(fit_B_row, combined_B_fit, x_vars = X_VARS, label = "Dataset_B ROW")
diag_E_row <- diagnose_row_glm(fit_E_row, combined_E_fit, x_vars = X_VARS, label = "Dataset_E ROW")

# Weighted class balance (THIS usually exposes the problem)
with(combined_B_fit, table(S))
with(combined_B_fit, tapply(w_fit, S, sum))                 # total weight by class
with(combined_B_fit, sum(w_fit * S) / sum(w_fit))           # weighted mean of S

g <- get_glm(fit_B_row)

# 1) What intercept should roughly be (from weighted mean of S)?
p_bar <- with(combined_B_fit, sum(w_fit * S) / sum(w_fit))
c(p_bar = p_bar, logit_p_bar = log(p_bar/(1-p_bar)))

# 2) Does the model reproduce the weighted mean? (it SHOULD, approximately)
c(
  obs_wmean_S = with(combined_B_fit, weighted.mean(S, w_fit)),
  fit_wmean_p = with(combined_B_fit, weighted.mean(fitted(g), w_fit))
)

# 3) Are fitted probabilities basically constant?
summary(fitted(g))
sd(fitted(g))

# 4) Weight pathologies (single gigantic weights can wreck IRLS)
summary(combined_B_fit$w_fit)
quantile(combined_B_fit$w_fit, c(.5,.9,.99,.999,1))
max(combined_B_fit$w_fit)

# 5) Show the non-identifiability explicitly
alias(g)$Complete

# ---- Cell level ----

cell_B    <- make_cell_counts(combined_B_fit, x_vars = X_VARS, s_var = "S", weight_var = "w_fit")
fit_B_grp <- fit_selection_logit_grouped(cell_B, RHS_X)

cell_E    <- make_cell_counts(combined_E_fit, x_vars = X_VARS, s_var = "S", weight_var = "w_fit")
fit_E_grp <- fit_selection_logit_grouped(cell_E, RHS_X)

msg("B grp glm: converged=", fit_B_grp$converged, " iter=", fit_B_grp$iter)
msg("E grp glm: converged=", fit_E_grp$converged, " iter=", fit_E_grp$iter)


# ============================================================
# STEP 5) Overlap plots (optional)
# ============================================================
if (isTRUE(SAVE_PLOTS)) {
  msg("\nSTEP 5 | Overlap plots")

  p_all_B_row <- predict_p_hat(fit_B_row, combined_B_fit)
  p_all_B_grp <- predict_p_hat(fit_B_grp, combined_B_fit)
  print(plot_overlap_p(p_all_B_row, combined_B_fit$S, "Overlap — Dataset_B (ROW)"))
  print(plot_overlap_p(p_all_B_grp, combined_B_fit$S, "Overlap — Dataset_B (GROUPED)"))

  p_all_E_row <- predict_p_hat(fit_E_row, combined_E_fit)
  p_all_E_grp <- predict_p_hat(fit_E_grp, combined_E_fit)
  print(plot_overlap_p(p_all_E_row, combined_E_fit$S, "Overlap — Dataset_E (ROW)"))
  print(plot_overlap_p(p_all_E_grp, combined_E_fit$S, "Overlap — Dataset_E (GROUPED)"))
}

# ============================================================
# STEP 6) Predict p_hat on supported nonprob and compute IPW (B and E)
# ============================================================
msg("\nSTEP 6 | Predict p_hat + compute IPW (ROW vs GROUPED)")

# ---- B ----
p_B_row <- predict_p_hat(fit_B_row, df_B_support)
p_B_grp <- predict_p_hat(fit_B_grp, df_B_support)

msg("B p_hat ROW quantiles:"); print(quantile(p_B_row, c(0,.001,.01,.5,.99,.999,1), na.rm=TRUE))
msg("B p_hat GRP quantiles:"); print(quantile(p_B_grp, c(0,.001,.01,.5,.99,.999,1), na.rm=TRUE))

w_B_row <- ipw_from_p(p_B_row)
w_B_grp <- ipw_from_p(p_B_grp)

dfB_ipw_row <- attach_ipw(df_B_support, p_B_row, w_B_row)
dfB_ipw_grp <- attach_ipw(df_B_support, p_B_grp, w_B_grp)

# ---- E ----
p_E_row <- predict_p_hat(fit_E_row, df_E_support)
p_E_grp <- predict_p_hat(fit_E_grp, df_E_support)

msg("E p_hat ROW quantiles:"); print(quantile(p_E_row, c(0,.001,.01,.5,.99,.999,1), na.rm=TRUE))
msg("E p_hat GRP quantiles:"); print(quantile(p_E_grp, c(0,.001,.01,.5,.99,.999,1), na.rm=TRUE))

w_E_row <- ipw_from_p(p_E_row)
w_E_grp <- ipw_from_p(p_E_grp)

dfE_ipw_row <- attach_ipw(df_E_support, p_E_row, w_E_row)
dfE_ipw_grp <- attach_ipw(df_E_support, p_E_grp, w_E_grp)


# ============================================================
# STEP 7) Weight diagnostics (B and E, ROW vs GROUPED)
# ============================================================
msg("\nSTEP 7 | Weight diagnostics")

msg("\nDataset_B weights (ROW)"); print(weight_summary(dfB_ipw_row$w_ipw))
msg("ESS B (ROW) = ", effective_sample_size(dfB_ipw_row$w_ipw))

msg("\nDataset_B weights (GROUPED)"); print(weight_summary(dfB_ipw_grp$w_ipw))
msg("ESS B (GROUPED) = ", effective_sample_size(dfB_ipw_grp$w_ipw))

msg("\nDataset_E weights (ROW)"); print(weight_summary(dfE_ipw_row$w_ipw))
msg("ESS E (ROW) = ", effective_sample_size(dfE_ipw_row$w_ipw))

msg("\nDataset_E weights (GROUPED)"); print(weight_summary(dfE_ipw_grp$w_ipw))
msg("ESS E (GROUPED) = ", effective_sample_size(dfE_ipw_grp$w_ipw))

# ============================================================
# STEP 8) AROP (B and E, ROW vs GROUPED)
# ============================================================
msg("\nSTEP 8 | Compute AROP (ROW vs GROUPED)")

arop_B_row <- bind_rows(
  arop_unweighted(dfB_ipw_row, INCOME_VAR, POVERTY_FRAC),
  arop_weighted(dfB_ipw_row, INCOME_VAR, "w_ipw", POVERTY_FRAC)
) |> mutate(dataset="Dataset_B", method="ROW") |> select(dataset, method, everything())

arop_B_grp <- bind_rows(
  arop_unweighted(dfB_ipw_grp, INCOME_VAR, POVERTY_FRAC),
  arop_weighted(dfB_ipw_grp, INCOME_VAR, "w_ipw", POVERTY_FRAC)
) |> mutate(dataset="Dataset_B", method="GROUPED") |> select(dataset, method, everything())

arop_E_row <- bind_rows(
  arop_unweighted(dfE_ipw_row, INCOME_VAR, POVERTY_FRAC),
  arop_weighted(dfE_ipw_row, INCOME_VAR, "w_ipw", POVERTY_FRAC)
) |> mutate(dataset="Dataset_E", method="ROW") |> select(dataset, method, everything())

arop_E_grp <- bind_rows(
  arop_unweighted(dfE_ipw_grp, INCOME_VAR, POVERTY_FRAC),
  arop_weighted(dfE_ipw_grp, INCOME_VAR, "w_ipw", POVERTY_FRAC)
) |> mutate(dataset="Dataset_E", method="GROUPED") |> select(dataset, method, everything())

arop_all <- bind_rows(arop_B_row, arop_B_grp, arop_E_row, arop_E_grp)

msg("\nAROP table:"); print(arop_all)

# ============================================================
# STEP 9) Balance tables (Ref vs Nonprob vs IPW), no extra functions defined
# ============================================================
msg("\nSTEP 9 | Balance tables (Ref vs Nonprob vs IPW)")

weighted_prop_1d <- function(df, var, weight_var = NULL) {
  v <- df[[var]]
  ok <- !is.na(v)
  if (is.null(weight_var)) {
    out <- df |> filter(ok) |> count(.data[[var]], name = "n") |>
      mutate(prop = n / sum(n))
  } else {
    w <- df[[weight_var]]
    ok <- ok & is.finite(w) & !is.na(w) & w > 0
    out <- df |> filter(ok) |> count(.data[[var]], wt = .data[[weight_var]], name = "w_sum") |>
      mutate(prop = w_sum / sum(w_sum))
  }
  names(out)[1] <- var
  out |> select(all_of(var), prop)
}

compare_marginal_simple <- function(ref_df, nonprob_df, var, ref_w = "design_weight", ipw_w = "w_ipw") {
  pr_ref <- weighted_prop_1d(ref_df, var, ref_w) |> rename(ref_prop = prop)
  pr_np  <- weighted_prop_1d(nonprob_df, var, NULL) |> rename(np_unw_prop = prop)
  pr_ipw <- weighted_prop_1d(nonprob_df, var, ipw_w) |> rename(np_ipw_prop = prop)

  out <- pr_ref |>
    full_join(pr_np,  by = var) |>
    full_join(pr_ipw, by = var) |>
    mutate(
      ref_prop    = tidyr::replace_na(ref_prop, 0),
      np_unw_prop = tidyr::replace_na(np_unw_prop, 0),
      np_ipw_prop = tidyr::replace_na(np_ipw_prop, 0),
      diff_unw = np_unw_prop - ref_prop,
      diff_ipw = np_ipw_prop - ref_prop
    ) |>
    arrange(desc(ref_prop))

  summ <- tibble::tibble(
    variable = var,
    level = "MAX_ABS_DIFF",
    ref_prop = NA_real_,
    np_unw_prop = NA_real_,
    np_ipw_prop = NA_real_,
    diff_unw = max(abs(out$diff_unw), na.rm = TRUE),
    diff_ipw = max(abs(out$diff_ipw), na.rm = TRUE)
  )

  out2 <- out |>
    mutate(variable = var, level = as.character(.data[[var]])) |>
    select(variable, level, ref_prop, np_unw_prop, np_ipw_prop, diff_unw, diff_ipw)

  bind_rows(out2, summ)
}

make_balance_table <- function(ref_df, df_ipw, ref_w = "design_weight", ipw_w = "w_ipw") {
  bind_rows(
    compare_marginal_simple(ref_df, df_ipw, "REG", ref_w, ipw_w),
    compare_marginal_simple(ref_df, df_ipw, "UEP", ref_w, ipw_w),
    compare_marginal_simple(ref_df, df_ipw, "DOU", ref_w, ipw_w)
  ) |>
    mutate(across(c(ref_prop, np_unw_prop, np_ipw_prop, diff_unw, diff_ipw), ~ round(.x, 6)))
}

balance_B_row <- make_balance_table(ref, dfB_ipw_row)
balance_B_grp <- make_balance_table(ref, dfB_ipw_grp)
balance_E_row <- make_balance_table(ref, dfE_ipw_row)
balance_E_grp <- make_balance_table(ref, dfE_ipw_grp)

msg("\nBalance — Dataset_B (ROW):"); print(balance_B_row)
msg("\nBalance — Dataset_B (GROUPED):"); print(balance_B_grp)
msg("\nBalance — Dataset_E (ROW):"); print(balance_E_row)
msg("\nBalance — Dataset_E (GROUPED):"); print(balance_E_grp)

# ============================================================
# Save outputs
# ============================================================
saveRDS(list(
  dfB_ipw_row = dfB_ipw_row,
  dfB_ipw_grp = dfB_ipw_grp,
  dfE_ipw_row = dfE_ipw_row,
  dfE_ipw_grp = dfE_ipw_grp,
  fit_B_row = fit_B_row, fit_B_grp = fit_B_grp,
  fit_E_row = fit_E_row, fit_E_grp = fit_E_grp,
  arop_all = arop_all,
  balance_B_row = balance_B_row, balance_B_grp = balance_B_grp,
  balance_E_row = balance_E_row, balance_E_grp = balance_E_grp,
  overlap_B = ov_B$summary, overlap_E = ov_E$summary
), file.path(OUT_DIR, "ipw_stepwise_row_vs_grouped.rds"))

readr::write_csv(arop_all, file.path(OUT_DIR, "arop_row_vs_grouped.csv"))

msg("\nSaved: ",
    file.path(OUT_DIR, "ipw_stepwise_row_vs_grouped.rds"),
    " and ",
    file.path(OUT_DIR, "arop_row_vs_grouped.csv"))

