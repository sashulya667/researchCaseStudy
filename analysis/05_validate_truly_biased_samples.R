# ============================================================
# analysis/validate_truly_biased_samples.R
# Validate methods on synthetic biased samples created by 04_create_biased_sample.
#
# Reads for each dataset:
# - <DATASET>__truth.csv
# - <DATASET>__ref_srs.csv
# - <DATASET>__np_bias_<VAR>__<strength>.csv  (mild/medium/strong)
#
# Compares:
# - UNWEIGHTED  (computed on SAME support as adjusted methods)
# - IPW_GROUPED_BOOL
# - IPW_CDF
# - DR_CDF
# ============================================================

devtools::load_all()
rm(list = ls())

library(tidyverse)
library(here)

msg <- function(...) cat(paste0(..., "\n"))

# -------------------------
# SETTINGS
# -------------------------
X_VARS <- c("REG","UEP","DOU")
INCOME_VAR <- "EDI"
POVERTY_FRAC <- 0.6

RHS_X_SEL <- ~ REG + UEP + DOU
RHS_Y     <- ~ REG + UEP + DOU

N_GRID <- 121
GRID_TRIM_Q <- c(0.01, 0.99)

DIR_IN <- here("outputs", "sim_biased_samples")  # change if needed
OUT_DIR <- here("outputs", "validate_truly_biased_samples")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

msg("DIR_IN  = ", DIR_IN)
msg("OUT_DIR = ", OUT_DIR)

# -------------------------
# File helpers
# -------------------------
path_truth <- function(dataset_name) file.path(DIR_IN, paste0(dataset_name, "__truth.csv"))
path_ref   <- function(dataset_name) file.path(DIR_IN, paste0(dataset_name, "__ref_srs.csv"))

list_np_bias_files <- function(dataset_name) {
  files <- list.files(
    DIR_IN,
    pattern = paste0("^", dataset_name, "__np_bias_.*\\.csv$"),
    full.names = TRUE
  )
  files <- files[grepl("__(mild|medium|strong)\\.csv$", files)]
  sort(files)
}

parse_bias_from_filename <- function(f) {
  b <- basename(f)
  parts <- strsplit(b, "__", fixed = TRUE)[[1]]
  bias_var <- sub("^np_bias_", "", parts[2])
  bias_strength <- sub("\\.csv$", "", parts[3])
  tibble(bias_var = bias_var, bias_strength = bias_strength)
}

# -------------------------
# Unweighted baseline
# -------------------------
compute_unweighted_arop <- function(df, income_var = INCOME_VAR, poverty_frac = POVERTY_FRAC) {
  y <- df[[income_var]]
  y <- y[is.finite(y) & !is.na(y)]
  med <- median(y)
  z <- poverty_frac * med
  arop <- mean(y < z)  # strict < for consistency
  tibble(
    method = "UNWEIGHTED",
    median = med,
    poverty_line = z,
    arop = arop
  )
}

# -------------------------
# Apply the SAME overlap trimming as in runners (harmonize first)
# Returns: list(df_support, n_dropped)
# -------------------------
trim_to_common_support_like_runners <- function(nonprob_df, ref_df,
                                                x_vars = X_VARS,
                                                ref_weight_var = "design_weight") {
  # harmonize factor levels exactly like runners
  h <- harmonize_x_levels(nonprob_df, ref_df, x_vars = x_vars)
  nonprob_h <- h$nonprob
  ref_h     <- h$ref

  train <- make_ipw_training_data(nonprob_h, ref_h, ref_weight_var = ref_weight_var)
  combined <- train$combined

  ov <- compute_overlap(combined, x_vars = x_vars, s_var = "S", w_var = "w_fit")
  bad_cells <- ov$bad_no_ref %>%
    dplyr::select(dplyr::all_of(x_vars)) %>%
    dplyr::distinct()

  drop_obj <- drop_bad_cells(nonprob_h, bad_cells, x_vars = x_vars)

  list(df_support = drop_obj$df, n_dropped = drop_obj$n_dropped)
}

# -------------------------
# One dataset runner
# -------------------------
validate_one_dataset <- function(dataset_name) {
  msg("\n==============================")
  msg("VALIDATE: ", dataset_name)
  msg("==============================")

  truth_file <- path_truth(dataset_name)
  ref_file   <- path_ref(dataset_name)

  if (!file.exists(truth_file)) stop("Missing truth file: ", truth_file, call. = FALSE)
  if (!file.exists(ref_file))   stop("Missing ref file: ", ref_file, call. = FALSE)

  truth <- readr::read_csv(truth_file, show_col_types = FALSE)
  ref   <- readr::read_csv(ref_file,   show_col_types = FALSE)

  stopifnot(all(c("median","poverty_line","arop") %in% names(truth)))

  files_np <- list_np_bias_files(dataset_name)
  if (length(files_np) == 0) stop("No biased nonprob files found for ", dataset_name, " in ", DIR_IN, call. = FALSE)

  msg("Nonprob files to validate:")
  walk2(seq_along(files_np), files_np, ~ msg(.x, " -> ", basename(.y)))

  out_all <- purrr::map_dfr(files_np, function(f_np) {
    nonprob <- readr::read_csv(f_np, show_col_types = FALSE)

    # sanity checks
    missing_x <- setdiff(X_VARS, names(nonprob))
    if (length(missing_x) > 0) {
      stop("Missing X vars in file: ", basename(f_np), " missing: ", paste(missing_x, collapse = ", "), call. = FALSE)
    }
    if (!INCOME_VAR %in% names(nonprob)) {
      stop("Missing income var ", INCOME_VAR, " in file: ", basename(f_np), call. = FALSE)
    }

    # --- trim for UNWEIGHTED to match target support ---
    trim <- trim_to_common_support_like_runners(nonprob, ref, x_vars = X_VARS, ref_weight_var = "design_weight")
    nonprob_support <- trim$df_support
    n_dropped <- trim$n_dropped

    est_unw <- compute_unweighted_arop(nonprob_support) %>%
      mutate(n_dropped_np = n_dropped)

    # --- adjusted methods (they do their own trimming internally) ---
    est_tbl <- run_all_methods(
      nonprob_df = nonprob,
      ref_df = ref,
      x_vars = X_VARS,
      income_var = INCOME_VAR,
      poverty_frac = POVERTY_FRAC,
      rhs_x_sel = RHS_X_SEL,
      rhs_y = RHS_Y,
      ref_weight_var = "design_weight",
      n_grid = N_GRID,
      grid_trim_q = GRID_TRIM_Q
    ) %>%
      mutate(n_dropped_np = n_dropped)

    est_tbl <- bind_rows(est_unw, est_tbl)

    info <- parse_bias_from_filename(f_np)

    est_tbl %>%
      mutate(
        dataset = dataset_name,
        file = basename(f_np),
        bias_var = info$bias_var,
        bias_strength = info$bias_strength,
        truth_median = truth$median[1],
        truth_poverty_line = truth$poverty_line[1],
        truth_arop = truth$arop[1],
        err_arop = arop - truth_arop,
        abs_err_arop = abs(err_arop),
        err_median = median - truth_median,
        abs_err_median = abs(err_median)
      ) %>%
      select(
        dataset, bias_var, bias_strength, method,
        median, poverty_line, arop,
        truth_median, truth_poverty_line, truth_arop,
        err_median, abs_err_median, err_arop, abs_err_arop,
        n_dropped_np,
        file
      )
  })

  out_all
}

# -------------------------
# RUN
# -------------------------
res_B <- validate_one_dataset("Dataset_B")
res_E <- validate_one_dataset("Dataset_E")
res_all <- bind_rows(res_B, res_E)

# ---- Diagnostics: how much got dropped? ----
msg("\nSupport trimming diagnostics (n_dropped_np):")
res_all %>%
  group_by(dataset, bias_var, bias_strength) %>%
  summarise(
    dropped_mean = mean(n_dropped_np, na.rm = TRUE),
    dropped_max  = max(n_dropped_np, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(dropped_max)) %>%
  print(n = 100)

# -------------------------
# Tables for paper
# -------------------------

# overall MAE (proportion)
tab_mae_overall <- res_all %>%
  group_by(method) %>%
  summarise(
    MAE_AROP = mean(abs_err_arop, na.rm = TRUE),
    MAE_median = mean(abs_err_median, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(MAE_AROP)

# overall MAE (pp)
tab_mae_overall_pp <- tab_mae_overall %>%
  mutate(MAE_AROP_pp = 100 * MAE_AROP) %>%
  select(method, MAE_AROP_pp, MAE_median)

# MAE by bias var/strength (pp), wide
tab_mae_by_strength_pp <- res_all %>%
  group_by(bias_var, bias_strength, method) %>%
  summarise(MAE_AROP_pp = 100 * mean(abs_err_arop, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = method, values_from = MAE_AROP_pp) %>%
  arrange(bias_var, factor(bias_strength, levels = c("mild","medium","strong")))

# -------------------------
# Optional: nicer labels for LaTeX
# -------------------------
pretty_method <- function(m) {
  recode(m,
         "UNWEIGHTED" = "Unweighted",
         "IPW_GROUPED_BOOL" = "IPW (grouped)",
         "IPW_CDF" = "IPW-CDF",
         "DR_CDF" = "DR-CDF")
}

tab_mae_overall_pp_pretty <- tab_mae_overall_pp %>%
  mutate(Method = pretty_method(method)) %>%
  select(Method, MAE_AROP_pp, MAE_median) %>%
  arrange(MAE_AROP_pp)

tab_mae_by_strength_pp_pretty <- tab_mae_by_strength_pp %>%
  rename(BiasVar = bias_var, Strength = bias_strength) %>%
  mutate(
    BiasVar = BiasVar,
    Strength = Strength,
    `Unweighted` = `UNWEIGHTED`,
    `IPW (grouped)` = `IPW_GROUPED_BOOL`,
    `IPW-CDF` = `IPW_CDF`,
    `DR-CDF` = `DR_CDF`
  ) %>%
  select(BiasVar, Strength, `Unweighted`, `IPW (grouped)`, `IPW-CDF`, `DR-CDF`)

# -------------------------
# Print to console
# -------------------------
msg("\nTable 6 (overall MAE, AROP in pp):")
print(tab_mae_overall_pp_pretty)

msg("\nTable 7 (MAE by bias var/strength, AROP in pp):")
print(tab_mae_by_strength_pp_pretty)

# -------------------------
# Save outputs
# -------------------------
readr::write_csv(res_all, file.path(OUT_DIR, "validation_all.csv"))
readr::write_csv(tab_mae_overall, file.path(OUT_DIR, "paper_mae_overall.csv"))
readr::write_csv(tab_mae_overall_pp, file.path(OUT_DIR, "paper_mae_overall_pp.csv"))
readr::write_csv(tab_mae_overall_pp_pretty, file.path(OUT_DIR, "paper_mae_overall_pp_pretty.csv"))
readr::write_csv(tab_mae_by_strength_pp, file.path(OUT_DIR, "paper_mae_by_biasvar_strength_pp.csv"))
readr::write_csv(tab_mae_by_strength_pp_pretty, file.path(OUT_DIR, "paper_mae_by_biasvar_strength_pp_pretty.csv"))

msg("\nSaved outputs in: ", OUT_DIR)
