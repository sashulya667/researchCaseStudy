# ============================================================
# analysis/validate_truly_biased_samples.R
# Validate methods on synthetic biased samples created by 04_create_biased_sample.
#
# Reads for each dataset:
# - <DATASET>__truth.csv
# - <DATASET>__ref_srs.csv
# - <DATASET>__np_bias_<VAR>__<strength>.csv  (mild/medium/strong)
# and compares:
# - IPW_GROUPED_BOOL
# - IPW_CDF
# - DR_CDF
# ============================================================

devtools::load_all()
rm(list = ls())

library(tidyverse)
library(here)

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

# IMPORTANT: set this to the folder where you SEE those generated CSVs
# from your screenshot.
DIR_IN <- here("outputs", "sim_biased_samples")  # change if needed

OUT_DIR <- here("outputs", "validate_truly_biased_samples")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

msg("DIR_IN = ", DIR_IN)
msg("OUT_DIR = ", OUT_DIR)

# -------------------------
# File helpers
# -------------------------
path_truth <- function(dataset_name) file.path(DIR_IN, paste0(dataset_name, "__truth.csv"))
path_ref   <- function(dataset_name) file.path(DIR_IN, paste0(dataset_name, "__ref_srs.csv"))

# Only biased nonprob microdata files (mild/medium/strong), exclude meta/truth/ref.
list_np_bias_files <- function(dataset_name) {
  files <- list.files(
    DIR_IN,
    pattern = paste0("^", dataset_name, "__np_bias_.*\\.csv$"),
    full.names = TRUE
  )
  files <- files[grepl("__(mild|medium|strong)\\.csv$", files)]
  sort(files)
}

# Parse bias_var and bias_strength from filename:
# Dataset_B__np_bias_DOU__medium.csv  -> bias_var=DOU, bias_strength=medium
parse_bias_from_filename <- function(f) {
  b <- basename(f)
  # split by "__"
  parts <- strsplit(b, "__", fixed = TRUE)[[1]]
  # parts: [1]=Dataset_B, [2]=np_bias_DOU, [3]=medium.csv
  bias_var <- sub("^np_bias_", "", parts[2])
  bias_strength <- sub("\\.csv$", "", parts[3])
  tibble(bias_var = bias_var, bias_strength = bias_strength)
}

# -------------------------
# Validation for ONE dataset
# -------------------------
validate_one_dataset <- function(dataset_name) {
  msg("\n==============================")
  msg("VALIDATE: ", dataset_name)
  msg("==============================")

  # Read truth + ref
  truth_file <- path_truth(dataset_name)
  ref_file <- path_ref(dataset_name)

  if (!file.exists(truth_file)) stop("Missing truth file: ", truth_file, call. = FALSE)
  if (!file.exists(ref_file))   stop("Missing ref file: ", ref_file, call. = FALSE)

  truth <- readr::read_csv(truth_file, show_col_types = FALSE)
  ref   <- readr::read_csv(ref_file, show_col_types = FALSE)

  # list biased nonprob
  files_np <- list_np_bias_files(dataset_name)
  if (length(files_np) == 0) stop("No biased nonprob files found for ", dataset_name, " in ", DIR_IN, call. = FALSE)

  msg("Nonprob files to validate (ONLY biased samples):")
  walk2(seq_along(files_np), files_np, ~ msg(.x, " -> ", basename(.y)))

  # run each biased sample
  out_all <- purrr::map_dfr(files_np, function(f_np) {
    nonprob <- readr::read_csv(f_np, show_col_types = FALSE)

    # hard guard: if file is wrong, fail early with a clear message
    missing_x <- setdiff(X_VARS, names(nonprob))
    if (length(missing_x) > 0) {
      stop("File does not look like biased nonprob microdata (missing X vars): ",
           basename(f_np), " missing: ", paste(missing_x, collapse = ", "),
           call. = FALSE)
    }

    # run methods (from R/runners.R)
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
    )

    info <- parse_bias_from_filename(f_np)

    # truth has: median, poverty_line, arop (from your generator)
    # Make sure truth columns exist
    stopifnot(all(c("median","poverty_line","arop") %in% names(truth)))

    est_tbl |>
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
      ) |>
      select(dataset, bias_var, bias_strength, method,
             median, poverty_line, arop,
             truth_median, truth_poverty_line, truth_arop,
             err_median, abs_err_median, err_arop, abs_err_arop,
             file)
  })

  out_all
}

# -------------------------
# RUN BOTH DATASETS
# -------------------------
res_B <- validate_one_dataset("Dataset_B")
res_E <- validate_one_dataset("Dataset_E")

res_all <- bind_rows(res_B, res_E)

msg("\nTop methods by abs_err_arop (overall):")
print(res_all |> arrange(abs_err_arop) |> select(dataset, bias_var, bias_strength, method, arop, truth_arop, abs_err_arop) |> head(20))

# Summaries: average absolute error by method
summ_arop <- res_all |>
  group_by(method) |>
  summarise(
    mean_abs_err_arop = mean(abs_err_arop, na.rm = TRUE),
    mean_abs_err_median = mean(abs_err_median, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(mean_abs_err_arop)

msg("\nSummary (mean absolute errors):")
print(summ_arop)

# Save
readr::write_csv(res_all, file.path(OUT_DIR, "validation_all.csv"))
readr::write_csv(summ_arop, file.path(OUT_DIR, "validation_summary_by_method.csv"))

msg("\nSaved:")
msg(file.path(OUT_DIR, "validation_all.csv"))
msg(file.path(OUT_DIR, "validation_summary_by_method.csv"))
