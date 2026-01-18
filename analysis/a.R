# ----------------------------
# Libraries
# ----------------------------
library(dplyr)
library(tidyr)
library(readr)

# ----------------------------
# Paths / inputs (EDIT THESE)
# ----------------------------
DATA_DIR   <- "/Users/sasha/Desktop/University/3rd semester/RCS/final/researchCaseStudy/data"
AMAP_RDATA <- "/Users/sasha/Downloads/AMAP_v0.2.1.RData"  # optional, not needed here

REF_W_COL <- "design_weight"
Y_COL     <- "EDI"
X_COLS    <- c("REG", "UEP", "DOU")

EPS_CLIP  <- 1e-6

# ----------------------------
# Load data (your loader)
# ----------------------------
load_csv <- function(fname) {
  path <- file.path(DATA_DIR, fname)
  stopifnot(file.exists(path))
  readr::read_csv(path, show_col_types = FALSE)
}

df_B_raw <- load_csv("Dataset_B.csv")
df_E_raw <- load_csv("Dataset_E.csv")
ref_raw  <- load_csv("reference_sample_probabilistic.csv")

# ----------------------------
# Helpers
# ----------------------------
assert_cols <- function(df, cols, df_name="df") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop(df_name, " is missing columns: ", paste(miss, collapse=", "))
}

make_joint_key <- function(df, x_cols) {
  df %>% mutate(.cell = paste(!!!rlang::syms(x_cols), sep = " | "))
}

preprocess_common_support <- function(np, ref, x_cols, y_col, ref_w_col) {
  assert_cols(np,  c(x_cols, y_col), "NP")
  assert_cols(ref, c(x_cols, ref_w_col), "REF")

  # 1) drop missing in X in both
  np1  <- np  %>% filter(if_all(all_of(x_cols), ~ !is.na(.)))
  ref1 <- ref %>% filter(if_all(all_of(x_cols), ~ !is.na(.)))

  # 2) common support: keep only NP cells where REF has positive total design weight
  ref_cells <- ref1 %>%
    make_joint_key(x_cols) %>%
    group_by(.cell) %>%
    summarise(ref_w_sum = sum(.data[[ref_w_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(ref_w_sum > 0)

  np2 <- np1 %>%
    make_joint_key(x_cols) %>%
    inner_join(ref_cells %>% select(.cell), by = ".cell") %>%
    select(-.cell)

  list(np = np2, ref = ref1)
}

fit_ipw_row_level <- function(np, ref, x_cols, ref_w_col, eps_clip = 1e-6) {
  # Force consistent factor levels across stacked data
  stacked <- bind_rows(
    np  %>% mutate(.S = 1, .w = 1),
    ref %>% mutate(.S = 0, .w = .data[[ref_w_col]])
  )

  # Convert X to factors (VERY important for categorical cells)
  for (xc in x_cols) stacked[[xc]] <- as.factor(stacked[[xc]])

  # Row-level weighted logit (the "classical" fit that can collapse)
  fml <- as.formula(paste0(".S ~ ", paste(paste0("factor(", x_cols, ")"), collapse = " + ")))
  mod <- glm(fml, data = stacked, family = binomial(), weights = .w)

  # Predict propensities on NP rows
  np2 <- np
  for (xc in x_cols) np2[[xc]] <- factor(np2[[xc]], levels = levels(stacked[[xc]]))

  p_hat <- as.numeric(predict(mod, newdata = np2, type = "response"))
  eta   <- as.numeric(predict(mod, newdata = np2, type = "link"))

  # Diagnostics: propensity collapse check
  diag <- list(
    p_summary   = summary(p_hat),
    p_quantile  = quantile(p_hat, probs = c(0, .001, .01, .1, .5, .9, .99, .999, 1), na.rm = TRUE),
    eta_summary = summary(eta)
  )

  # Weights: r(x)=(1-p)/p with clipping
  p_clip <- pmin(pmax(p_hat, eps_clip), 1 - eps_clip)
  r      <- (1 - p_clip) / p_clip

  # Normalize to sum to 1 (as in your paper)
  w_norm <- r / sum(r)

  # Weight diagnostics + ESS
  ess <- (sum(w_norm)^2) / sum(w_norm^2)
  w_diag <- list(
    r_summary = summary(r),
    w_summary = summary(w_norm),
    w_sd      = sd(w_norm),
    ess       = ess
  )

  list(model = mod, p_hat = p_hat, eta = eta, r = r, w = w_norm, diag = diag, w_diag = w_diag)
}

fit_ipw_grouped_cell <- function(np, ref, x_cols, ref_w_col, eps_clip = 1e-6) {
  # Build joint cells
  np_c <- np %>%
    make_joint_key(x_cols) %>%
    count(.cell, name = "y1")

  ref_c <- ref %>%
    make_joint_key(x_cols) %>%
    group_by(.cell) %>%
    summarise(y0 = sum(.data[[ref_w_col]], na.rm = TRUE), .groups = "drop")

  cells <- full_join(np_c, ref_c, by = ".cell") %>%
    mutate(y1 = replace_na(y1, 0),
           y0 = replace_na(y0, 0)) %>%
    # keep supported cells (y0>0 if you already enforced common support; still safe)
    filter(y0 > 0)

  # Recover X columns by splitting .cell back (keeps it simple)
  cells <- cells %>%
    separate(.cell, into = x_cols, sep = " \\| ", remove = FALSE)

  for (xc in x_cols) cells[[xc]] <- as.factor(cells[[xc]])

  fml <- as.formula(paste0("cbind(y1, y0) ~ ", paste(paste0("factor(", x_cols, ")"), collapse = " + ")))
  mod <- glm(fml, data = cells, family = binomial())

  # Map fitted p(cell) back to NP rows
  p_cell <- cells %>%
    mutate(p_hat = as.numeric(predict(mod, type = "response"))) %>%
    select(.cell, p_hat)

  np2 <- np %>%
    make_joint_key(x_cols) %>%
    left_join(p_cell, by = ".cell")

  # weights for NP
  p_clip <- pmin(pmax(np2$p_hat, eps_clip), 1 - eps_clip)
  r      <- (1 - p_clip) / p_clip
  w_norm <- r / sum(r)
  ess <- (sum(w_norm)^2) / sum(w_norm^2)

  list(model = mod,
       cell_table = cells,
       p_hat = np2$p_hat,
       r = r,
       w = w_norm,
       diag = list(p_summary = summary(np2$p_hat),
                   p_quantile = quantile(np2$p_hat, probs = c(0,.01,.1,.5,.9,.99,1), na.rm=TRUE)),
       w_diag = list(r_summary = summary(r),
                     w_summary = summary(w_norm),
                     w_sd = sd(w_norm),
                     ess = ess))
}

# ----------------------------
# RUN: choose dataset
# ----------------------------
run_one <- function(np_raw, ref_raw, name = "Dataset_B") {
  cat("\n============================\n")
  cat("RUN:", name, "\n")
  cat("============================\n")

  pp <- preprocess_common_support(np_raw, ref_raw, X_COLS, Y_COL, REF_W_COL)
  np  <- pp$np
  ref <- pp$ref

  cat("NP n:", nrow(np), " | REF n:", nrow(ref), "\n")

  # Row-level IPW (expected to collapse in your setting)
  cat("\n--- Row-level IPW fit ---\n")
  row <- fit_ipw_row_level(np, ref, X_COLS, REF_W_COL, eps_clip = EPS_CLIP)

  print(row$diag$p_summary)
  print(row$diag$p_quantile)
  print(row$diag$eta_summary)

  cat("\nWeights after clipping/normalization:\n")
  print(row$w_diag$w_summary)
  cat("w sd:", row$w_diag$w_sd, "\n")
  cat("ESS:", row$w_diag$ess, "\n")

  # Quick “collapse flag”
  p_unique <- length(unique(round(row$p_hat, 16)))
  cat("Unique p_hat (rounded 1e-16):", p_unique, "\n")

  # Grouped-cell IPW (stable)
  cat("\n--- Grouped-cell IPW fit ---\n")
  grp <- fit_ipw_grouped_cell(np, ref, X_COLS, REF_W_COL, eps_clip = EPS_CLIP)

  print(grp$diag$p_summary)
  print(grp$diag$p_quantile)

  cat("\nGrouped weights:\n")
  print(grp$w_diag$w_summary)
  cat("w sd:", grp$w_diag$w_sd, "\n")
  cat("ESS:", grp$w_diag$ess, "\n")

  invisible(list(row = row, grp = grp))
}

# ----------------------------
# Execute
# ----------------------------
res_B <- run_one(df_B_raw, ref_raw, "Dataset_B")
# res_E <- run_one(df_E_raw, ref_raw, "Dataset_E")
