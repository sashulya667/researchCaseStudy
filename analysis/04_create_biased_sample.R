# ============================================================
# analysis/04_create_biased_sample.R
# Simulation: treat dataset as pseudo-population, generate
# (1) probability reference sample (SRSWOR) with design weights
# (2) biased "nonprob" samples biased on REG or UEP or DOU
#     with 3 bias strengths
# ============================================================

devtools::load_all()
rm(list = ls())

library(tidyverse)
library(here)

msg <- function(...) cat(paste0(..., "\n"))

# -------------------------
# Settings
# -------------------------
X_VARS <- c("REG", "UEP", "DOU")
INCOME_VAR <- "EDI"
POVERTY_FRAC <- 0.6

# sizes: reference sample and "nonprob" biased sample
N_REF  <- 10013
N_NP   <- 10000        # expected size for biased nonprob sample (Poisson -> approx)
MAX_PI <- 0.95         # cap inclusion probs for Poisson sampling

# bias strengths (bigger => stronger tilt away from population)
BIAS_LEVELS <- c(mild = 0.6, medium = 1.2, strong = 2.0)

OUT_DIR <- here("outputs", "sim_biased_samples")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

# -------------------------
# Helpers
# -------------------------

# compute "truth" AROP on full pseudo-population
truth_arop <- function(df, income_var = INCOME_VAR, poverty_frac = POVERTY_FRAC) {
  y <- df[[income_var]]
  y <- y[is.finite(y) & !is.na(y)]
  med <- median(y)
  z <- poverty_frac * med
  arop <- mean(y < z)
  tibble(median = med, poverty_line = z, arop = arop)
}

# SRSWOR reference sample with known design weights
make_reference_srs <- function(pop_df, n_ref = N_REF, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(pop_df)
  idx <- sample.int(N, size = n_ref, replace = FALSE)
  ref <- pop_df[idx, , drop = FALSE] |> mutate(design_weight = N / n_ref)
  ref
}

# create level tilt scores centered around 0 (so strength=0 gives no bias)
# for k levels: scores are equally spaced from -1 to +1
level_scores <- function(f) {
  lev <- levels(f)
  k <- length(lev)
  s <- seq(-1, 1, length.out = k)
  names(s) <- lev
  s
}

# find scale c for Poisson sampling with pi_i = min(MAX_PI, c * w_i),
# such that sum(pi_i) ~= n_target
solve_scale_for_poisson <- function(w, n_target, max_pi = MAX_PI) {
  w <- as.numeric(w)
  stopifnot(all(w >= 0), is.finite(n_target), n_target > 0)

  f <- function(c) sum(pmin(max_pi, c * w)) - n_target

  # bracket c
  c_lo <- 0
  c_hi <- 1
  while (f(c_hi) < 0) {
    c_hi <- c_hi * 2
    if (c_hi > 1e9) break
  }
  uniroot(f, interval = c(c_lo, c_hi), tol = 1e-10)$root
}

# biased Poisson sample: bias on one categorical var only
# strength controls tilt: pi higher for high-score levels, lower for low-score
make_biased_nonprob <- function(pop_df, var, n_target = N_NP, strength = 1.0,
                                max_pi = MAX_PI, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  f <- pop_df[[var]]
  stopifnot(is.factor(f))

  s <- level_scores(f)
  score_i <- unname(s[as.character(f)])

  # tilt weights
  w_i <- exp(strength * score_i)

  # scale to hit expected size
  c_scale <- solve_scale_for_poisson(w_i, n_target = n_target, max_pi = max_pi)
  pi_i <- pmin(max_pi, c_scale * w_i)

  take <- rbinom(nrow(pop_df), size = 1, prob = pi_i) == 1
  samp <- pop_df[take, , drop = FALSE] |>
    mutate(
      pi_true = pi_i[take],
      bias_var = var,
      bias_strength = strength
    )

  list(sample = samp, pi = pi_i, c_scale = c_scale, n = nrow(samp))
}

# marginal proportions helper (unweighted)
prop_1d <- function(df, var) {
  df |> count(.data[[var]], name = "n") |> mutate(prop = n / sum(n)) |> select(!!var, prop)
}

max_abs_diff_marginal <- function(pop_df, samp_df, var) {
  p_pop <- prop_1d(pop_df, var) |> rename(pop = prop)
  p_smp <- prop_1d(samp_df, var) |> rename(smp = prop)
  full_join(p_pop, p_smp, by = var) |>
    mutate(across(c(pop, smp), ~replace_na(.x, 0)),
           diff = smp - pop) |>
    summarise(max_abs_diff = max(abs(diff))) |>
    pull(max_abs_diff)
}

# run full generation for one pseudo-population dataset
generate_for_dataset <- function(pop_df, dataset_name) {
  msg("\n==============================")
  msg("DATASET: ", dataset_name, " (pseudo-population)")
  msg("==============================")

  # truth
  tru <- truth_arop(pop_df) |> mutate(dataset = dataset_name) |> select(dataset, everything())
  msg("Truth AROP (population):"); print(tru)

  # reference sample
  ref <- make_reference_srs(pop_df, n_ref = N_REF, seed = 100)
  msg("Reference SRS size: ", nrow(ref), " | design_weight = N/n =", unique(ref$design_weight))

  # save reference
  write_csv(ref, file.path(OUT_DIR, paste0(dataset_name, "__ref_srs.csv")))
  write_csv(tru, file.path(OUT_DIR, paste0(dataset_name, "__truth.csv")))

  # biased nonprob samples
  all_meta <- list()
  for (bias_var in X_VARS) {
    for (nm in names(BIAS_LEVELS)) {
      strength <- unname(BIAS_LEVELS[[nm]])
      seed <- 1000 + as.integer(factor(bias_var)) * 100 + which(names(BIAS_LEVELS) == nm)

      obj <- make_biased_nonprob(
        pop_df = pop_df,
        var = bias_var,
        n_target = N_NP,
        strength = strength,
        max_pi = MAX_PI,
        seed = seed
      )

      samp <- obj$sample
      mad_bias <- max_abs_diff_marginal(pop_df, samp, bias_var)

      meta <- tibble(
        dataset = dataset_name,
        type = "nonprob_biased",
        bias_var = bias_var,
        bias_level = nm,
        bias_strength = strength,
        n_target = N_NP,
        n_realized = obj$n,
        max_pi = MAX_PI,
        c_scale = obj$c_scale,
        max_abs_diff_on_bias_var = mad_bias
      )
      all_meta[[paste(bias_var, nm, sep = "_")]] <- meta

      # save sample
      out_name <- paste0(dataset_name, "__np_bias_", bias_var, "__", nm, ".csv")
      write_csv(samp, file.path(OUT_DIR, out_name))

      msg("Saved ", out_name,
          " | n=", obj$n,
          " | bias_var=", bias_var,
          " | level=", nm,
          " | max_abs_diff(", bias_var, ")=", round(mad_bias, 4))
    }
  }

  meta_tbl <- bind_rows(all_meta)
  write_csv(meta_tbl, file.path(OUT_DIR, paste0(dataset_name, "__np_bias_meta.csv")))

  list(truth = tru, ref = ref, meta = meta_tbl)
}

# -------------------------
# Load data
# -------------------------
load_data <- function(filename) read.csv(file.path(here("data"), filename))

df_B_raw <- load_data("Dataset_B.csv")
df_E_raw <- load_data("Dataset_E.csv")

msg("Loaded rows: B=", nrow(df_B_raw), " | E=", nrow(df_E_raw))

# -------------------------
# Harmonize factors + drop missing X
# -------------------------
# Use union of levels between B and E (no ref needed here)
lev_REG <- sort(unique(c(df_B_raw$REG, df_E_raw$REG)))
lev_UEP <- sort(unique(c(df_B_raw$UEP, df_E_raw$UEP)))
lev_DOU <- sort(unique(c(df_B_raw$DOU, df_E_raw$DOU)))

prep_pop <- function(df) {
  df |>
    mutate(
      REG = factor(REG, levels = lev_REG),
      UEP = factor(UEP, levels = lev_UEP),
      DOU = factor(DOU, levels = lev_DOU)
    ) |>
    filter(!if_any(all_of(X_VARS), is.na))
}

df_B <- prep_pop(df_B_raw)
df_E <- prep_pop(df_E_raw)

msg("After drop_missing_x: B=", nrow(df_B), " | E=", nrow(df_E))

# -------------------------
# Generate simulation samples
# -------------------------
res_B <- generate_for_dataset(df_B, "Dataset_B")
res_E <- generate_for_dataset(df_E, "Dataset_E")

msg("\nDONE. Outputs in: ", OUT_DIR)
