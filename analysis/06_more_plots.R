# ============================================================
# analysis/01c_paper_artifacts.R
# One script that produces ALL paper-ready tables + figures.
# Outputs go to: outputs/paper_figures/
# ============================================================

devtools::load_all()
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
  library(readr)
  library(dplyr)
  library(tidyr)
})

HAS_PATCHWORK   <- requireNamespace("patchwork", quietly = TRUE)
HAS_KABLE_EXTRA <- requireNamespace("kableExtra", quietly = TRUE)
HAS_KNITR       <- requireNamespace("knitr", quietly = TRUE)

# -------------------------
# Settings
# -------------------------
DATA_DIR <- here("data")
OUT_DIR  <- here("outputs", "paper_figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SIM_DIR <- here("outputs", "sim_biased_samples")  # from analysis/04_create_biased_sample.R

INCOME_VAR <- "EDI"
X_VARS <- c("REG","UEP","DOU")
REF_WEIGHT <- "design_weight"

# These should match analysis/04_create_biased_sample.R
N_REF  <- 10013
N_NP   <- 10000
MAX_PI <- 0.95
BIAS_LEVELS <- c(mild = 0.6, medium = 1.2, strong = 2.0)

# -------------------------
# Helpers
# -------------------------
stop_if_missing_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

load_csv <- function(fname) {
  path <- file.path(DATA_DIR, fname)
  stop_if_missing_file(path)
  readr::read_csv(path, show_col_types = FALSE)
}

write_table_tex <- function(df, path, caption) {
  if (HAS_KABLE_EXTRA) {
    kableExtra::kbl(df, format = "latex", booktabs = TRUE, digits = 4, caption = caption) |>
      kableExtra::kable_styling(latex_options = c("hold_position")) |>
      writeLines(path)
  } else {
    if (!HAS_KNITR) stop("Need knitr for LaTeX table output. Install: install.packages('knitr')", call. = FALSE)
    cat("% Install kableExtra for nicer LaTeX output.\n", file = path)
    cat(knitr::kable(df, format = "latex", booktabs = TRUE, digits = 4, caption = caption),
        file = path, append = TRUE)
  }
}

order_levels_safely <- function(level_chr) {
  x <- unique(level_chr)
  x_num <- suppressWarnings(as.numeric(x))
  if (all(!is.na(x_num))) as.character(sort(unique(x_num))) else x
}

signed_log <- function(x) sign(x) * log1p(abs(x))

# -------------------------
# Load + harmonize + drop missing X
# -------------------------
df_B_raw <- load_csv("Dataset_B.csv")
df_E_raw <- load_csv("Dataset_E.csv")
ref_raw  <- load_csv("reference_sample_probabilistic.csv")

make_levels <- function(...) sort(unique(unlist(list(...))))
lev_REG <- make_levels(df_B_raw$REG, df_E_raw$REG, ref_raw$REG)
lev_UEP <- make_levels(df_B_raw$UEP, df_E_raw$UEP, ref_raw$UEP)
lev_DOU <- make_levels(df_B_raw$DOU, df_E_raw$DOU, ref_raw$DOU)

harmonize <- function(df) {
  df |>
    mutate(
      REG = factor(REG, levels = lev_REG),
      UEP = factor(UEP, levels = lev_UEP),
      DOU = factor(DOU, levels = lev_DOU)
    )
}

drop_missing_x <- function(df, x_vars = X_VARS) df |> filter(!if_any(all_of(x_vars), is.na))

df_Bx <- drop_missing_x(harmonize(df_B_raw))
df_Ex <- drop_missing_x(harmonize(df_E_raw))
refx  <- drop_missing_x(harmonize(ref_raw))

# Quick sanity print (optional)
cat("\nDropped rows due to missing X:\n")
cat("Dataset_B:", nrow(df_B_raw) - nrow(df_Bx), "\n")
cat("Dataset_E:", nrow(df_E_raw) - nrow(df_Ex), "\n")
cat("Reference :", nrow(ref_raw)  - nrow(refx),  "\n\n")

# ============================================================
# TABLE 1: Data overview
# ============================================================
income_summary <- function(y) {
  y <- y[is.finite(y) & !is.na(y)]
  tibble(
    income_median = median(y),
    income_p25    = unname(quantile(y, 0.25)),
    income_p75    = unname(quantile(y, 0.75)),
    income_iqr    = IQR(y)
  )
}

missing_summary <- function(df, vars) {
  tibble(var = vars) |>
    mutate(
      n_missing = purrr::map_int(var, ~sum(is.na(df[[.x]]))),
      p_missing = n_missing / nrow(df)
    )
}

weight_summary <- function(df, w) {
  if (!(w %in% names(df))) return(tibble(weight_mass = NA, w_min = NA, w_median = NA, w_max = NA))
  x <- df[[w]]
  x <- x[is.finite(x) & !is.na(x)]
  tibble(
    weight_mass = sum(x),
    w_min = min(x),
    w_median = median(x),
    w_max = max(x)
  )
}

source_overview <- function(df, name, weight_var = NULL) {
  miss <- missing_summary(df, c(INCOME_VAR, X_VARS)) |>
    pivot_wider(names_from = var, values_from = c(n_missing, p_missing))

  levs <- tibble(
    n_levels_REG = nlevels(df$REG),
    n_levels_UEP = nlevels(df$UEP),
    n_levels_DOU = nlevels(df$DOU)
  )

  inc <- income_summary(df[[INCOME_VAR]])

  wsum <- if (!is.null(weight_var)) weight_summary(df, weight_var) else
    tibble(weight_mass = NA, w_min = NA, w_median = NA, w_max = NA)

  tibble(source = name, n_rows = nrow(df)) |>
    bind_cols(wsum) |>
    bind_cols(levs) |>
    bind_cols(inc) |>
    bind_cols(miss)
}

tab_data_overview <- bind_rows(
  source_overview(df_Bx, "Dataset_B (nonprob)"),
  source_overview(df_Ex, "Dataset_E (nonprob)"),
  source_overview(refx,  "Reference (prob)", weight_var = REF_WEIGHT)
)

write_csv(tab_data_overview, file.path(OUT_DIR, "table_data_overview.csv"))
write_table_tex(tab_data_overview, file.path(OUT_DIR, "table_data_overview.tex"),
                "Data overview (nonprob samples and probability reference sample).")

# ============================================================
# TABLE 2: Simulation design + truth
# ============================================================
truth_arop <- function(pop_df) {
  y <- pop_df[[INCOME_VAR]]
  y <- y[is.finite(y) & !is.na(y)]
  med <- median(y)
  z <- 0.6 * med
  tibble(median = med, poverty_line = z, arop = mean(y < z))
}

tab_sim_design <- bind_rows(
  tibble(dataset = "Dataset_B", N_pseudo_pop = nrow(df_Bx),
         n_ref = N_REF, n_np_target = N_NP, pi_max = MAX_PI,
         bias_mild = unname(BIAS_LEVELS["mild"]),
         bias_medium = unname(BIAS_LEVELS["medium"]),
         bias_strong = unname(BIAS_LEVELS["strong"])) |>
    bind_cols(truth_arop(df_Bx)),
  tibble(dataset = "Dataset_E", N_pseudo_pop = nrow(df_Ex),
         n_ref = N_REF, n_np_target = N_NP, pi_max = MAX_PI,
         bias_mild = unname(BIAS_LEVELS["mild"]),
         bias_medium = unname(BIAS_LEVELS["medium"]),
         bias_strong = unname(BIAS_LEVELS["strong"])) |>
    bind_cols(truth_arop(df_Ex))
)

write_csv(tab_sim_design, file.path(OUT_DIR, "table_sim_design.csv"))
write_table_tex(tab_sim_design, file.path(OUT_DIR, "table_sim_design.tex"),
                "Pseudo-population validation design and ground truth (per dataset).")

# ============================================================
# TABLE 3: Realized biased sample sizes (summary)
# ============================================================
bias_levels <- c("mild","medium","strong")

read_sim_n <- function(dataset_name, bias_var, bias_level) {
  f <- file.path(SIM_DIR, paste0(dataset_name, "__np_bias_", bias_var, "__", bias_level, ".csv"))
  if (!file.exists(f)) return(NULL)
  df <- readr::read_csv(f, show_col_types = FALSE)
  tibble(dataset = dataset_name, bias_var = bias_var, bias_level = bias_level, n_realized = nrow(df))
}

tab_sim_sizes <- bind_rows(
  purrr::map_dfr(X_VARS, \(v) purrr::map_dfr(bias_levels, \(bl) read_sim_n("Dataset_B", v, bl))),
  purrr::map_dfr(X_VARS, \(v) purrr::map_dfr(bias_levels, \(bl) read_sim_n("Dataset_E", v, bl)))
)

if (is.null(tab_sim_sizes) || nrow(tab_sim_sizes) == 0) {
  message("No saved simulation samples found in: ", SIM_DIR)
} else {
  tab_sim_sizes_summary <- tab_sim_sizes |>
    group_by(dataset) |>
    summarise(min_n = min(n_realized), median_n = median(n_realized), max_n = max(n_realized), .groups = "drop")

  write_csv(tab_sim_sizes, file.path(OUT_DIR, "table_sim_realized_sizes_full.csv"))
  write_csv(tab_sim_sizes_summary, file.path(OUT_DIR, "table_sim_realized_sizes_summary.csv"))

  write_table_tex(tab_sim_sizes_summary, file.path(OUT_DIR, "table_sim_realized_sizes_summary.tex"),
                  "Realized sizes of biased nonprobability samples (summary across scenarios).")
}

# ============================================================
# FIG A: Marginals (B vs E vs Reference-weighted)
# ============================================================
prop_1d_unweighted <- function(df, var) {
  df |> count(.data[[var]], name = "n") |> mutate(prop = n / sum(n)) |>
    transmute(level = as.character(.data[[var]]), prop)
}
prop_1d_weighted <- function(df, var, w) {
  df |> group_by(.data[[var]]) |> summarise(w_sum = sum(.data[[w]], na.rm = TRUE), .groups="drop") |>
    mutate(prop = w_sum / sum(w_sum)) |> transmute(level = as.character(.data[[var]]), prop)
}
make_marginal_long <- function(df, source_name, is_ref = FALSE) {
  map_dfr(X_VARS, \(v) {
    p <- if (is_ref) prop_1d_weighted(df, v, REF_WEIGHT) else prop_1d_unweighted(df, v)
    p |> mutate(var = v, source = source_name)
  })
}

marg_long <- bind_rows(
  make_marginal_long(df_Bx, "Dataset_B"),
  make_marginal_long(df_Ex, "Dataset_E"),
  make_marginal_long(refx,  "Reference (design-weighted)", is_ref = TRUE)
) |>
  group_by(var) |>
  mutate(level = factor(level, levels = order_levels_safely(level))) |>
  ungroup() |>
  mutate(var = factor(var, levels = c("DOU","REG","UEP")))

p_marg <- ggplot(marg_long, aes(x = level, y = prop, fill = source)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  facet_wrap(~ var, scales = "free_x", ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x=NULL, y="Proportion", fill=NULL,
       title="Marginal composition of auxiliary variables",
       subtitle="Dataset_B vs Dataset_E vs Reference") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

ggsave(file.path(OUT_DIR, "fig_marginals_B_E_ref.png"), p_marg, width=9, height=9, dpi=300)

# ============================================================
# FIG B: Income distributions (signed-log hist, B vs E)
# ============================================================
inc_long <- bind_rows(
  df_Bx |> transmute(dataset="Dataset_B", slog = signed_log(.data[[INCOME_VAR]])),
  df_Ex |> transmute(dataset="Dataset_E", slog = signed_log(.data[[INCOME_VAR]]))
)

p_inc <- ggplot(inc_long, aes(x = slog)) +
  geom_histogram(bins = 140) +
  facet_wrap(~ dataset, ncol = 2, scales = "free_y") +
  labs(title = "Income distributions (signed-log scale)",
       subtitle = "Dataset_B vs Dataset_E",
       x = "signed-log(EDI)", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "fig_income_hist_signedlog_B_E.png"), p_inc, width=12, height=5, dpi=300)

cat("\nDONE. Paper artifacts written to:\n", OUT_DIR, "\n")

###################################









# ============================================================
# analysis/05_make_sim_biasdiff_figures.R
# Best visual simulation sanity-check plot:
#   bars = (biased sample proportion - pseudo-population proportion)
#   faceted by (bias variable) x (bias strength)
# Reads biased samples written by analysis/04_create_biased_sample.R
# ============================================================

devtools::load_all()
rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(scales)
})

HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)

# -------------------------
# Settings (match your pipeline)
# -------------------------
DATA_DIR <- here("data")
SIM_DIR  <- here("outputs", "sim_biased_samples")

# Where to save figures for Overleaf: copy these into ./figures/
OUT_DIR  <- here("outputs", "paper_figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

X_VARS <- c("REG", "UEP", "DOU")
BIAS_LEVELS <- c("mild", "medium", "strong")

# -------------------------
# Helpers
# -------------------------
stop_if_missing_file <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}

load_csv <- function(fname) {
  path <- file.path(DATA_DIR, fname)
  stop_if_missing_file(path)
  readr::read_csv(path, show_col_types = FALSE)
}

drop_missing_x <- function(df, x_vars = X_VARS) {
  df |> filter(!if_any(all_of(x_vars), is.na))
}

# stable ordering of levels on x-axis (numeric if possible)
order_levels_safely <- function(level_chr) {
  lev <- unique(level_chr)
  lev_num <- suppressWarnings(as.numeric(lev))
  if (all(!is.na(lev_num))) as.character(sort(unique(lev_num))) else lev
}

read_sim_sample <- function(dataset_name, bias_var, bias_level) {
  f <- file.path(SIM_DIR, paste0(dataset_name, "__np_bias_", bias_var, "__", bias_level, ".csv"))
  if (!file.exists(f)) return(NULL)
  readr::read_csv(f, show_col_types = FALSE)
}

props_1d <- function(df, var) {
  df |>
    count(.data[[var]], name = "n") |>
    mutate(prop = n / sum(n)) |>
    transmute(level = as.character(.data[[var]]), prop)
}

make_biasdiff_long <- function(pop_df, dataset_name) {
  # pop proportions for each var
  out <- map_dfr(X_VARS, function(v) {
    p_pop <- props_1d(pop_df, v) |> rename(pop_prop = prop)

    map_dfr(BIAS_LEVELS, function(bl) {
      samp <- read_sim_sample(dataset_name, v, bl)
      if (is.null(samp)) return(tibble())

      p_s <- props_1d(samp, v) |> rename(samp_prop = prop)

      full_join(p_pop, p_s, by = "level") |>
        mutate(
          pop_prop  = replace_na(pop_prop, 0),
          samp_prop = replace_na(samp_prop, 0),
          diff      = samp_prop - pop_prop,
          dataset   = dataset_name,
          bias_var  = v,
          bias_level = bl
        )
    })
  })

  out |>
    mutate(
      bias_level = factor(bias_level, levels = BIAS_LEVELS),
      bias_var   = factor(bias_var, levels = c("REG", "UEP", "DOU")),
      sign = ifelse(diff >= 0, "Over-represented", "Under-represented")
    ) |>
    group_by(bias_var) |>
    mutate(level = factor(level, levels = order_levels_safely(level))) |>
    ungroup()
}

plot_biasdiff <- function(biasdiff_df, dataset_name) {
  ggplot(biasdiff_df, aes(x = level, y = diff, fill = sign)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_col(width = 0.78) +
    facet_grid(bias_var ~ bias_level, scales = "free_x") +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = NULL,
      y = expression(Delta~"proportion (biased - population)"),
      fill = NULL,
      title = paste0("Simulation sanity check: induced marginal bias (", dataset_name, ")"),
      subtitle = "Bars show how the biased nonprob sample deviates from the pseudo-population for the biased variable"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "top"
    )
}

# -------------------------
# Load pseudo-populations (same as your simulation: drop missing X)
# -------------------------
df_B <- load_csv("Dataset_B.csv") |> drop_missing_x()
df_E <- load_csv("Dataset_E.csv") |> drop_missing_x()

# -------------------------
# Build and plot
# -------------------------
bd_B <- make_biasdiff_long(df_B, "Dataset_B")
bd_E <- make_biasdiff_long(df_E, "Dataset_E")

if (nrow(bd_B) == 0) stop("No simulation files found for Dataset_B in: ", SIM_DIR)
if (nrow(bd_E) == 0) stop("No simulation files found for Dataset_E in: ", SIM_DIR)

pB <- plot_biasdiff(bd_B, "Dataset_B")
pE <- plot_biasdiff(bd_E, "Dataset_E")

ggsave(file.path(OUT_DIR, "fig_sim_biasdiff_B.png"), pB, width = 12, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "fig_sim_biasdiff_E.png"), pE, width = 12, height = 7, dpi = 300)

if (HAS_PATCHWORK) {
  pBE <- pB / pE + patchwork::plot_layout(heights = c(1, 1))
  ggsave(file.path(OUT_DIR, "fig_sim_biasdiff_B_and_E.png"), pBE, width = 12, height = 14, dpi = 300)
}

cat("\nSaved:\n",
    file.path(OUT_DIR, "fig_sim_biasdiff_B.png"), "\n",
    file.path(OUT_DIR, "fig_sim_biasdiff_E.png"), "\n")
