# ----------------------------
# Libraries
# ----------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(readr)
library(scales)

# ----------------------------
# Paths / inputs (EDIT THESE)
# ----------------------------
DATA_DIR   <- "/Users/sasha/Desktop/University/3rd semester/RCS/final/researchCaseStudy/data"
AMAP_RDATA <- "/Users/sasha/Downloads/AMAP_v0.2.1.RData"

REF_W_COL <- "design_weight"  # <- from your ref file
Y_COL     <- "EDI"

# ----------------------------
# Load AMELIA map
# ----------------------------
load(AMAP_RDATA)  # loads object AMAP
AMAP_sf <- st_as_sf(AMAP) %>% select(REG, geometry)

# ----------------------------
# Load data
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
# Preprocessing (same as before)
# ----------------------------
AUX_VARS <- c("REG", "UEP", "DOU")

drop_missing_aux <- function(df) df %>% filter(!if_any(all_of(AUX_VARS), is.na))

get_ref_support_cells <- function(ref_df) {
  ref_df %>%
    group_by(across(all_of(AUX_VARS))) %>%
    summarise(ref_mass = sum(.data[[REF_W_COL]], na.rm = TRUE), .groups = "drop") %>%
    filter(ref_mass > 0) %>%
    select(all_of(AUX_VARS))
}

enforce_common_support <- function(np_df, ref_support_cells) {
  np_df %>% semi_join(ref_support_cells, by = AUX_VARS)
}

df_B <- drop_missing_aux(df_B_raw)
df_E <- drop_missing_aux(df_E_raw)
ref  <- drop_missing_aux(ref_raw)

ref_support_cells <- get_ref_support_cells(ref)
df_B_cs <- enforce_common_support(df_B, ref_support_cells)
df_E_cs <- enforce_common_support(df_E, ref_support_cells)

# ----------------------------
# FIGURE 1: REG map (clean, official-style)
# ----------------------------
p_map_reg <- ggplot(AMAP_sf) +
  geom_sf(aes(fill = factor(REG)), color = "white", linewidth = 0.15) +
  labs(title = "AMELIA geography: regions (REG)", fill = "REG") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

# ----------------------------
# FIGURE 2: Auxiliary imbalance vs reference
#   - Output in PERCENTAGE POINTS (pp): 0.15 -> 15 pp
#   - Axis shows integers like -15, -10, 0, 10, 15
# ----------------------------
marginal_share_np <- function(df, var, dataset_label) {
  df %>%
    count(var_level = .data[[var]], name = "mass") %>%
    mutate(
      prop = mass / sum(mass),
      var = var,
      dataset = dataset_label
    ) %>%
    select(var, var_level, dataset, prop)
}

marginal_share_ref <- function(ref_df, var) {
  ref_df %>%
    group_by(var_level = .data[[var]]) %>%
    summarise(mass = sum(.data[[REF_W_COL]], na.rm = TRUE), .groups = "drop") %>%
    mutate(
      prop = mass / sum(mass),
      var = var
    ) %>%
    select(var, var_level, prop)
}

# Build imbalance table: (NP share - REF share) * 100 = percentage points
imbalance_one_dataset <- function(np_df, dataset_label) {
  out <- lapply(AUX_VARS, function(v) {
    np_p  <- marginal_share_np(np_df, v, dataset_label)
    ref_p <- marginal_share_ref(ref, v)

    np_p %>%
      left_join(ref_p, by = c("var", "var_level"), suffix = c("_np", "_ref")) %>%
      mutate(diff_pp = 100 * (prop_np - prop_ref)) %>%     # <-- percentage points
      select(var, var_level, dataset, diff_pp)
  })
  bind_rows(out)
}

imb_df <- bind_rows(
  imbalance_one_dataset(df_B_cs, "Dataset_B"),
  imbalance_one_dataset(df_E_cs, "Dataset_E")
)

p_imbalance <- ggplot(imb_df, aes(x = factor(var_level), y = diff_pp, fill = dataset)) +
  geom_col(width = 0.85, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ var, nrow = 1, scales = "free_x") +
  scale_y_continuous(labels = function(x) paste0(x, " pp")) +
  labs(
    title = "Auxiliary-variable imbalance relative to the reference sample",
    subtitle = "Bars show (nonprobability share − reference share) in percentage points",
    x = "Level",
    y = "Difference",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

# ----------------------------
# Print & Save
# ----------------------------
print(p_map_reg)
print(p_imbalance)

ggsave("fig_D1_map_reg.png",           p_map_reg,   width = 12, height = 4.5, dpi = 300)
ggsave("fig_D2_imbalance_aux_pp.png",  p_imbalance, width = 12, height = 4.8, dpi = 300)








# -----------------------------------------------------



# ============================
# Libraries
# ============================
library(dplyr)
library(ggplot2)
library(sf)
library(readr)
library(tidyr)
library(patchwork)
library(scales)

# ============================
# Paths (EDIT)
# ============================
DATA_DIR   <- "/Users/sasha/Desktop/University/3rd semester/RCS/final/researchCaseStudy/data"
AMAP_RDATA <- "/Users/sasha/Downloads/AMAP_v0.2.1.RData"

REF_W_COL <- "design_weight"
Y_COL     <- "EDI"
AUX_VARS  <- c("REG", "UEP", "DOU")

# ============================
# Load map
# ============================
load(AMAP_RDATA)     # loads AMAP
AMAP_sf <- st_as_sf(AMAP) %>% select(REG, geometry)

# ============================
# Load data
# ============================
load_csv <- function(f) read_csv(file.path(DATA_DIR, f), show_col_types = FALSE)

df_B <- load_csv("Dataset_B.csv")
df_E <- load_csv("Dataset_E.csv")
ref  <- load_csv("reference_sample_probabilistic.csv")

# ============================
# Helper: REG coverage per sample
# ============================
reg_presence <- function(df, label) {
  df %>%
    distinct(REG) %>%
    mutate(present = label)
}

reg_map_df <- bind_rows(
  reg_presence(ref,  "REF"),
  reg_presence(df_B, "Dataset_B"),
  reg_presence(df_E, "Dataset_E")
)

AMAP_reg <- AMAP_sf %>%
  left_join(reg_map_df, by = "REG")

# ============================
# FIGURE 1 — REG maps side by side
# ============================
p_maps <- ggplot(AMAP_reg) +
  geom_sf(aes(fill = present), color = "white", linewidth = 0.15) +
  facet_wrap(~ present, nrow = 1) +
  scale_fill_manual(
    values = c("REF" = "#2166ac", "Dataset_B" = "#1b7837", "Dataset_E" = "#762a83"),
    na.value = "grey90"
  ) +
  labs(
    title = "Geographical coverage (REG) by sample",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# ============================
# FIGURE 2 — Auxiliary distributions
# ============================
aux_long <- bind_rows(
  ref %>% mutate(sample = "REF", weight = .data[[REF_W_COL]]),
  df_B %>% mutate(sample = "Dataset_B", weight = 1),
  df_E %>% mutate(sample = "Dataset_E", weight = 1)
) %>%
  pivot_longer(all_of(AUX_VARS), names_to = "var", values_to = "level")

p_aux <- ggplot(aux_long, aes(x = factor(level), weight = weight, fill = sample)) +
  geom_bar(position = "fill") +
  facet_wrap(~ var, scales = "free_x", nrow = 1) +
  scale_y_continuous(labels = percent) +
  labs(
    title = "Auxiliary-variable distributions",
    subtitle = "Reference weighted, NP samples unweighted",
    x = "Level",
    y = "Share",
    fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ============================
# FIGURE 3 — Income (EDI) distributions
# ============================
income_long <- bind_rows(
  ref %>% transmute(sample = "REF", EDI = .data[[Y_COL]], weight = .data[[REF_W_COL]]),
  df_B %>% transmute(sample = "Dataset_B", EDI = .data[[Y_COL]], weight = 1),
  df_E %>% transmute(sample = "Dataset_E", EDI = .data[[Y_COL]], weight = 1)
)

p_income <- ggplot(income_long, aes(x = EDI, weight = weight, color = sample)) +
  geom_density(linewidth = 1.1) +
  labs(
    title = "Income (EDI) distribution by sample",
    x = "EDI",
    y = "Density",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ============================
# Show figures
# ============================
p_maps
p_aux
p_income

# ============================
# Save (optional)
# ============================
ggsave("fig_maps_reg.png",   p_maps,   width = 12, height = 4, dpi = 300)
ggsave("fig_aux_dist.png",   p_aux,    width = 12, height = 4, dpi = 300)
ggsave("fig_income_dist.png",p_income, width = 8,  height = 4, dpi = 300)



reg_compare <- bind_rows(
  df_B %>% count(REG) %>% mutate(sample = "B"),
  df_E %>% count(REG) %>% mutate(sample = "E")
) %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = sample, values_from = prop) %>%
  mutate(diff_pp = 100 * (B - E))

ggplot(reg_compare, aes(x = factor(REG), y = diff_pp)) +
  geom_col(fill = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(
    title = "REG distribution difference (Dataset B − Dataset E)",
    x = "REG",
    y = "Difference (pp)"
  ) +
  theme_minimal(base_size = 12)


aux_diff <- function(var) {
  bind_rows(
    df_B %>% count(level = .data[[var]]) %>% mutate(sample = "B"),
    df_E %>% count(level = .data[[var]]) %>% mutate(sample = "E")
  ) %>%
    group_by(sample) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = sample, values_from = prop) %>%
    mutate(
      var = var,
      diff_pp = 100 * (B - E)
    )
}

aux_diff_df <- bind_rows(
  aux_diff("REG"),
  aux_diff("UEP"),
  aux_diff("DOU")
)

ggplot(aux_diff_df, aes(x = factor(level), y = diff_pp)) +
  geom_col(fill = "grey40") +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  facet_wrap(~ var, scales = "free_x", nrow = 1) +
  labs(
    title = "Auxiliary-variable differences (Dataset B − Dataset E)",
    x = "Level",
    y = "Difference (pp)"
  ) +
  theme_minimal(base_size = 12)


tv_distance <- function(var) {
  tab <- bind_rows(
    df_B %>% count(level = .data[[var]]) %>% mutate(sample = "B"),
    df_E %>% count(level = .data[[var]]) %>% mutate(sample = "E")
  ) %>%
    group_by(sample) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from = sample, values_from = prop)

  0.5 * sum(abs(tab$B - tab$E))
}

data.frame(
  variable = c("REG", "UEP", "DOU"),
  TV = c(
    tv_distance("REG"),
    tv_distance("UEP"),
    tv_distance("DOU")
  )
)


# 33333333---------------------



# ----------------------------
# REG shares per sample
# ----------------------------
reg_share <- function(df, label) {
  df %>%
    count(REG) %>%
    mutate(
      share = n / sum(n),
      sample = label
    ) %>%
    select(REG, share, sample)
}

reg_shares <- bind_rows(
  reg_share(df_B, "Dataset_B"),
  reg_share(df_E, "Dataset_E")
)

AMAP_reg_share <- AMAP_sf %>%
  left_join(reg_shares, by = "REG")

# ----------------------------
# Informative REG maps
# ----------------------------
p_maps <- ggplot(AMAP_reg_share) +
  geom_sf(aes(fill = share), color = "white", linewidth = 0.1) +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_viridis_c(
    option = "C",
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    title = "Regional concentration of non-probability samples",
    fill = "Share of sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p_maps



# ---------------


# ============================
# Auxiliary-difference plot
#   Dataset B − Dataset E
# ============================

plot_aux_diff_BE <- function(df_B, df_E, aux_vars,
                             ylim_pp = c(-3, 3),
                             neutral_band_pp = 1) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)

  aux_diff <- function(var) {
    bind_rows(
      df_B %>%
        filter(!is.na(.data[[var]])) %>%
        count(level = .data[[var]]) %>%
        mutate(sample = "B"),
      df_E %>%
        filter(!is.na(.data[[var]])) %>%
        count(level = .data[[var]]) %>%
        mutate(sample = "E")
    ) %>%
      group_by(sample) %>%
      mutate(prop = n / sum(n)) %>%
      ungroup() %>%
      select(-n) %>%
      pivot_wider(names_from = sample, values_from = prop) %>%
      mutate(
        var = var,
        diff_pp = 100 * (B - E)
      )
  }

  aux_diff_df <- bind_rows(lapply(aux_vars, aux_diff))

  ggplot(aux_diff_df, aes(x = factor(level), y = diff_pp, fill = diff_pp)) +
    geom_rect(
      xmin = -Inf, xmax = Inf,
      ymin = -neutral_band_pp, ymax = neutral_band_pp,
      fill = "grey90", alpha = 0.6
    ) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    facet_wrap(~ var, scales = "free_y", nrow = 1) +
    scale_fill_gradient2(
      low = "#b2182b",
      mid = "grey70",
      high = "#2166ac",
      midpoint = 0
    ) +
    labs(
      title = "Auxiliary-variable differences (Dataset B − Dataset E)",
      subtitle = "Grey band: differences within ±1 percentage point",
      x = "Level",
      y = "Difference (pp)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )

}

# ============================
# USAGE
# ============================

AUX_VARS <- c("REG", "UEP", "DOU")

p_aux_diff <- plot_aux_diff_BE(
  df_B = df_B,
  df_E = df_E,
  aux_vars = AUX_VARS
)

print(p_aux_diff)
