devtools::load_all()

rm(list = ls())

# ---- Settings ----
PACKAGE_NAME <- "researchCaseStudy"
SAVE_PLOTS <- TRUE    # set FALSE if you don't want files saved
BINS <- 200

# ---- Libraries ----
library(tidyverse)
library(car)

# ---- Helper to load from package inst/extdata ----
load_extdata <- function(filename) {
  path <- here::here("data")
  read.csv(file.path(path, filename))
}


# ---- Load datasets ----
df_B <- load_extdata("Dataset_B.csv")
df_E <- load_extdata("Dataset_E.csv")

reference_sample_probabilities = load_extdata("reference_sample_probabilistic.csv")

str(reference_sample_probabilities)

cat("Loaded Dataset_B:", nrow(df_B), "rows,", ncol(df_B), "cols\n")
cat("Loaded Dataset_E:", nrow(df_E), "rows,", ncol(df_E), "cols\n")
str(df_E)
str(df_B)
# ---- Check required columns ----
required_cols <- c("EDI", "REG", "UEP", "DOU")
missing_B <- setdiff(required_cols, names(df_B))
missing_E <- setdiff(required_cols, names(df_E))

if (length(missing_B) > 0) stop("Dataset_B missing: ", paste(missing_B, collapse = ", "))
if (length(missing_E) > 0) stop("Dataset_E missing: ", paste(missing_E, collapse = ", "))

# ---- Add signed-log EDI ----
df_B <- add_signed_log_income(df_B, income_var = "EDI", out_var = "signed_log_EDI")
df_E <- add_signed_log_income(df_E, income_var = "EDI", out_var = "signed_log_EDI")

# ============================================================
# 1) MISSINGNESS
# ============================================================

cat("\n===== Missingness summary: Dataset_B =====\n")
print(eda_missing_summary(df_B) %>% arrange(desc(missing)) %>% head(20))

cat("\n===== Missingness summary: Dataset_E =====\n")
print(eda_missing_summary(df_E) %>% arrange(desc(missing)) %>% head(20))

p <- plot_missingness(df_B, name = "Dataset_B", save = SAVE_PLOTS)
print(p)
p <- plot_missingness(df_E, name = "Dataset_E", save = SAVE_PLOTS)
print(p)

# ============================================================
# 2) INCOME DISTRIBUTIONS
# ============================================================

cat("\n===== Income sign stats: Dataset_B =====\n")
print(income_sign_stats(df_B, income_var = "EDI"))

cat("\n===== Income sign stats: Dataset_E =====\n")
print(income_sign_stats(df_E, income_var = "EDI"))

p <- plot_income_distribution(df_B, name = "Dataset_B", save = SAVE_PLOTS, bins = BINS)
print(p)
p <- plot_income_distribution(df_E, name = "Dataset_E", save = SAVE_PLOTS, bins = BINS)
print(p)

p <- plot_income_distribution_by(df_B, name = "Dataset_B", facet_var = "REG", save = SAVE_PLOTS, bins = BINS)
print(p)
p <- plot_income_distribution_by(df_E, name = "Dataset_E", facet_var = "REG", save = SAVE_PLOTS, bins = BINS)
print(p)

# ============================================================
# 3) REGION / DOU BOX PLOTS
# ============================================================

p <- plot_income_box_by(df_B, name = "Dataset_B", group_var = "REG", save = SAVE_PLOTS, income_var = "EDI")
print(p)
p <- plot_income_box_by(df_E, name = "Dataset_E", group_var = "REG", save = SAVE_PLOTS, income_var = "EDI")
print(p)

p <- plot_income_box_by(df_B, name = "Dataset_B", group_var = "DOU", save = SAVE_PLOTS, income_var = "EDI")
print(p)
p <- plot_income_box_by(df_E, name = "Dataset_E", group_var = "DOU", save = SAVE_PLOTS, income_var = "EDI")
print(p)

p <- plot_signed_log_box_by_two_groups(df_B, name = "Dataset_B", x_var = "DOU", fill_var = "REG", save = SAVE_PLOTS)
print(p)
p <- plot_signed_log_box_by_two_groups(df_E, name = "Dataset_E", x_var = "DOU", fill_var = "REG", save = SAVE_PLOTS)
print(p)

# ============================================================
# 4) CATEGORICAL PROPORTIONS (REG, DOU, UEP)
# ============================================================

cat("\n===== Proportions: Dataset_B =====\n")
print(categorical_proportions(df_B, "REG"))
print(categorical_proportions(df_B, "DOU"))
print(categorical_proportions(df_B, "UEP"))

cat("\n===== Proportions: Dataset_E =====\n")
print(categorical_proportions(df_E, "REG"))
print(categorical_proportions(df_E, "DOU"))
print(categorical_proportions(df_E, "UEP"))

p <- plot_categorical_proportions(df_B, "Dataset_B", "REG", save = SAVE_PLOTS); print(p)
p <- plot_categorical_proportions(df_B, "Dataset_B", "DOU", save = SAVE_PLOTS); print(p)
p <- plot_categorical_proportions(df_B, "Dataset_B", "UEP", save = SAVE_PLOTS); print(p)

p <- plot_categorical_proportions(df_E, "Dataset_E", "REG", save = SAVE_PLOTS); print(p)
p <- plot_categorical_proportions(df_E, "Dataset_E", "DOU", save = SAVE_PLOTS); print(p)
p <- plot_categorical_proportions(df_E, "Dataset_E", "UEP", save = SAVE_PLOTS); print(p)

# ============================================================
# 5) COMMON SUPPORT (REG × UEP × DOU)
# ============================================================

cat("\n===== Common support counts (sorted): Dataset_B =====\n")
support_B <- common_support_counts(df_B, vars = c("REG", "UEP", "DOU"))
print(head(support_B, 30))

cat("\nCounts below thresholds (Dataset_B):\n")
print(list(
  n_lt_100 = sum(support_B$n < 100),
  n_lt_30  = sum(support_B$n < 30),
  n_lt_10  = sum(support_B$n < 10)
))

cat("\n===== Common support counts (sorted): Dataset_E =====\n")
support_E <- common_support_counts(df_E, vars = c("REG", "UEP", "DOU"))
print(head(support_E, 30))

cat("\nCounts below thresholds (Dataset_E):\n")
print(list(
  n_lt_100 = sum(support_E$n < 100),
  n_lt_30  = sum(support_E$n < 30),
  n_lt_10  = sum(support_E$n < 10)
))

p <- plot_support_heatmap(df_B, name = "Dataset_B", save = SAVE_PLOTS)
print(p)
p <- plot_support_heatmap(df_E, name = "Dataset_E", save = SAVE_PLOTS)
print(p)

# ============================================================
# 6) CORRELATION + VIF CHECKS
# ============================================================

cat("\n===== Correlation matrix: Dataset_B (UEP, DOU, REG) =====\n")
print(covariate_correlation(df_B, vars = c("UEP", "DOU", "REG")))

cat("\n===== Correlation matrix: Dataset_E (UEP, DOU, REG) =====\n")
print(covariate_correlation(df_E, vars = c("UEP", "DOU", "REG")))

p <- plot_correlation_heatmap(df_B, name = "Dataset_B", vars = c("UEP", "DOU", "REG"), save = SAVE_PLOTS)
print(p)
p <- plot_correlation_heatmap(df_E, name = "Dataset_E", vars = c("UEP", "DOU", "REG"), save = SAVE_PLOTS)
print(p)

# VIF (optional; requires 'car')
vB <- vif_check(df_B, name = "Dataset_B", response = "signed_log_EDI", predictors = c("UEP", "DOU", "REG"))
print(vB$vif)
vE <- vif_check(df_E, name = "Dataset_E", response = "signed_log_EDI", predictors = c("UEP", "DOU", "REG"))
print(vE$vif)

##################### weight interpretations #########################

# ============================================================
# IPW DIAGNOSTICS USING REFERENCE SAMPLE
# ============================================================

# Drop negligible missingness
df_E_ipw <- df_E |> filter(!is.na(UEP))
ref_ipw <- reference_sample_probabilities |> filter(!is.na(UEP))

# ---- Estimate selection propensity ----
ps_obj <- estimate_selection_propensity(
  df = df_E_ipw,
  ref_df = ref_ipw,
  formula = source ~ factor(REG) + factor(UEP) + factor(DOU)
)

ps <- ps_obj$ps
source <- ps_obj$combined$source

# ---- Propensity overlap ----
p <- plot_propensity_histogram(
  ps = ps,
  source = source,
  name = "Dataset_E",
  save = SAVE_PLOTS
)
print(p)

# ---- IPW weights ----
w_raw <- compute_ipw_weights(ps, source)
w_trim <- trim_weights(w_raw, p = 0.99)

# ---- Weight distributions ----
p <- plot_weight_distribution(w_raw, "Dataset_E", "raw", save = SAVE_PLOTS)
print(p)

p <- plot_weight_distribution(w_trim, "Dataset_E", "trimmed", save = SAVE_PLOTS)
print(p)

# ---- ESS ----
ess <- ess_table(
  w_raw = w_raw,
  w_trimmed = w_trim,
  n = nrow(df_E_ipw)
)

cat("\n===== Effective Sample Size =====\n")
print(ess)

