# -------------------------
# Basic helpers
# -------------------------

assert_required_cols <- function(df, cols, fn_name = "function") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf(
      "%s: missing required columns: %s",
      fn_name, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

coerce_to_numeric <- function(df, cols) {
  for (c in cols) {
    if (c %in% names(df)) {
      df[[c]] <- as.numeric(df[[c]])
    }
  }
  df
}

coerce_to_factor <- function(df, cols) {
  for (c in cols) {
    if (c %in% names(df)) {
      df[[c]] <- as.factor(df[[c]])
    }
  }
  df
}

# -------------------------
# Loading convenience (optional)
# -------------------------

load_extdata_csv <- function(filename, package, ...) {
  path <- system.file("extdata", package = package)
  if (path == "") stop("Package extdata path not found. Is the package installed?", call. = FALSE)
  read.csv(file.path(path, filename), ...)
}

# -------------------------
# Missingness
# -------------------------

eda_missing_summary <- function(df) {
  data.frame(
    variable = names(df),
    missing = colSums(is.na(df)),
    missing_rate = colMeans(is.na(df)),
    row.names = NULL
  )
}

plot_missingness <- function(df, name, save = FALSE) {
  ms <- eda_missing_summary(df)
  p <- ggplot2::ggplot(ms, ggplot2::aes(x = stats::reorder(variable, missing_rate), y = missing_rate)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste("Missingness by variable —", name),
      x = NULL,
      y = "Missing rate"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_missingness.png"), save)
  p
}

# -------------------------
# Income transforms + summary stats
# -------------------------

add_signed_log_income <- function(df, income_var = "EDI", out_var = "signed_log_EDI") {
  assert_required_cols(df, income_var, "add_signed_log_income")
  df[[out_var]] <- sign(df[[income_var]]) * log(abs(df[[income_var]]) + 1)
  df
}

income_sign_stats <- function(df, income_var = "EDI") {
  assert_required_cols(df, income_var, "income_sign_stats")
  x <- df[[income_var]]
  data.frame(
    n = length(x),
    n_leq_0 = sum(x <= 0, na.rm = TRUE),
    n_lt_0 = sum(x < 0, na.rm = TRUE),
    share_leq_0 = mean(x <= 0, na.rm = TRUE),
    share_lt_0 = mean(x < 0, na.rm = TRUE),
    row.names = NULL
  )
}

# -------------------------
# Distributions (overall + facet)
# -------------------------

plot_income_distribution <- function(df, name, save = FALSE, signed_log_var = "signed_log_EDI", bins = 200) {
  assert_required_cols(df, signed_log_var, "plot_income_distribution")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[signed_log_var]])) +
    ggplot2::geom_histogram(bins = bins, fill = "grey30") +
    ggplot2::labs(
      title = paste("Signed-log income distribution —", name),
      x = signed_log_var,
      y = "Count"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_income_distribution.png"), save)
  p
}

plot_income_distribution_by <- function(df, name, facet_var, save = FALSE, signed_log_var = "signed_log_EDI", bins = 200, ncol = 2) {
  assert_required_cols(df, c(signed_log_var, facet_var), "plot_income_distribution_by")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[signed_log_var]])) +
    ggplot2::geom_histogram(bins = bins, fill = "grey30") +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_var)), ncol = ncol) +
    ggplot2::labs(
      title = paste("Signed-log income distribution by", facet_var, "—", name),
      x = signed_log_var,
      y = "Count"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_income_by_", facet_var, ".png"), save)
  p
}

# -------------------------
# Boxplots for income by categorical variable
# -------------------------

plot_income_box_by <- function(df, name, group_var, save = FALSE, income_var = "EDI") {
  assert_required_cols(df, c(group_var, income_var), "plot_income_box_by")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = as.factor(.data[[group_var]]), y = .data[[income_var]])) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::labs(
      title = paste(income_var, "by", group_var, "—", name),
      x = group_var,
      y = income_var
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_box_", group_var, ".png"), save)
  p
}

plot_signed_log_box_by_two_groups <- function(df, name, x_var, fill_var, save = FALSE, signed_log_var = "signed_log_EDI") {
  assert_required_cols(df, c(x_var, fill_var, signed_log_var), "plot_signed_log_box_by_two_groups")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = as.factor(.data[[x_var]]), y = .data[[signed_log_var]], fill = as.factor(.data[[fill_var]]))) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::labs(
      title = paste(signed_log_var, "by", x_var, "and", fill_var, "—", name),
      x = x_var,
      y = signed_log_var
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_signedlog_box_", x_var, "_fill_", fill_var, ".png"), save)
  p
}

# -------------------------
# Categorical proportions
# -------------------------

categorical_proportions <- function(df, var) {
  assert_required_cols(df, var, "categorical_proportions")
  tab <- table(df[[var]], useNA = "ifany")
  prop <- prop.table(tab)
  data.frame(
    level = names(tab),
    n = as.integer(tab),
    prop = as.numeric(prop),
    row.names = NULL
  )
}

plot_categorical_proportions <- function(df, name, var, save = FALSE) {
  cp <- categorical_proportions(df, var)

  p <- ggplot2::ggplot(cp, ggplot2::aes(x = as.factor(level), y = prop)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = paste("Proportions for", var, "—", name),
      x = var,
      y = "Proportion"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_prop_", var, ".png"), save)
  p
}

# -------------------------
# Common support / overlap checks (REG × UEP × DOU etc.)
# -------------------------

common_support_counts <- function(df, vars = c("REG", "UEP", "DOU")) {
  assert_required_cols(df, vars, "common_support_counts")
  dplyr::count(df, dplyr::across(dplyr::all_of(vars)), name = "n") |>
    dplyr::arrange(n)
}

flag_low_support <- function(counts_df, threshold = 30) {
  assert_required_cols(counts_df, "n", "flag_low_support")
  dplyr::mutate(counts_df, low_support = .data$n < threshold)
}

support_report <- function(df, name, vars = c("REG", "UEP", "DOU"), thresholds = c(100, 30, 10)) {
  counts <- common_support_counts(df, vars)
  out <- list(counts = counts)
  for (t in thresholds) {
    out[[paste0("n_lt_", t)]] <- sum(counts$n < t)
  }
  out$name <- name
  out$vars <- vars
  out
}

plot_support_heatmap <- function(df, name, save = FALSE, reg_var = "REG", x_var = "UEP", y_var = "DOU") {
  assert_required_cols(df, c(reg_var, x_var, y_var), "plot_support_heatmap")

  counts <- dplyr::count(df, .data[[reg_var]], .data[[x_var]], .data[[y_var]], name = "n")

  p <- ggplot2::ggplot(counts, ggplot2::aes(x = as.factor(.data[[x_var]]), y = as.factor(.data[[y_var]]), fill = n)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(stats::as.formula(paste("~", reg_var))) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(
      title = paste("Joint support:", reg_var, "×", x_var, "×", y_var, "—", name),
      x = x_var,
      y = y_var,
      fill = "Count"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_support_heatmap.png"), save)
  p
}

# -------------------------
# Correlation + multicollinearity (VIF)
# -------------------------

covariate_correlation <- function(df, vars = c("UEP", "DOU", "REG")) {
  assert_required_cols(df, vars, "covariate_correlation")
  tmp <- df
  tmp <- coerce_to_numeric(tmp, vars)
  stats::cor(tmp[, vars, drop = FALSE], use = "pairwise.complete.obs")
}

plot_correlation_heatmap <- function(df, name, vars = c("UEP", "DOU", "REG"), save = FALSE) {
  cm <- covariate_correlation(df, vars)
  m <- as.data.frame(as.table(cm))
  names(m) <- c("var1", "var2", "corr")

  p <- ggplot2::ggplot(m, ggplot2::aes(x = var1, y = var2, fill = corr)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2() +
    ggplot2::labs(
      title = paste("Correlation heatmap —", name),
      x = NULL, y = NULL, fill = "corr"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_cor_heatmap.png"), save)
  p
}

vif_check <- function(df, name, response = "signed_log_EDI", predictors = c("UEP", "DOU", "REG")) {
  assert_required_cols(df, c(response, predictors), "vif_check")

  if (!requireNamespace("car", quietly = TRUE)) {
    stop("Package 'car' is required for VIF. Install it first: install.packages('car')", call. = FALSE)
  }

  # Use factors as in your original approach
  formula_str <- paste0(response, " ~ ", paste(sprintf("factor(%s)", predictors), collapse = " + "))
  mod <- stats::lm(stats::as.formula(formula_str), data = df)
  v <- car::vif(mod)

  list(
    name = name,
    model = mod,
    vif = v
  )
}

# -------------------------
# Reference sample probability checks
# -------------------------

check_probabilities <- function(ref_df, prob_var) {
  assert_required_cols(ref_df, prob_var, "check_probabilities")
  p <- ref_df[[prob_var]]

  list(
    min = min(p, na.rm = TRUE),
    max = max(p, na.rm = TRUE),
    sum = sum(p, na.rm = TRUE),
    n_negative = sum(p < 0, na.rm = TRUE),
    n_zero = sum(p == 0, na.rm = TRUE),
    n_na = sum(is.na(p))
  )
}

plot_probability_distribution <- function(ref_df, name, prob_var, save = FALSE, bins = 100) {
  assert_required_cols(ref_df, prob_var, "plot_probability_distribution")

  p <- ggplot2::ggplot(ref_df, ggplot2::aes(x = .data[[prob_var]])) +
    ggplot2::geom_histogram(bins = bins, fill = "grey30") +
    ggplot2::labs(
      title = paste("Reference probabilities distribution —", name),
      x = prob_var,
      y = "Count"
    ) +
    ggplot2::theme_minimal()

  save_plot_if_needed(p, paste0(name, "_prob_dist_", prob_var, ".png"), save)
  p
}
