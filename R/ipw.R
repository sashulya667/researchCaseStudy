# ============================================================
# R/ipw.R
# Simple + correct IPW helpers for nonprob sample using a
# probability reference sample with design weights.
#
# Assumes:
# - factor levels are harmonized in analysis
# - missing X rows are removed in analysis
# ============================================================

# -------------------------
# Step 1: Build combined training data
# -------------------------
make_ipw_training_data <- function(nonprob_df,
                                   reference_df,
                                   ref_weight_var = "design_weight") {

  np <- dplyr::mutate(nonprob_df, S = 1L)
  rf <- dplyr::mutate(reference_df, S = 0L)

  combined <- dplyr::bind_rows(np, rf) |>
    dplyr::mutate(
      w_fit = dplyr::if_else(S == 0L, as.numeric(.data[[ref_weight_var]]), 1.0)
    )

  list(
    combined = combined,
    n_nonprob = nrow(np),
    n_ref = nrow(rf)
  )
}

# -------------------------
# Separate step: aggregate to cells
# -------------------------
make_cell_counts <- function(combined,
                             x_vars = c("REG", "UEP", "DOU"),
                             s_var = "S",
                             weight_var = "w_fit") {

  combined |>
    dplyr::group_by(dplyr::across(dplyr::all_of(x_vars))) |>
    dplyr::summarise(
      y1 = sum(.data[[s_var]] == 1L, na.rm = TRUE),
      y0 = sum(.data[[weight_var]] * (.data[[s_var]] == 0L), na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------------
# Step 2a: Fit selection model (row-level)
# -------------------------
fit_selection_logit_row <- function(combined, formula, weight_var = "w_fit") {
  combined2 <- combined
  combined2[[".__w__"]] <- as.numeric(combined2[[weight_var]])

  stats::glm(
    formula = formula,
    data = combined2,
    family = stats::binomial(),
    weights = .__w__
  )
}

# -------------------------
# Step 2b: Fit selection model (grouped cells)
# -------------------------
fit_selection_logit_grouped <- function(cell, rhs_formula) {
  rhs <- paste(deparse(rhs_formula[[2]]), collapse = "")
  f2 <- stats::as.formula(paste0("cbind(y1, y0) ~ ", rhs))

  stats::glm(
    formula = f2,
    data = cell,
    family = stats::binomial()
  )
}

# -------------------------
# Step 3: Predict propensities
# -------------------------
predict_p_hat <- function(fit, df_new) {
  if (!is.null(fit$xlevels)) {
    for (v in names(fit$xlevels)) {
      if (v %in% names(df_new)) {
        df_new[[v]] <- factor(df_new[[v]], levels = fit$xlevels[[v]])
      }
    }
  }
  as.numeric(stats::predict(fit, newdata = df_new, type = "response"))
}

# -------------------------
# Step 3b: Convert propensities to IPW
# -------------------------
ipw_from_p <- function(p_hat, eps = 1e-6) {
  # clipping prevents Inf/NaN weights when p_hat is 0/1 or numerically extreme
  p <- pmin(pmax(as.numeric(p_hat), eps), 1 - eps)
  (1-p) / p
}

attach_ipw <- function(nonprob_df, p_hat, w_ipw,
                       p_col = "p_hat", w_col = "w_ipw") {
  dplyr::mutate(
    nonprob_df,
    !!p_col := as.numeric(p_hat),
    !!w_col := as.numeric(w_ipw)
  )
}

# -------------------------
# Diagnostics
# -------------------------
weight_summary <- function(w) {
  w <- as.numeric(w)
  data.frame(
    n = length(w),
    min = min(w, na.rm = TRUE),
    p50 = stats::quantile(w, 0.50, na.rm = TRUE),
    mean = mean(w, na.rm = TRUE),
    p95 = stats::quantile(w, 0.95, na.rm = TRUE),
    p99 = stats::quantile(w, 0.99, na.rm = TRUE),
    max = max(w, na.rm = TRUE),
    sd = stats::sd(w, na.rm = TRUE),
    row.names = NULL
  )
}

effective_sample_size <- function(w) {
  w <- as.numeric(w)
  (sum(w, na.rm = TRUE)^2) / sum(w^2, na.rm = TRUE)
}

plot_overlap_p <- function(p_hat_all, S_all, title = "Propensity overlap") {
  df <- tibble::tibble(
    p = as.numeric(p_hat_all),
    source = factor(S_all, levels = c(0, 1), labels = c("Reference", "Nonprob"))
  )
  ggplot2::ggplot(df, ggplot2::aes(x = p, fill = source)) +
    ggplot2::geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
    ggplot2::labs(title = title, x = "p_hat = P(Nonprob | X)", y = "Count") +
    ggplot2::theme_minimal()
}

plot_weights <- function(w, title = "IPW weights") {
  ggplot2::ggplot(tibble::tibble(w = as.numeric(w)), ggplot2::aes(x = w)) +
    ggplot2::geom_histogram(bins = 120, fill = "grey30") +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = title, x = "Weight (log10 scale)", y = "Count") +
    ggplot2::theme_minimal()
}

trim_weights <- function(w, cap_q = 0.99) {
  cap <- stats::quantile(as.numeric(w), cap_q, na.rm = TRUE)
  pmin(as.numeric(w), cap)
}

# -------------------------
# AROP
# -------------------------
weighted_quantile <- function(x, w, p = 0.5) {
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & w > 0
  x <- x[ok]; w <- w[ok]
  if (length(x) == 0) stop("weighted_quantile: no valid data", call. = FALSE)

  o <- order(x)
  x <- x[o]; w <- w[o]
  cw <- cumsum(w)
  x[which(cw >= p * sum(w))[1]]
}

arop_unweighted <- function(df, income_var = "EDI", poverty_frac = 0.6) {
  x <- df[[income_var]]
  x <- x[is.finite(x) & !is.na(x)]
  med <- stats::median(x)
  line <- poverty_frac * med
  arop <- mean(x < line)
  tibble::tibble(estimator = "unweighted", median = med, poverty_line = line, arop = arop)
}

arop_weighted <- function(df, income_var = "EDI", weight_var, poverty_frac = 0.6) {
  x <- df[[income_var]]
  w <- df[[weight_var]]
  med <- weighted_quantile(x, w, p = 0.5)
  line <- poverty_frac * med

  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & w > 0
  poor <- as.numeric(x < line)
  arop <- sum(w[ok] * poor[ok]) / sum(w[ok])
  tibble::tibble(estimator = paste0("weighted_", weight_var), median = med, poverty_line = line, arop = arop)
}

# -------------------------
# Overlap diagnostics
# -------------------------
compute_overlap <- function(combined, x_vars = c("REG", "UEP", "DOU"),
                            s_var = "S", w_var = "w_fit") {

  cnt_np <- combined |>
    dplyr::filter(.data[[s_var]] == 1) |>
    dplyr::count(dplyr::across(dplyr::all_of(x_vars)), name = "n_np")

  cnt_ref <- combined |>
    dplyr::filter(.data[[s_var]] == 0) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(x_vars))) |>
    dplyr::summarise(N_ref = sum(.data[[w_var]], na.rm = TRUE), .groups = "drop")

  overlap <- dplyr::full_join(cnt_np, cnt_ref, by = x_vars) |>
    dplyr::mutate(
      n_np = tidyr::replace_na(n_np, 0L),
      N_ref = tidyr::replace_na(N_ref, 0),
      bad_no_ref_support = (n_np > 0 & N_ref == 0),
      bad_no_np_support  = (n_np == 0 & N_ref > 0)
    )

  list(
    overlap = overlap,
    bad_no_ref = overlap |> dplyr::filter(bad_no_ref_support),
    bad_no_np  = overlap |> dplyr::filter(bad_no_np_support),
    summary = overlap |> dplyr::summarise(
      cells_total = dplyr::n(),
      cells_no_ref_support = sum(bad_no_ref_support),
      cells_no_np_support  = sum(bad_no_np_support)
    )
  )
}

drop_bad_cells <- function(df, bad_cells, x_vars = c("REG", "UEP", "DOU")) {
  if (nrow(bad_cells) == 0) return(list(df = df, n_dropped = 0L))
  before <- nrow(df)
  df2 <- df |> dplyr::anti_join(bad_cells |> dplyr::select(dplyr::all_of(x_vars)), by = x_vars)
  list(df = df2, n_dropped = before - nrow(df2))
}

msg <- function(...) cat(paste0(..., "\n"))


# ============================================================
# ROW-GLM NUMERICAL FAILURE DIAGNOSTICS
# ============================================================

get_glm <- function(obj) {
  if (inherits(obj, "glm")) return(obj)
  if (!is.null(obj$fit) && inherits(obj$fit, "glm")) return(obj$fit)
  if (!is.null(obj$model) && inherits(obj$model, "glm")) return(obj$model)
  stop("Can't find glm object inside fit_*_row. Expected obj, obj$fit, or obj$model to be a glm.")
}

diagnose_row_glm <- function(fit_obj, combined_df, x_vars, s_var = "S", w_var = "w_fit",
                             eps = 1e-12, top_k = 10, label = "ROW") {
  g <- get_glm(fit_obj)

  cat("\n====================\n")
  cat("DIAG:", label, "\n")
  cat("====================\n")

  # ---- A) Boundary fitted probabilities (the actual warning quantified)
  p <- as.numeric(g$fitted.values)
  eta <- as.numeric(g$linear.predictors)

  cat("\n[A] Fitted p-hat boundary check\n")
  print(tibble::tibble(
    p_min = min(p, na.rm = TRUE),
    p_q001 = unname(quantile(p, 0.001, na.rm = TRUE)),
    p_q01  = unname(quantile(p, 0.01,  na.rm = TRUE)),
    p_med  = unname(quantile(p, 0.5,   na.rm = TRUE)),
    p_q99  = unname(quantile(p, 0.99,  na.rm = TRUE)),
    p_q999 = unname(quantile(p, 0.999, na.rm = TRUE)),
    p_max = max(p, na.rm = TRUE),
    n_p_lt_eps = sum(p < eps, na.rm = TRUE),
    n_p_gt_1m  = sum(p > 1 - eps, na.rm = TRUE)
  ))

  cat("\n[A2] Linear predictor range (η): big |η| => p numerically 0/1\n")
  print(tibble::tibble(
    eta_min = min(eta, na.rm = TRUE),
    eta_q01 = unname(quantile(eta, 0.01,  na.rm = TRUE)),
    eta_med = unname(quantile(eta, 0.5,   na.rm = TRUE)),
    eta_q99 = unname(quantile(eta, 0.99,  na.rm = TRUE)),
    eta_max = max(eta, na.rm = TRUE),
    n_abs_eta_gt_20 = sum(abs(eta) > 20, na.rm = TRUE),
    n_abs_eta_gt_30 = sum(abs(eta) > 30, na.rm = TRUE)
  ))

  # ---- B) Coefficient blow-up / instability symptoms
  cat("\n[B] Coefficient magnitude + SE (large/NA is typical under separation)\n")
  co <- summary(g)$coefficients
  co_tbl <- tibble::tibble(
    term = rownames(co),
    estimate = co[,1],
    se = co[,2],
    z = co[,3],
    p_value = co[,4]
  ) |>
    dplyr::mutate(abs_est = abs(estimate), is_na = is.na(estimate) | is.na(se)) |>
    dplyr::arrange(dplyr::desc(abs_est))

  print(dplyr::summarise(
    co_tbl,
    n_terms = dplyr::n(),
    n_na = sum(is_na),
    max_abs_est = max(abs_est, na.rm = TRUE),
    n_abs_est_gt_10 = sum(abs_est > 10, na.rm = TRUE),
    n_abs_est_gt_20 = sum(abs_est > 20, na.rm = TRUE)
  ))
  cat("\nTop terms by |estimate|:\n")
  print(dplyr::slice_head(co_tbl, n = 15))

  # ---- C) Rank deficiency / aliasing (another common numeric issue)
  cat("\n[C] Model matrix rank / aliasing\n")
  X <- model.matrix(g)
  qrX <- qr(X)
  cat("ncol(X) =", ncol(X), "  rank(X) =", qrX$rank, "  (rank < ncol => non-identifiability)\n")
  al <- tryCatch(alias(g)$Complete, error = function(e) NULL)
  if (!is.null(al)) {
    cat("Aliased (non-estimable) coefficients detected.\n")
  } else {
    cat("No complete aliasing detected (doesn't rule out quasi-separation).\n")
  }

  # ---- D) Cell-level proof: some cells imply extreme selection odds => p ~ 0/1
  cat("\n[D] Cell-level selection odds diagnostics (REG×UEP×DOU)\n")

  # recompute p-hat on the training rows, then aggregate by cell
  combined_df2 <- combined_df |>
    dplyr::mutate(p_hat_row = as.numeric(predict(g, type = "response")))

  cell <- combined_df2 |>
    dplyr::group_by(dplyr::across(dplyr::all_of(x_vars))) |>
    dplyr::summarise(
      n_total = dplyr::n(),
      n_np = sum(.data[[s_var]] == 1),
      ref_w_sum = sum(dplyr::if_else(.data[[s_var]] == 0, .data[[w_var]], 0), na.rm = TRUE),
      # pseudo "empirical" selection rate (weighted ref + unweighted nonprob)
      p_emp = n_np / (n_np + ref_w_sum),
      p_hat = dplyr::first(p_hat_row),
      p_hat_min = min(p_hat_row),
      p_hat_max = max(p_hat_row),
      p_hat_sd  = sd(p_hat_row),
      .groups = "drop"
    ) |>
    # stabilized logit of empirical rate to avoid Inf
    dplyr::mutate(
      p_emp_stab = pmin(pmax(p_emp, eps), 1 - eps),
      logit_emp = log(p_emp_stab / (1 - p_emp_stab)),
      abs_logit_emp = abs(logit_emp),
      near_boundary_hat = (p_hat < eps) | (p_hat > 1 - eps)
    ) |>
    dplyr::arrange(dplyr::desc(abs_logit_emp))

  cat("\nCells with most extreme empirical log-odds (these drive η to ±∞):\n")
  print(dplyr::slice_head(cell, n = top_k))

  cat("\nHow many cells have p-hat numerically at boundary?\n")
  print(dplyr::summarise(cell,
                         n_cells = dplyr::n(),
                         n_cells_near_boundary_hat = sum(near_boundary_hat, na.rm = TRUE),
                         max_abs_logit_emp = max(abs_logit_emp, na.rm = TRUE),
                         max_p_hat_sd_within_cell = max(p_hat_sd, na.rm = TRUE)
  ))

  # ---- E) Direct consequence: implied IPW weight explosion (if ipw = (1-p)/p)
  cat("\n[E] Consequence for IPW weights from ROW p-hat (w = (1-p)/p)\n")
  w <- (1 - p) / p
  w <- w[is.finite(w) & !is.na(w) & w > 0]
  print(tibble::tibble(
    w_min = min(w),
    w_q50 = unname(quantile(w, 0.5)),
    w_q90 = unname(quantile(w, 0.9)),
    w_q99 = unname(quantile(w, 0.99)),
    w_q999 = unname(quantile(w, 0.999)),
    w_max = max(w),
    ess = (sum(w)^2) / sum(w^2)
  ))

  invisible(list(glm = g, cell = cell))
}

