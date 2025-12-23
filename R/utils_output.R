ensure_output_dir <- function(subdir = "figures") {
  path <- file.path("outputs", subdir)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  invisible(path)
}

save_plot_if_needed <- function(plot, filename, save) {
  if (save) {
    ensure_output_dir("figures")
    ggsave(
      filename = file.path("outputs", "figures", filename),
      plot = plot,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
}
