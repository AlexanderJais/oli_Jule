# Publication figures: one look for every panel (sans 7 pt, hairline axes, colour-blind-safe colours
# checked for colour-vision deficiency), panels composed with letters, and every figure saved as vector
# PDF plus 600-dpi PNG and TIFF (LZW) at journal widths (single column 88 mm, double column 180 mm).

fig_col <- c(ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781", grid = "#e1e0d9", light = "#c3c2b7",
             platelet = "#2a78d6", neuro = "#eb6834", third = "#1baf7a",           # categorical slots 1-3
             low = "#2a78d6", mid = "#f0efec", high = "#e34948")                   # diverging: blue - grey - red

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(text = element_text(colour = fig_col[["ink"]]),
          axis.text = element_text(colour = fig_col[["ink2"]]),
          axis.line = element_line(colour = fig_col[["ink2"]], linewidth = 0.3),
          axis.ticks = element_line(colour = fig_col[["ink2"]], linewidth = 0.3),
          panel.grid = element_blank(),
          legend.key.size = unit(3, "mm"), legend.margin = margin(0, 0, 0, 0),
          legend.text = element_text(colour = fig_col[["ink2"]]), legend.title = element_text(colour = fig_col[["ink2"]]),
          strip.background = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
          plot.title = element_text(face = "bold", size = rel(1), hjust = 0, margin = margin(b = 2)),
          plot.subtitle = element_text(colour = fig_col[["ink2"]], size = rel(0.9), hjust = 0, margin = margin(b = 3)),
          plot.margin = margin(4, 6, 4, 4))
}

#' Hairline row guides for forest plots (rows on the y axis).
row_guides <- function() theme(panel.grid.major.y = element_line(colour = fig_col[["grid"]], linewidth = 0.25))

#' Several ggplots as one figure with panel letters (A, B, ...).
compose <- function(plots, ncol = 2, rel_widths = 1, rel_heights = 1, labels = "AUTO", align = "none")
  cowplot::plot_grid(plotlist = plots, ncol = ncol, labels = labels, label_size = 9, label_fontface = "bold",
                     rel_widths = rel_widths, rel_heights = rel_heights, align = align, axis = "tblr")

#' Save a figure as PDF (vector), PNG and TIFF (600 dpi) under output/<dir>/<name>.*
save_figure <- function(p, cfg, dir, name, width_mm = 180, height_mm = 150) {
  base <- out_path(cfg, dir, name)
  pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggsave(paste0(base, ".pdf"), p, width = width_mm, height = height_mm, units = "mm", device = pdf_dev, bg = "white")
  ggsave(paste0(base, ".png"), p, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white")
  ggsave(paste0(base, ".tiff"), p, width = width_mm, height = height_mm, units = "mm", dpi = 600, bg = "white", compression = "lzw")
  invisible(base)
}

#' "rho = 0.45 (95% CI 0.13 to 0.68), p = 0.0071" from a spearman1() row.
rho_txt <- function(r, ci = TRUE)
  if (!nrow(r) || is.na(r$rho)) "n/a" else
    paste0(sprintf("rho = %.2f", r$rho), if (ci) sprintf(" (95%% CI %.2f to %.2f)", r$ci_low, r$ci_high) else "", sprintf(", p = %s", fmt_p(r$p)))
