# Helpers for the per-visit analysis (visits_v1_v4/): figures for slides - 16:9, large fonts, saved as PNG (300 dpi) and PDF,
# so they can be pasted into PowerPoint as they are.

#' Volcano for a slide: colour = higher / lower at p < 0.05, black ring = FDR < fdr, the n_label proteins with the lowest p
#' labelled. r: de_limma() result.
volcano_slide <- function(r, colours, n_label = 20, fdr = 0.05, up = "higher", down = "lower", xlab = "log2 fold change", base = 16) {
  lv <- c(up, down, "p >= 0.05")
  v <- r |> filter(!is.na(logFC), !is.na(p)) |>
    mutate(cls = factor(case_when(p < 0.05 & logFC > 0 ~ up, p < 0.05 & logFC < 0 ~ down, TRUE ~ "p >= 0.05"), lv))
  lab <- v |> slice_min(p, n = n_label, with_ties = FALSE)
  ggplot(v, aes(logFC, -log10(p))) +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey60") +
    geom_point(aes(colour = cls), size = 1.8) +
    geom_point(data = \(d) filter(d, FDR < fdr), shape = 21, size = 3.4, stroke = 0.7, colour = "black") +
    ggrepel::geom_text_repel(data = lab, aes(label = Assay), size = base / 4.2, max.overlaps = Inf, min.segment.length = 0,
                             segment.colour = "grey55", segment.size = 0.3, box.padding = 0.4, seed = 1) +
    scale_colour_manual(values = set_names(c(colours$above, colours$below, "grey78"), lv), name = NULL, drop = FALSE) +
    labs(x = xlab, y = "-log10 p") +
    theme_bw(base_size = base) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot")
}

#' GSEA dot plot for a slide: the n gene sets with the lowest p, ordered by NES; filled = FDR < fdr, open = not.
gsea_slide <- function(g, colours, fdr = 0.05, n = 20, up = "higher", down = "lower", base = 14) {
  if (is.null(g) || !nrow(g)) return(NULL)
  top <- g |> slice_min(p, n = n, with_ties = FALSE) |> arrange(NES) |>
    mutate(label = paste0(str_trunc(if_else(gene_set == toupper(gene_set), str_to_sentence(gene_set), gene_set), 55), " (", collection, ")"),
           label = factor(label, unique(label)), dir = factor(if_else(NES > 0, up, down), c(up, down)),
           sig = factor(if_else(FDR < fdr, sprintf("FDR < %.2g", fdr), sprintf("FDR >= %.2g", fdr)), c(sprintf("FDR < %.2g", fdr), sprintf("FDR >= %.2g", fdr))))
  ggplot(top, aes(NES, label)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_point(aes(colour = dir, size = measured_proteins_in_set, shape = sig), stroke = 1.1) +
    scale_colour_manual(values = set_names(c(colours$above, colours$below), c(up, down)), name = NULL, drop = FALSE) +
    scale_shape_manual(values = set_names(c(16, 1), levels(top$sig)), name = NULL, drop = FALSE) +
    scale_size_area(max_size = 7, name = "proteins") +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 4)), shape = guide_legend(order = 2, override.aes = list(size = 4)),
           size = guide_legend(order = 3)) +                # keys of levels without a gene set stay visible
    labs(x = "normalised enrichment score (NES)", y = NULL) +
    theme_bw(base_size = base) + theme(panel.grid.minor = element_blank(), legend.position = "right", plot.title.position = "plot")
}

#' Levels over the visits for a slide (as the step 13 time course): mean +/- SE at the lesion site (red) and in non-lesional
#' skin (blue), grey band = healthy skin (mean +/- SE), the paired test of each visit written above it.
#' d: Assay, visit_num, site, mean_NPX, se; hb: Assay, mean_NPX, se; pl: Assay, visit_num, label.
levels_slide <- function(d, hb, pl, base = 16) {
  pd <- position_dodge(width = 0.15)
  ggplot(d, aes(visit_num, mean_NPX, colour = site, group = site)) +
    geom_rect(data = hb, aes(xmin = -Inf, xmax = Inf, ymin = mean_NPX - se, ymax = mean_NPX + se), inherit.aes = FALSE, fill = "grey88") +
    geom_hline(data = hb, aes(yintercept = mean_NPX), colour = "grey50", linetype = 2) +
    geom_line(position = pd, linewidth = 0.9) +
    geom_linerange(aes(ymin = mean_NPX - se, ymax = mean_NPX + se), position = pd, linewidth = 0.9, na.rm = TRUE) +
    geom_point(aes(shape = site), position = pd, size = 3) +
    geom_text(data = pl, aes(visit_num, Inf, label = label), inherit.aes = FALSE, vjust = 1.4, size = base / 4.5, colour = "grey20") +
    scale_colour_manual(values = c(`lesion site` = "firebrick", `non-lesional` = "steelblue"), name = NULL) +
    scale_shape_manual(values = c(`lesion site` = 16, `non-lesional` = 17), name = NULL) +
    scale_x_continuous(breaks = sort(unique(d$visit_num)), labels = paste0("V", sort(unique(d$visit_num))), expand = expansion(add = 0.3)) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    facet_wrap(~Assay, scales = "free_y") +
    labs(x = NULL, y = "NPX (dISF)") +
    theme_bw(base_size = base) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot")
}

#' Heatmap for a slide: rows x visits, fill = value (diverging around 0), * = star (e.g. p < 0.05). d: row, col, value, star.
heat_slide <- function(d, colours, legend = "log2FC", base = 14) {
  lim <- max(abs(d$value), na.rm = TRUE)
  ggplot(d, aes(col, row, fill = value)) + geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = if_else(coalesce(star, FALSE), "*", "")), size = base / 2.2, vjust = 0.75) +
    scale_fill_gradient2(low = colours$below, mid = "white", high = colours$above, midpoint = 0, limits = c(-lim, lim), name = legend, na.value = "grey92") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = base) + theme(panel.grid = element_blank(), plot.title.position = "plot")
}

#' NES of gene sets over the visits for a slide, with p and FDR written at each point. d: set, visit_num, NES, p, FDR.
nes_slide <- function(d, colours, base = 16) {
  ggplot(d, aes(visit_num, NES)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_line(colour = "grey50", linewidth = 0.8) +
    geom_point(aes(colour = NES > 0), size = 4) +
    geom_text(aes(label = sprintf("p = %s\nFDR = %s", formatC(p, digits = 2, format = "g"), formatC(FDR, digits = 2, format = "g"))),
              vjust = -0.6, size = base / 4.8, lineheight = 0.9) +
    scale_colour_manual(values = c(`TRUE` = colours$above, `FALSE` = colours$below), labels = c(`TRUE` = "higher at the lesion site", `FALSE` = "lower at the lesion site"), name = NULL) +
    scale_x_continuous(breaks = sort(unique(d$visit_num)), labels = paste0("V", sort(unique(d$visit_num))), expand = expansion(add = 0.4)) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) + facet_wrap(~set) +
    labs(x = NULL, y = "NES (lesion site vs non-lesional)") +
    theme_bw(base_size = base) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot")
}
