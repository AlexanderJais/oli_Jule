# Helpers for the relapse follow-up (run_relapse_followup.R): checks of the relapse vs non-relapse comparison, marker
# and score comparisons with an exact relabelling p, cytokine heatmaps over the visits, volcano with every p < 0.05 labelled.

#' OlinkIDs of the assays of a gene (complexes such as "IL12A_IL12B" count for each gene); the first per gene if several.
assay_ids <- function(info, genes) {
  hit <- info |> mutate(gene = str_split(Assay, "_")) |> unnest(gene) |> filter(gene %in% genes) |>
    group_by(gene) |> slice_head(n = 1) |> ungroup()
  hit[match(intersect(genes, hit$gene), hit$gene), ] |> select(gene, OlinkID, Assay) |> distinct(OlinkID, .keep_all = TRUE)
}

#' Exact p of the difference in means: share of all ways to choose the relapse group with a difference at least as large
#' (two-sided: in absolute value; "greater": higher in relapse). Patients without a value are left out.
exact_p_diff <- function(v, is_rel, alternative = c("two.sided", "greater")) {
  alternative <- match.arg(alternative)
  ok <- !is.na(v); v <- v[ok]; y <- is_rel[ok]
  if (sum(y) < 2 || sum(!y) < 2) return(NA_real_)
  obs <- mean(v[y]) - mean(v[!y])
  d <- apply(utils::combn(length(v), sum(y)), 2, \(i) mean(v[i]) - mean(v[-i]))
  if (alternative == "greater") mean(d >= obs - 1e-9) else mean(abs(d) >= abs(obs) - 1e-9)
}

#' Relapse vs non-relapse for one per-patient value (a protein or a score, log2 scale): group means, log2FC
#' (relapse - non-relapse), Welch t-test, exact relabelling p, AUC, and whether the difference keeps its sign when
#' any one patient is left out (FALSE = one patient drives it).
compare_groups <- function(v, is_rel) {
  a <- v[is_rel & !is.na(v)]; b <- v[!is_rel & !is.na(v)]
  if (length(a) < 2 || length(b) < 2)
    return(tibble(n_relapse = length(a), n_nonrelapse = length(b), mean_relapse = NA_real_, mean_nonrelapse = NA_real_,
                  log2FC = NA_real_, p_t = NA_real_, p_exact = NA_real_, AUC = NA_real_, same_sign_without_any_one_patient = NA))
  d <- mean(a) - mean(b)
  loo <- map_dbl(seq_along(v), \(i) { x <- v[-i]; y <- is_rel[-i]; mean(x[y], na.rm = TRUE) - mean(x[!y], na.rm = TRUE) })
  tibble(n_relapse = length(a), n_nonrelapse = length(b), mean_relapse = mean(a), mean_nonrelapse = mean(b), log2FC = d,
         p_t = tryCatch(t.test(a, b)$p.value, error = \(e) NA_real_), p_exact = exact_p_diff(v, is_rel),
         AUC = mean(outer(a, b, ">") + 0.5 * outer(a, b, "==")),
         same_sign_without_any_one_patient = all(sign(loo) == sign(d), na.rm = TRUE))
}

#' GSEA result with the original gene set key (pathway), so sets can be re-tested; same names as run_gsea().
gsea_keyed <- function(res, sets, min_size = 10, max_size = 500) {
  g <- run_gsea(res, sets, min_size, max_size)
  if (!nrow(g)) return(g)
  key <- tibble(pathway = names(sets), collection = str_extract(pathway, "^[A-Z]+"),
                gene_set = if_else(pathway == toupper(pathway), str_replace_all(str_remove(pathway, "^[A-Z]+_"), "_", " "),
                                   str_remove(pathway, "^KEGG_"))) |> distinct(collection, gene_set, .keep_all = TRUE)
  g |> left_join(key, by = c("collection", "gene_set"))
}

#' Mean log2FC of the measured members of each gene set (does the NES sign agree with the proteins in the set?).
set_mean_logfc <- function(res, sets, pathways) {
  fc <- res |> filter(!is.na(logFC)) |> transmute(gene = str_split(Assay, "_"), logFC) |> unnest(gene) |>
    group_by(gene) |> summarise(logFC = mean(logFC), .groups = "drop")
  map_dbl(pathways, \(p) { x <- if (is.na(p)) numeric() else fc$logFC[fc$gene %in% sets[[p]]]; if (length(x)) mean(x) else NA_real_ })
}

#' Volcano plot with every protein below p_label labelled (colour: higher / lower in relapse at p < p_label; ring:
#' FDR < fdr). Returns list(plot, data) - data is exactly what is drawn.
volcano_all_p <- function(r, title, colours, p_label = 0.05, fdr = 0.05, gsea_used = TRUE) {
  v <- r |> mutate(minus_log10_p = -log10(p),
                   class = case_when(p < p_label & logFC > 0 ~ sprintf("higher in relapse (p < %.2g)", p_label),
                                     p < p_label & logFC < 0 ~ sprintf("lower in relapse (p < %.2g)", p_label), TRUE ~ "not significant"),
                   labelled = coalesce(p < p_label, FALSE), FDR_significant = coalesce(FDR < fdr, FALSE))
  lab <- v |> filter(labelled)
  pal <- set_names(c(colours$above, colours$below, "grey75"), c(sprintf("higher in relapse (p < %.2g)", p_label),
                                                                 sprintf("lower in relapse (p < %.2g)", p_label), "not significant"))
  n_lab <- nrow(lab)
  p <- ggplot(v, aes(logFC, minus_log10_p)) +
    geom_hline(yintercept = -log10(p_label), linetype = 2, colour = "grey60", linewidth = 0.3) +
    geom_point(aes(colour = class), size = 1.1) +
    geom_point(data = \(d) filter(d, FDR_significant), shape = 21, size = 2.4, colour = "black", stroke = 0.4) +
    ggrepel::geom_text_repel(data = lab, aes(label = Assay), size = if (n_lab > 150) 1.6 else if (n_lab > 60) 2 else 2.6,
                             max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60", segment.size = 0.2,
                             max.time = 5, max.iter = 1e5, seed = 1) +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "log2 fold change (relapse - non-relapse)", y = "-log10 p", title = title,
         subtitle = sprintf("%d proteins tested; %d with p < %.2g labelled (%d higher, %d lower in relapse); black ring: FDR < %.2f (%d). Dashed line: p = %.2g",
                            nrow(v), n_lab, p_label, sum(lab$logFC > 0), sum(lab$logFC < 0), fdr, sum(v$FDR_significant, na.rm = TRUE), p_label)) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot", legend.position = "bottom")
  list(plot = p, data = v |> transmute(Assay, OlinkID, UniProt, log2FC = logFC, p, minus_log10_p, FDR, FDR_significant, labelled, class,
                                      t_used_for_GSEA = t, across(starts_with("mean_NPX_")), across(starts_with("pct_above_LOD_"))) |>
         arrange(p))
}

#' Group-mean heatmap of proteins over the visits: y = protein, x = visit, one panel per group (and matrix).
#' Fill = mean NPX minus the protein's mean over all cells of its matrix; text = mean NPX and n; * = most values below LOD.
heatmap_group_means <- function(cm, colours) {
  ggplot(cm, aes(visit, Assay, fill = centred)) + geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f%s\nn=%d", mean_NPX, if_else(pct_above_LOD < 50, "*", ""), n)), size = 2.6, lineheight = 0.9) +
    scale_fill_gradient2(low = colours$below, mid = "white", high = colours$above, midpoint = 0, name = "NPX - protein mean\n(log2)") +
    facet_grid(matrix ~ group, scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), plot.title.position = "plot", strip.text = element_text(face = "bold"))
}
