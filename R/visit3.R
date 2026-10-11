# Helpers for the dedicated visit 3 analysis (run_visit3.R): differential proteins per matrix, volcano plots, GSEA.

#' limma for one comparison of two sample groups.
#' d: long Olink data (SampleID, OlinkID, Assay, UniProt, value, below_lod); samples: SampleID, SubjectID, grp (factor:
#' reference level first). A protein is tested if it is above LOD in >= min_f of the samples of at least one group.
#' paired = TRUE adds the person as a factor (both samples of one person are compared with each other).
#' logFC = second group minus first group (log2).
de_limma <- function(d, samples, paired = FALSE, min_f = 0.5, fdr = 0.05) {
  lv <- levels(samples$grp)
  d <- d |> filter(SampleID %in% samples$SampleID) |> left_join(samples |> select(SampleID, grp), by = "SampleID")
  grp_stats <- d |> group_by(OlinkID, grp) |>
    summarise(mean = mean(value, na.rm = TRUE), pct = 100 * mean(!below_lod, na.rm = TRUE), .groups = "drop")
  tested <- grp_stats |> group_by(OlinkID) |> summarise(k = any(pct >= 100 * min_f, na.rm = TRUE)) |> filter(k) |> pull(OlinkID)
  w <- d |> filter(OlinkID %in% .env$tested) |> select(OlinkID, SampleID, value) |>   # .env: npx_clean has a column 'keep'
    pivot_wider(names_from = SampleID, values_from = value, values_fn = \(v) mean(v, na.rm = TRUE))
  M <- as.matrix(w[, -1]); rownames(M) <- w$OlinkID
  s <- samples[match(colnames(M), samples$SampleID), ]
  design <- if (paired) model.matrix(~ factor(SubjectID) + grp, s) else model.matrix(~ grp, s)
  fit <- limma::eBayes(limma::lmFit(M, design))
  tt <- limma::topTable(fit, coef = ncol(design), number = Inf, sort.by = "none")
  wide <- grp_stats |> mutate(grp = as.character(grp)) |>
    pivot_wider(names_from = grp, values_from = c(mean, pct), names_glue = "{.value}_{grp}")
  info <- d |> distinct(OlinkID, Assay, UniProt)
  tibble(OlinkID = rownames(tt), logFC = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val) |>
    left_join(info, by = "OlinkID") |> left_join(wide, by = "OlinkID") |>
    mutate(significant = FDR < fdr, direction = case_when(!significant ~ "not significant", logFC > 0 ~ paste("higher in", lv[2]),
                                                          TRUE ~ paste("higher in", lv[1]))) |>
    rename_with(\(x) str_replace(x, "^mean_", "mean_NPX_"), starts_with("mean_")) |>
    rename_with(\(x) str_replace(x, "^pct_", "pct_above_LOD_"), starts_with("pct_")) |>
    relocate(Assay, UniProt, .after = OlinkID) |> relocate(significant, direction, .after = FDR) |>
    mutate(across(where(is.double), \(v) signif(v, 4))) |> arrange(p)
}

#' Volcano plot; the n_label most regulated proteins are labelled (significant ones by |logFC|, then the lowest p).
volcano_plot <- function(res, title, subtitle, colours, n_label = 25) {
  res <- res |> mutate(cl = case_when(significant & logFC > 0 ~ "up", significant ~ "down", TRUE ~ "ns"))
  lab <- res |> arrange(desc(significant), if_else(significant, -abs(logFC), p)) |> slice_head(n = n_label)
  ggplot(res, aes(logFC, -log10(p))) +
    geom_point(aes(colour = cl), size = 1.1) +
    ggrepel::geom_text_repel(data = lab, aes(label = Assay), size = 2.6, max.overlaps = Inf, min.segment.length = 0,
                             segment.colour = "grey60", segment.size = 0.3) +
    scale_colour_manual(values = c(up = colours$above, down = colours$below, ns = "grey75"), guide = "none") +
    labs(x = "log2 fold change", y = "-log10 p", title = title, subtitle = subtitle) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
}

#' Gene sets: MSigDB collections (config enrichment$collections) plus KEGG (current file if downloaded, else MSigDB legacy).
v3_gene_sets <- function(cfg) {
  ms <- map(cfg$enrichment$collections %||% c("H", "C2:CP:REACTOME", "C5:GO:BP"), \(cl) {
    parts <- str_split_fixed(cl, ":", 2)
    g <- if (parts[2] == "") msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1])
         else msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = parts[2])
    split(g$gene_symbol, g$gs_name)
  }) |> unlist(recursive = FALSE)
  kegg <- read_kegg(cfg$paths$kegg)
  if (is.null(kegg)) {
    k <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
    kegg <- split(k$gene_symbol, k$gs_name)
  }
  c(ms, kegg)
}

#' GSEA (fgsea) on the limma t statistics; complexes like "IL12A_IL12B" count for each gene.
run_gsea <- function(res, sets, min_size = 10, max_size = 500) {
  st <- res |> filter(!is.na(t)) |> transmute(gene = str_split(Assay, "_"), t) |> unnest(gene) |>
    group_by(gene) |> slice_max(abs(t), n = 1, with_ties = FALSE) |> ungroup()
  out <- suppressWarnings(fgsea::fgseaMultilevel(sets, setNames(st$t, st$gene), minSize = min_size, maxSize = max_size, eps = 0))
  if (!nrow(out)) return(tibble())
  as_tibble(out) |>
    transmute(collection = str_extract(pathway, "^[A-Z]+"),
              gene_set = if_else(pathway == toupper(pathway), str_replace_all(str_remove(pathway, "^[A-Z]+_"), "_", " "), str_remove(pathway, "^KEGG_")),
              measured_proteins_in_set = size, NES = round(NES, 2), p = signif(pval, 4), FDR = signif(padj, 4),
              leading_edge = map_chr(leadingEdge, paste, collapse = ", ")) |> arrange(p)
}

#' Dot plot of the top gene sets of one GSEA.
gsea_plot <- function(g, title, colours, fdr = 0.05, n = 25) {
  sig <- g |> filter(FDR < fdr)
  top <- (if (nrow(sig)) sig else g) |> slice_head(n = n) |>
    mutate(label = str_trunc(paste0(if_else(gene_set == toupper(gene_set), str_to_sentence(gene_set), gene_set), " (", collection, ")"), 75),
           label = factor(label, rev(unique(label))))
  ggplot(top, aes(NES, label)) + geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(aes(size = measured_proteins_in_set, colour = NES > 0)) +
    scale_colour_manual(values = c(`TRUE` = colours$above, `FALSE` = colours$below), labels = c(`TRUE` = "up", `FALSE` = "down"), name = NULL) +
    scale_size_area(max_size = 6, name = "measured\nproteins") +
    labs(x = "normalised enrichment score (NES)", y = NULL,
         title = if (nrow(sig)) sprintf("%s: %d gene sets with FDR < %.2f (top %d)", title, nrow(sig), fdr, nrow(top))
                 else sprintf("%s: no gene set with FDR < %.2f - the %d lowest p values", title, fdr, nrow(top))) +
    theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
}
