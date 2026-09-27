# 07 - Pathway enrichment (GSEA) for every model contrast
# Proteins are ranked by the moderated t statistic; gene sets from MSigDB (msigdbr) plus a
# custom AD/Th2 set. Out: output/enrichment/gsea_results.csv (+ dot plots)

source("R/utils.R")
cfg <- load_config()
clear_outputs(cfg, "enrichment")

files <- list.files(file.path(cfg$paths$output, "models"), pattern = "_results\\.csv$", full.names = TRUE)
if (!length(files)) stop("No model results found - run scripts 04 and 05 first.")
res <- map(files, \(f) read_csv(f, show_col_types = FALSE) |> mutate(file = basename(f))) |> bind_rows()

collections <- cfg$enrichment$collections
sets <- load_gene_sets(collections)
sets[["CUSTOM_AD_TH2_AXIS"]] <- c("CCL17", "CCL22", "CCL18", "CCL26", "CCL11", "CCL13", "CCL24", "IL13",
                                  "IL4", "IL5", "IL31", "IL4R", "IL13RA2", "TSLP", "POSTN")
msg("%d gene sets from %s + custom Th2 set", length(sets), paste(collections, collapse = ", "))

gsea <- res |>
  filter(!is.na(t)) |>
  group_by(file, model, contrast) |>
  group_modify(\(r, key) {
    # Assay names are gene symbols; complexes like "IL12A_IL12B" count for each gene
    st <- r |> transmute(gene = str_split(Assay, "_"), t) |> unnest(gene) |>
      group_by(gene) |> slice_max(abs(t), n = 1, with_ties = FALSE) |> ungroup()
    stats <- setNames(st$t, st$gene)
    out <- suppressWarnings(fgsea::fgseaMultilevel(sets, stats, minSize = min(cfg$enrichment$min_size, 5),
                                                   maxSize = cfg$enrichment$max_size, eps = 0))
    if (!nrow(out)) return(tibble())
    as_tibble(out) |> mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"),
                             small_set = size < cfg$enrichment$min_size)
  }) |> ungroup() |> arrange(file, model, contrast, pval)

save_csv(gsea, cfg, "enrichment", "gsea_results.csv")
msg("%d contrast x gene-set tests; %d with padj < %.2f", nrow(gsea), sum(gsea$padj < cfg$stats$fdr, na.rm = TRUE), cfg$stats$fdr)

top <- gsea |> filter(padj < cfg$stats$fdr) |> group_by(model, contrast) |> slice_min(padj, n = 10, with_ties = FALSE) |> ungroup()
if (nrow(top)) {
  p <- ggplot(top, aes(NES, reorder(str_trunc(pathway, 50), NES), size = size, colour = -log10(padj))) +
    geom_point() + facet_wrap(~paste(model, contrast), scales = "free_y", ncol = 2) +
    labs(title = "Top enriched gene sets", x = "normalised enrichment score", y = NULL)
  save_plot(p, cfg, "enrichment", "top_gene_sets.png", width = 12, height = 2 + 1.2 * n_distinct(paste(top$model, top$contrast)))
}
