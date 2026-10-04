# 09 - Aim 1: descriptive proteome profile of dISF
# Which proteins are reliably measurable in dermal ISF, how this differs by skin state,
# what biology the detectable ISF proteome covers, which proteins vary most (clustered heatmap) and a
# heatmap of the strongest lesional vs non-lesional proteins of step 04.
# NPX is a relative, protein-specific scale, so levels are NOT compared between proteins;
# the profile uses detection rates, variability and pathway composition instead.
# Out: output/isf_profile/*  (incl. isf_profile.xlsx)

source("R/utils.R")
cfg   <- load_config()
clear_outputs(cfg, "isf_profile")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
min_f <- cfg$qc$min_detect_frac

isf <- clean |> filter(matrix == "ISF") |>
  mutate(state_group = case_when(group == "HC" ~ "healthy (HC)",
                                 group == "CPUO" ~ paste("CPUO", state),
                                 TRUE ~ paste("AD", state)))

# ---- detection per protein and skin state --------------------------------------------------------
det_state <- isf |>
  group_by(OlinkID, Assay, UniProt, state_group) |>
  summarise(frac = mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = state_group, values_from = frac, names_prefix = "detected: ")
det_all <- isf |>
  group_by(OlinkID, Assay) |>
  summarise(n_samples = n_distinct(SampleID), frac_detected = mean(!below_lod, na.rm = TRUE),
            sd_npx = sd(value, na.rm = TRUE), .groups = "drop")
ad_cols <- c("detected: AD lesional", "detected: AD non-lesional", "detected: healthy (HC)")
profile <- det_all |> left_join(det_state, by = c("OlinkID", "Assay")) |>
  mutate(
    detection_class = case_when(frac_detected >= 0.9 ~ "robust (>= 90%)",
                                frac_detected >= min_f ~ sprintf("detected (%.0f-90%%)", 100 * min_f),
                                frac_detected > 0.1 ~ "sporadic (10-50%)",
                                TRUE ~ "not detected (<= 10%)")) |>
  arrange(desc(frac_detected))
# detectable in lesional AD skin but not in non-lesional or healthy skin
profile$lesion_restricted <- if (all(ad_cols %in% names(profile)))
  profile[[ad_cols[1]]] >= min_f & profile[[ad_cols[2]]] < min_f & profile[[ad_cols[3]]] < min_f else NA
save_csv(profile, cfg, "isf_profile", "isf_detection_profile.csv")
msg("dISF: %s", paste(names(table(profile$detection_class)), table(profile$detection_class), collapse = ", "))
msg("%d proteins detectable only in lesional AD skin", sum(profile$lesion_restricted, na.rm = TRUE))

by_state <- isf |> group_by(state_group, SampleID) |>
  summarise(n_detected = sum(!below_lod, na.rm = TRUE), .groups = "drop")
p <- ggplot(by_state, aes(reorder(state_group, n_detected, median), n_detected)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  coord_flip() + labs(title = "Proteins above LOD per dISF sample", x = NULL, y = "proteins detected")
save_plot(p, cfg, "isf_profile", "detected_per_sample.png", width = 8, height = 4)

# ---- pathway composition of the detectable ISF proteome (over-representation vs all assays) ----
detected_genes <- profile |> filter(frac_detected >= min_f) |> pull(Assay) |> str_split("_") |> unlist() |> unique()
universe <- clean |> distinct(Assay) |> pull(Assay) |> str_split("_") |> unlist() |> unique()
sets <- map(cfg$enrichment$collections, \(cl) {
  parts <- str_split_fixed(cl, ":", 2)
  g <- if (parts[2] == "") msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1])
       else msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = parts[2])
  split(g$gene_symbol, g$gs_name)
}) |> unlist(recursive = FALSE)
ora <- fgsea::fora(sets, genes = detected_genes, universe = universe,
                   minSize = cfg$enrichment$min_size, maxSize = cfg$enrichment$max_size) |>
  as_tibble() |> mutate(overlapGenes = map_chr(overlapGenes, paste, collapse = ";")) |>
  arrange(pval)
save_csv(ora, cfg, "isf_profile", "isf_detected_pathways.csv")
msg("Detectable dISF proteome: %d of %d panel genes; %d gene sets over-represented (padj < %.2f)",
    length(detected_genes), length(universe), sum(ora$padj < cfg$stats$fdr), cfg$stats$fdr)

# ---- heatmaps (clustered, with sample information bars) --------------------------------------------
# A: the 50 most variable dISF proteins (no group information used to pick them) - which proteins move together,
#    and do samples group by skin state or by patient?
# B: the strongest lesional vs non-lesional proteins of step 04 - the disease signal itself.
# Values: PCNormalizedNPX z-scored per protein over the samples shown (colours capped at +/-3; the workbook has the
# uncapped values). Rows (and in A the columns) are ordered by hierarchical clustering (Euclidean distance, Ward).
m <- wide$ISF
assay_of <- clean |> distinct(OlinkID, Assay) |> (\(d) setNames(d$Assay, d$OlinkID))()
sinfo <- isf |> distinct(SampleID, state_group, SubjectID, group, visit, plate)
ann_cols <- function(ids) {
  a <- sinfo[match(ids, sinfo$SampleID), ]
  data.frame(`skin state` = a$state_group,
             patient = if_else(a$group == "AD", a$SubjectID, if_else(a$group == "HC", "healthy volunteer", "CPUO")),
             visit = coalesce(a$visit, "-"), plate = a$plate, row.names = ids, check.names = FALSE)
}
ann_colours <- function(ann) {
  pal <- \(lv, cols) setNames(rep_len(cols, length(lv)), lv)
  st <- sort(unique(ann$`skin state`))
  pt <- sort(unique(ann$patient)); ad <- setdiff(pt, c("healthy volunteer", "CPUO"))
  list(`skin state` = setNames(c(`AD lesional` = "firebrick", `AD ex-lesional` = "orange", `AD non-lesional` = "steelblue",
                                 `healthy (HC)` = "darkgreen", `CPUO lesional` = "purple", `CPUO non-lesional` = "plum")[st], st),
       patient = c(setNames(grDevices::hcl.colors(max(length(ad), 2), "Dark 3")[seq_along(ad)], ad),
                   `healthy volunteer` = "grey80", CPUO = "grey45")[pt],
       visit = setNames(grDevices::hcl.colors(length(unique(ann$visit)), "Blues 3", rev = TRUE), sort(unique(ann$visit))),
       plate = setNames(grDevices::hcl.colors(length(unique(ann$plate)), "Set 2"), sort(unique(ann$plate))))
}
zscore <- \(x) t(scale(t(x)))
draw_heatmap <- function(z, ann_col, file, title, cluster_cols = TRUE, gaps_col = NULL, ann_row = NULL, ann_row_colours = NULL) {
  zc <- pmax(pmin(z, 3), -3); zc[is.na(zc)] <- 0                        # display / clustering only
  rownames(zc) <- make.unique(assay_of[rownames(z)] %|% rownames(z))
  if (!is.null(ann_row)) rownames(ann_row) <- rownames(zc)
  cols <- ann_colours(ann_col)
  ph <- pheatmap::pheatmap(zc, color = grDevices::colorRampPalette(c("steelblue", "white", "firebrick"))(101),
                           breaks = seq(-3, 3, length.out = 102), clustering_method = "ward.D2",
                           cluster_rows = nrow(zc) > 2, cluster_cols = cluster_cols && ncol(zc) > 2, gaps_col = gaps_col,
                           annotation_col = ann_col, annotation_row = ann_row,
                           annotation_colors = c(cols, ann_row_colours), show_colnames = FALSE, fontsize_row = 7,
                           border_color = NA, main = title, silent = TRUE)
  f <- out_path(cfg, "isf_profile", file)
  grDevices::png(f, width = 14, height = 10, units = "in", res = 150); grid::grid.draw(ph$gtable); grDevices::dev.off()
  list(rows = if (inherits(ph$tree_row, "hclust")) rownames(z)[ph$tree_row$order] else rownames(z),
       cols = if (inherits(ph$tree_col, "hclust")) colnames(z)[ph$tree_col$order] else colnames(z))
}
`%|%` <- function(a, b) ifelse(is.na(a), b, a)
heat_tabs <- list()
export_heat <- function(prefix, z, raw, ord, extra = NULL) {
  ann <- tibble(OlinkID = ord$rows, Assay = assay_of[ord$rows], row = seq_along(ord$rows)) |>
    (\(d) if (is.null(extra)) d else left_join(d, extra, by = "OlinkID"))()
  heat_tabs[[paste0(prefix, "_samples")]] <<- tibble(column = seq_along(ord$cols), SampleID = ord$cols) |>
    left_join(sinfo |> select(SampleID, skin_state = state_group, SubjectID, visit, plate), by = "SampleID")
  heat_tabs[[paste0(prefix, "_zscores")]] <<- bind_cols(ann, as_tibble(z[ord$rows, ord$cols, drop = FALSE]))
  heat_tabs[[paste0(prefix, "_NPX")]] <<- bind_cols(ann, as_tibble(raw[ord$rows, ord$cols, drop = FALSE]))
}

if (!is.null(m) && nrow(m) > 2 && ncol(m) > 2) {
  sds <- apply(m, 1, sd, na.rm = TRUE)
  top <- names(sort(sds, decreasing = TRUE))[seq_len(min(50, nrow(m)))]
  z <- zscore(m[top, , drop = FALSE])
  ord <- draw_heatmap(z, ann_cols(colnames(z)), "top_variable_heatmap.png",
                      "50 most variable dISF proteins (z-score per protein; rows and samples clustered)")
  export_heat("variable", z, m, ord, tibble(OlinkID = top, rank_by_SD = seq_along(top), SD_NPX = sds[top]))
}

res_path <- file.path(cfg$paths$output, "models", "ISF_results.csv")
if (!is.null(m) && file.exists(res_path)) {
  lnl <- read_csv(res_path, show_col_types = FALSE) |>
    filter(model == "states_all_visits", contrast == "AD_L_vs_NL", OlinkID %in% rownames(m)) |> arrange(P.Value)
  n_each <- 25
  sel <- bind_rows(lnl |> filter(logFC > 0) |> slice_head(n = n_each), lnl |> filter(logFC < 0) |> slice_head(n = n_each))
  ids <- sinfo |> filter(state_group %in% c("AD lesional", "AD ex-lesional", "AD non-lesional", "healthy (HC)")) |>
    mutate(o = match(state_group, c("AD lesional", "AD ex-lesional", "AD non-lesional", "healthy (HC)"))) |>
    arrange(o, SubjectID, visit) |> filter(SampleID %in% colnames(m))
  if (nrow(sel) >= 3 && nrow(ids) > 2) {
    z <- zscore(m[sel$OlinkID, ids$SampleID, drop = FALSE])
    ann_row <- data.frame(direction = if_else(sel$logFC > 0, "higher in lesional", "lower in lesional"),
                          significant = if_else(sel$adj.P.Val < cfg$stats$fdr, sprintf("FDR < %g", cfg$stats$fdr), "not significant"))
    row_cols <- list(direction = c(`higher in lesional` = "firebrick", `lower in lesional` = "steelblue"),
                     significant = setNames(c("black", "grey85"), c(sprintf("FDR < %g", cfg$stats$fdr), "not significant")))
    row_cols$significant <- row_cols$significant[unique(ann_row$significant)]
    row_cols$direction <- row_cols$direction[unique(ann_row$direction)]
    ord <- draw_heatmap(z, ann_cols(colnames(z)), "lesional_vs_nonlesional_heatmap.png",
                        sprintf("Strongest lesional vs non-lesional dISF proteins (step 04, all visits): top %d up and %d down by p-value",
                                sum(sel$logFC > 0), sum(sel$logFC < 0)),
                        cluster_cols = FALSE, gaps_col = cumsum(table(factor(ids$o, levels = 1:4)))[1:3], ann_row = ann_row,
                        ann_row_colours = row_cols)
    export_heat("lesional", z, m, ord, sel |> transmute(OlinkID, logFC_L_vs_NL = logFC, P.Value, FDR = adj.P.Val))
  }
} else msg("Step 04 results not found - lesional vs non-lesional heatmap skipped.")

writexl::write_xlsx(c(list(detection_profile = profile, detected_pathways = ora), heat_tabs),
                    out_path(cfg, "isf_profile", "isf_profile.xlsx"))
