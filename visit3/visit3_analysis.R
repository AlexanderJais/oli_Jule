# Dedicated analysis of one visit (default V3; run with run_visit3.R), MicroAD samples only.
#   dISF: lesion site vs healthy skin, lesion site vs non-lesional skin (same patients, paired)
#   serum: AD vs healthy (serum has no skin site, so there is no lesional vs non-lesional comparison in serum)
# Differential proteins (limma) per matrix, volcano plots, GSEA, and what dISF shows that serum does not.
# In:  output of steps 01-02 (data/npx_clean.rds); gene sets from MSigDB and KEGG (tools/download_kegg.R, optional)
# Out: output/visit3/ - V3_DEGs_dISF.xlsx, V3_DEGs_serum.xlsx, V3_GSEA.xlsx, V3_dISF_vs_serum.xlsx, volcano/, gsea/,
#      V3_dISF_vs_serum.pdf

source("R/utils.R")
source("R/qc_overview.R")
source("R/visit3.R")
cfg <- load_config()
v3 <- cfg$visit3 %||% list()
vis <- v3$visit %||% "V3"; n_lab <- v3$top_labels %||% 25
fdr <- cfg$stats$fdr %||% 0.05; min_f <- cfg$qc$min_detect_frac %||% 0.5
colours <- cfg$qc_overview$colours %||% list(above = "#69005F", below = "#FF506E")
out <- "visit3"
clear_outputs(cfg, out)
set.seed(1)

clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
smp <- clean |> filter(cohort == "MicroAD") |> distinct(SampleID, SubjectID, matrix, group, visit, site, state)
grp <- list(
  isf_L  = smp |> filter(matrix == "ISF", group == "AD", visit == vis, site == "L"),
  isf_NL = smp |> filter(matrix == "ISF", group == "AD", visit == vis, site == "NL"),
  isf_H  = smp |> filter(matrix == "ISF", group == "HC"),
  ser_AD = smp |> filter(matrix == "Serum", group == "AD", visit == vis),
  ser_H  = smp |> filter(matrix == "Serum", group == "HC"))
small <- names(grp)[map_int(grp, nrow) < 3]
if (length(small)) stop(sprintf("Too few samples at %s for: %s (need >= 3 each).", vis, paste(small, collapse = ", ")), call. = FALSE)
l_states <- grp$isf_L |> count(state) |> mutate(t = paste(n, state)) |> pull(t) |> paste(collapse = ", ")
h_visits <- paste(sort(unique(grp$isf_H$visit)), collapse = ", ")
samples_txt <- c(
  sprintf("dISF lesion site at %s: %d samples of %d patients (clinical state of the lesion site: %s)", vis, nrow(grp$isf_L), n_distinct(grp$isf_L$SubjectID), l_states),
  sprintf("dISF non-lesional skin at %s: %d samples", vis, nrow(grp$isf_NL)),
  sprintf("dISF healthy skin: %d samples of healthy volunteers (visit %s - healthy volunteers are sampled once)", nrow(grp$isf_H), h_visits),
  sprintf("serum AD at %s: %d samples; serum healthy: %d samples (visit %s)", vis, nrow(grp$ser_AD), nrow(grp$ser_H),
          paste(sort(unique(grp$ser_H$visit)), collapse = ", ")))
walk(samples_txt, \(x) msg(x))

# ---- differential proteins ------------------------------------------------------------------------------------------
mk <- function(a, b, la, lb) bind_rows(a |> mutate(grp = la), b |> mutate(grp = lb)) |> mutate(grp = factor(grp, c(la, lb)))
cmp <- list(
  dISF_healthy_vs_lesional     = list(matrix = "dISF", sheet = "healthy_vs_lesional", s = mk(grp$isf_H, grp$isf_L, "healthy", "lesional"), paired = FALSE,
                                      title = "dISF: lesion site vs healthy skin", fc = "lesional - healthy"),
  dISF_lesional_vs_nonlesional = list(matrix = "dISF", sheet = "lesional_vs_nonlesional", s = mk(grp$isf_NL, grp$isf_L, "non-lesional", "lesional"), paired = TRUE,
                                      title = "dISF: lesion site vs non-lesional skin (same patients)", fc = "lesional - non-lesional"),
  serum_healthy_vs_AD          = list(matrix = "serum", sheet = "healthy_vs_AD", s = mk(grp$ser_H, grp$ser_AD, "healthy", "AD"), paired = FALSE,
                                      title = "Serum: AD vs healthy", fc = "AD - healthy"))
res <- map(cmp, \(k) de_limma(clean, k$s, paired = k$paired, min_f = min_f, fdr = fdr))
summ <- imap(cmp, \(k, nm) {
  r <- res[[nm]]; lv <- levels(k$s$grp)
  tibble(matrix = k$matrix, comparison = k$sheet, log2FC = k$fc, n_samples = sprintf("%s %d, %s %d", lv[1], sum(k$s$grp == lv[1]), lv[2], sum(k$s$grp == lv[2])),
         proteins_tested = nrow(r), significant = sum(r$significant), higher = sum(r$significant & r$logFC > 0), lower = sum(r$significant & r$logFC < 0))
}) |> bind_rows()
print(summ)
deg_readme <- function(m) tibble(item = c("samples", "comparisons", "test", "columns"),
  note = c(paste(samples_txt, collapse = " | "),
           paste(sprintf("%s: log2FC = %s", summ$comparison[summ$matrix == m], summ$log2FC[summ$matrix == m]), collapse = "; "),
           sprintf("limma (PCNormalizedNPX); proteins above LOD in >= %.0f%% of the samples of at least one of the two groups; significant = FDR (Benjamini-Hochberg) < %.2f%s",
                   100 * min_f, fdr, if (m == "dISF") "; lesional vs non-lesional is paired (person as factor)" else ""),
           "mean_NPX_<group> and pct_above_LOD_<group> per group; direction = in which group the protein is higher"))
xl_list <- function(m) c(list(README = deg_readme(m), summary = summ |> filter(matrix == m)),
                         set_names(res[names(cmp)[map_chr(cmp, "matrix") == m]], map_chr(keep(cmp, \(k) k$matrix == m), "sheet")))
writexl::write_xlsx(xl_list("dISF"), out_path(cfg, out, paste0(vis, "_DEGs_dISF.xlsx")))
writexl::write_xlsx(xl_list("serum"), out_path(cfg, out, paste0(vis, "_DEGs_serum.xlsx")))

# ---- figures: volcano plots ---------------------------------------------------------------------------------------------
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
save_fig <- function(p, file, width = 11.69, height = 8.27) for (type in c("pdf", "png")) {
  open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
  print(set_family(p, fi$family))
  close_device(fi)
}
iwalk(cmp, \(k, nm) {
  r <- res[[nm]]; s <- summ[summ$comparison == k$sheet, ]
  save_fig(volcano_plot(r, sprintf("%s, %s: %d of %d proteins differ (FDR < %.2f; %d higher, %d lower)", k$title, vis, s$significant, s$proteins_tested, fdr, s$higher, s$lower),
                        sprintf("%s. log2FC = %s. Labelled: the %d most regulated proteins (significant ones by fold change, then lowest p)", s$n_samples, k$fc, n_lab),
                        colours, n_lab), file.path("volcano", paste0(vis, "_volcano_", nm)))
})

# ---- GSEA ---------------------------------------------------------------------------------------------------------------
sets <- v3_gene_sets(cfg)
kegg_src <- if (file.exists(cfg$paths$kegg %||% "")) attr(read_kegg(cfg$paths$kegg), "release") else "KEGG legacy (MSigDB)"
gsea <- map(res, \(r) run_gsea(r, sets, min_size = cfg$enrichment$min_size %||% 10, max_size = cfg$enrichment$max_size %||% 500))
walk2(names(gsea), gsea, \(nm, g) msg("GSEA %s: %d gene sets tested, %d with FDR < %.2f", nm, nrow(g), sum(g$FDR < fdr), fdr))
writexl::write_xlsx(c(list(README = tibble(item = c("method", "gene sets", "sheets"), note = c(
  "GSEA (fgsea) on the limma t statistics of each comparison: all tested proteins ranked; NES > 0 = gene set higher in the second group (lesional / AD)",
  sprintf("MSigDB %s and KEGG (%s); sets with %d-%d measured proteins", paste(cfg$enrichment$collections, collapse = ", "), kegg_src,
          cfg$enrichment$min_size %||% 10, cfg$enrichment$max_size %||% 500),
  "one sheet per comparison"))), gsea), out_path(cfg, out, paste0(vis, "_GSEA.xlsx")))
iwalk(gsea, \(g, nm) if (nrow(g)) save_fig(gsea_plot(g, paste0(cmp[[nm]]$title, ", ", vis), colours, fdr), file.path("gsea", paste0(vis, "_GSEA_", nm))))

# ---- what does dISF show that serum does not? (AD vs healthy in both matrices) -----------------------------------------
d <- res$dISF_healthy_vs_lesional; s <- res$serum_healthy_vs_AD
prot <- full_join(d |> select(OlinkID, Assay_d = Assay, logFC_dISF = logFC, FDR_dISF = FDR, sig_dISF = significant),
                  s |> select(OlinkID, Assay_s = Assay, logFC_serum = logFC, FDR_serum = FDR, sig_serum = significant), by = "OlinkID") |>
  mutate(Assay = coalesce(Assay_d, Assay_s), measurable_dISF = !is.na(FDR_dISF), measurable_serum = !is.na(FDR_serum),
         category = case_when(sig_dISF %in% TRUE & sig_serum %in% TRUE & sign(logFC_dISF) == sign(logFC_serum) ~ "changed in both (same direction)",
                              sig_dISF %in% TRUE & sig_serum %in% TRUE ~ "changed in both (opposite direction)",
                              sig_dISF %in% TRUE & !measurable_serum ~ "changed in dISF only - not measurable in serum",
                              sig_dISF %in% TRUE ~ "changed in dISF only - measurable in serum",
                              sig_serum %in% TRUE & !measurable_dISF ~ "changed in serum only - not measurable in dISF",
                              sig_serum %in% TRUE ~ "changed in serum only - measurable in dISF",
                              TRUE ~ "not changed")) |>
  select(OlinkID, Assay, category, logFC_dISF, FDR_dISF, logFC_serum, FDR_serum, measurable_dISF, measurable_serum) |>
  arrange(category, pmin(FDR_dISF, FDR_serum, na.rm = TRUE))
gd <- gsea$dISF_healthy_vs_lesional; gs <- gsea$serum_healthy_vs_AD
gsets <- full_join(gd |> select(collection, gene_set, NES_dISF = NES, FDR_dISF = FDR, measured_in_dISF = measured_proteins_in_set),
                   gs |> select(collection, gene_set, NES_serum = NES, FDR_serum = FDR, measured_in_serum = measured_proteins_in_set),
                   by = c("collection", "gene_set")) |>
  mutate(category = case_when(FDR_dISF < fdr & FDR_serum < fdr ~ "significant in both", FDR_dISF < fdr ~ "significant in dISF only",
                              FDR_serum < fdr ~ "significant in serum only", TRUE ~ "not significant"),
         extra_proteins_measured_in_dISF = coalesce(measured_in_dISF, 0) - coalesce(measured_in_serum, 0)) |>
  arrange(factor(category, c("significant in dISF only", "significant in both", "significant in serum only", "not significant")),
          pmin(FDR_dISF, FDR_serum, na.rm = TRUE))
pc <- count(prot, category)
gc <- count(gsets, category)
nn <- \(tb, x) sum(tb$n[tb$category %in% x])
key <- tibble(question = c(
  "proteins measurable (tested)", "proteins changed, lesion site vs healthy (dISF) / AD vs healthy (serum)",
  "changed in both matrices", "changed only in dISF", "   of these: not measurable in serum at all", "changed only in serum",
  "gene sets significant (GSEA)", "gene sets significant in both", "gene sets significant only in dISF", "gene sets significant only in serum"),
  dISF = c(nrow(d), sum(d$significant), nn(pc, c("changed in both (same direction)", "changed in both (opposite direction)")),
           nn(pc, c("changed in dISF only - not measurable in serum", "changed in dISF only - measurable in serum")),
           nn(pc, "changed in dISF only - not measurable in serum"), NA, sum(gd$FDR < fdr), nn(gc, "significant in both"),
           nn(gc, "significant in dISF only"), NA),
  serum = c(nrow(s), sum(s$significant), nn(pc, c("changed in both (same direction)", "changed in both (opposite direction)")), NA, NA,
            nn(pc, c("changed in serum only - not measurable in dISF", "changed in serum only - measurable in dISF")),
            sum(gs$FDR < fdr), nn(gc, "significant in both"), NA, nn(gc, "significant in serum only")))
print(key, n = Inf)
writexl::write_xlsx(list(README = tibble(item = c("comparison", "samples", "proteins", "gene_sets", "key_numbers"), note = c(
  sprintf("AD vs healthy at %s in both matrices: dISF lesion site vs healthy skin, serum AD vs healthy (FDR < %.2f)", vis, fdr),
  paste(samples_txt, collapse = " | "),
  "every protein: changed in which matrix, and whether it is measurable in the other (above LOD in >= 50% of at least one group)",
  "every gene set: GSEA in both matrices; measured_in_* = proteins of the set measured in that matrix",
  "the counts behind the comparison")),
  key_numbers = key, proteins = prot, gene_sets = gsets), out_path(cfg, out, paste0(vis, "_dISF_vs_serum.xlsx")))

p1 <- ggplot(pc |> filter(category != "not changed") |> mutate(category = factor(category, rev(sort(unique(category))))), aes(n, category)) +
  geom_col(aes(fill = str_detect(category, "dISF only")), width = 0.7) + geom_text(aes(label = n), hjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = colours$above, `FALSE` = "grey65"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "proteins", y = NULL, title = sprintf("AD vs healthy at %s: which proteins change in which matrix?", vis),
       subtitle = sprintf("dISF: lesion site vs healthy skin (%d of %d tested proteins changed); serum: AD vs healthy (%d of %d). FDR < %.2f",
                          sum(d$significant), nrow(d), sum(s$significant), nrow(s), fdr)) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
both <- prot |> filter(measurable_dISF, measurable_serum)
p2 <- ggplot(both, aes(logFC_serum, logFC_dISF)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = case_when(str_detect(category, "both") ~ "both", str_detect(category, "dISF only") ~ "dISF only",
                                    str_detect(category, "serum only") ~ "serum only", TRUE ~ "not changed")), size = 1) +
  scale_colour_manual(values = c(both = "black", `dISF only` = colours$above, `serum only` = colours$below, `not changed` = "grey80"), name = "changed in") +
  labs(x = "log2FC serum (AD - healthy)", y = "log2FC dISF (lesion site - healthy)",
       title = sprintf("Proteins measurable in both matrices (%d): fold change in dISF vs serum (r = %.2f)", nrow(both),
                       suppressWarnings(cor(both$logFC_dISF, both$logFC_serum, use = "complete.obs")))) +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
top_sets <- gsets |> filter(category != "not significant") |> slice_head(n = 30) |>
  mutate(label = str_trunc(paste0(if_else(gene_set == toupper(gene_set), str_to_sentence(gene_set), gene_set), " (", collection, ")"), 70),
         label = factor(label, rev(unique(label))))
p3 <- if (nrow(top_sets)) ggplot(top_sets |> select(label, category, NES_dISF, NES_serum, FDR_dISF, FDR_serum) |>
                                    pivot_longer(c(NES_dISF, NES_serum), names_to = "m", values_to = "NES") |>
                                    mutate(matrix = if_else(m == "NES_dISF", "dISF", "serum"), sig = if_else(matrix == "dISF", FDR_dISF, FDR_serum) < fdr),
                                  aes(NES, label, colour = matrix, shape = coalesce(sig, FALSE))) +
  geom_vline(xintercept = 0, colour = "grey70") + geom_point(size = 2.4) +
  scale_colour_manual(values = c(dISF = colours$above, serum = colours$below)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = sprintf("FDR < %.2f", fdr), `FALSE` = "not significant"), name = NULL) +
  labs(x = "NES (> 0: higher in AD / lesion site)", y = NULL,
       title = sprintf("Gene sets: %d significant only in dISF, %d in both, %d only in serum (top %d shown)",
                       nn(gc, "significant in dISF only"), nn(gc, "significant in both"), nn(gc, "significant in serum only"), nrow(top_sets))) +
  theme_bw(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot") else NULL
open_device(out_path(cfg, out, paste0(vis, "_dISF_vs_serum.pdf")), "pdf", fi)
for (p in compact(list(p1, p2, p3))) print(set_family(p, fi$family))
close_device(fi)
msg("Visit analysis %s: %s", vis, file.path(cfg$paths$output, out))
