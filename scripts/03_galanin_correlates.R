# 03 - Aim 2: what does galanin correlate with?
#   - proteins: Olink GAL vs every measurable Olink protein (Spearman; partial Spearman adjusted for
#     age, sex and Olink plate; BH-FDR over the proteins), the same for the ELISA, and the proteins
#     found with both measurements
#   - proteins stored and released together with galanin (neuroendocrine dense-core vesicles; config:
#     neuroendocrine_proteins): pre-specified, so tested even if below LOD; FDR within this set
#   - gene-set enrichment (GSEA) of the correlation ranking
#   - clinical parameters: Olink GAL and the ELISA vs every analysed clinical variable
# Out: output/2_galanin_correlates/

source("R/utils.R")
source("R/leip.R")
cfg <- load_config()
set.seed(cfg$seed %||% 1)
d <- leip_load(cfg)
S <- d$S
fdr_cut <- cfg$fdr %||% 0.05
clear_outputs(cfg, "2_galanin_correlates")
Q <- "2 What does galanin correlate with?"
ans <- answers_new()
measures <- compact(list(`Olink GAL` = if (d$has_gal) d$npx, `galanin ELISA` = if (d$has_elisa) d$elisa))
if (!length(measures)) {
  answer(ans, Q, "galanin", "not possible", "Neither Olink GAL nor the ELISA is available.")
  answers_save(ans, cfg, "2_galanin_correlates", "answers.csv")
  quit(save = "no")
}
cov_of <- \(m) if (m == "galanin ELISA") setdiff(d$cov, "plate") else d$cov   # the ELISA was not run on the Olink plates

# ---- proteins ------------------------------------------------------------------------------------------------------------
Mp <- d$M[, setdiff(colnames(d$M), d$gal), drop = FALSE]
prot <- imap(measures, \(v, m) {
  u <- spearman_vs(v, Mp)
  a <- spearman_vs(v, Mp, S[d$cov])                                        # proteins are Olink values: plate stays in
  tibble(measure = m, OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
         n_adj = a$n, rho_adj = a$rho, p_adj = a$p) |>
    mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH"))
}) |> bind_rows() |>
  left_join(d$det |> select(OlinkID, Assay, protein, frac_above_lod), by = "OlinkID") |>
  relocate(measure, OlinkID, Assay, protein) |> arrange(measure, p)
wide <- prot |> select(measure, OlinkID, protein, rho, p, fdr) |>
  pivot_wider(names_from = measure, values_from = c(rho, p, fdr), names_glue = "{.value}_{measure}") |>
  rename_with(\(x) str_replace_all(x, " ", "_"))
if (all(c("rho_Olink_GAL", "rho_galanin_ELISA") %in% names(wide)))
  wide <- wide |> mutate(both_p05_same_direction = coalesce(p_Olink_GAL < 0.05 & p_galanin_ELISA < 0.05 &
                                                            sign(rho_Olink_GAL) == sign(rho_galanin_ELISA), FALSE)) |>
    arrange(desc(both_p05_same_direction), p_Olink_GAL)

# ---- pre-specified: proteins released together with galanin ----------------------------------------------------------------------
ne <- tibble(gene = unlist(cfg$neuroendocrine_proteins)) |> mutate(OlinkID = map_chr(gene, \(g) find_assay(d$det, g)[1]))
ne_found <- ne |> filter(!is.na(OlinkID), OlinkID != d$gal)
ne_tab <- if (nrow(ne_found)) imap(measures, \(v, m) {
  Yn <- d$Y[, ne_found$OlinkID, drop = FALSE]
  u <- spearman_vs(v, Yn, min_n = 5); a <- spearman_vs(v, Yn, S[cov_of(m)], min_n = 5)
  tibble(measure = m, gene = ne_found$gene, OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
         rho_adj = a$rho, p_adj = a$p) |> mutate(fdr_in_set = p.adjust(p, "BH"))
}) |> bind_rows() |> left_join(d$det |> select(OlinkID, frac_above_lod, measurable), by = "OlinkID") else tibble()
ne_missing <- ne$gene[is.na(ne$OlinkID)]

# ---- gene sets ------------------------------------------------------------------------------------------------------------------------
gsea <- NULL
if (isTRUE(cfg$gsea$run %||% TRUE) && ncol(Mp) >= 50) {
  gsea <- tryCatch(imap(measures, \(v, m) {
    r <- prot |> filter(measure == m)
    rank_gsea(r$Assay, r$rho, cfg$gsea$collections %||% c("H", "C2:CP:REACTOME", "C5:GO:BP"), cfg$gsea$min_size %||% 10,
              cfg$gsea$max_size %||% 500) |> mutate(measure = m, .before = 1)
  }) |> bind_rows(), error = \(e) { msg("Gene-set enrichment skipped: %s", conditionMessage(e)); NULL })
}

# ---- clinical parameters -------------------------------------------------------------------------------------------------------------
params <- d$pinfo$parameter[d$pinfo$analysed]
clin <- imap(measures, \(v, m) {
  map(setdiff(params, if (m == "galanin ELISA") "galanin_elisa"), \(pm) {
    u <- spearman1(S[[pm]], v, min_n = 10)
    a <- spearman1(S[[pm]], v, { z <- S[setdiff(cov_of(m), pm)]; if (ncol(z)) z }, min_n = 10)
    tibble(measure = m, parameter = pm, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
           rho_adj = a$rho, p_adj = a$p)
  }) |> bind_rows() |> mutate(fdr = p.adjust(p, "BH"))
}) |> bind_rows() |>
  left_join(d$pinfo |> select(parameter, label, type, key), by = "parameter") |>
  relocate(measure, parameter, label, type, key)

# ---- answers ---------------------------------------------------------------------------------------------------------------------------
for (m in names(measures)) {
  r <- prot |> filter(measure == m); s <- r |> filter(fdr < fdr_cut); pos <- r |> filter(rho > 0); neg <- r |> filter(rho < 0)
  answer(ans, Q, sprintf("Proteins correlated with the %s", m),
         sprintf("%d of %d proteins at FDR < %g (%d positive, %d negative); adjusted for %s: %d", nrow(s), sum(!is.na(r$p)), fdr_cut,
                 sum(s$rho > 0), sum(s$rho < 0), covs_label(d$cov), sum(r$fdr_adj < fdr_cut, na.rm = TRUE)),
         sprintf("Strongest positive: %s. Strongest negative: %s. By chance alone about %d proteins reach p < 0.05.",
                 top_str(pos$protein, pos$rho, pos$p, 10), top_str(neg$protein, neg$rho, neg$p, 5), round(0.05 * sum(!is.na(r$p)))))
}
if ("both_p05_same_direction" %in% names(wide)) {
  b <- wide |> filter(both_p05_same_direction)
  answer(ans, Q, "Proteins found with both measurements (p < 0.05 for Olink GAL and the ELISA, same direction)",
         sprintf("%d proteins", nrow(b)),
         if (nrow(b)) paste(head(sprintf("%s (Olink %+.2f, ELISA %+.2f)", b$protein, b$rho_Olink_GAL, b$rho_galanin_ELISA), 15), collapse = ", ")
         else "none")
}
if (nrow(ne_tab)) for (m in unique(ne_tab$measure)) {
  x <- ne_tab |> filter(measure == m) |> arrange(p); h <- x |> filter(p < 0.05)
  answer(ans, Q, sprintf("Proteins released together with galanin (neuroendocrine vesicles) - %s", m),
         if (nrow(h)) sprintf("%d of %d at p < 0.05 (%d at FDR < %g within the set)", nrow(h), nrow(x), sum(x$fdr_in_set < fdr_cut), fdr_cut)
         else sprintf("none of %d at p < 0.05", nrow(x)),
         paste0(paste(sprintf("%s %+.2f (p = %s)", x$gene, x$rho, fmt_p(x$p)), collapse = "; "),
                if (length(ne_missing)) sprintf(". Not on the Olink panel: %s.", paste(ne_missing, collapse = ", ")) else "."))
}
if (!is.null(gsea)) for (m in unique(gsea$measure)) {
  g <- gsea |> filter(measure == m, padj < fdr_cut)
  answer(ans, Q, sprintf("Gene sets (GSEA) - %s", m), if (nrow(g)) sprintf("%d gene sets at padj < %g", nrow(g), fdr_cut) else "no gene set at padj < 0.05",
         if (nrow(g)) paste(head(sprintf("%s (NES %+.1f)", g$pathway, g$NES), 8), collapse = "; ") else "")
}
for (m in unique(clin$measure)) {
  x <- clin |> filter(measure == m) |> arrange(p); h <- x |> filter(p < 0.05)
  answer(ans, Q, sprintf("Clinical parameters - %s", m),
         if (nrow(h)) sprintf("%d of %d at p < 0.05 (%d at FDR < %g)", nrow(h), nrow(x), sum(x$fdr < fdr_cut), fdr_cut) else sprintf("none of %d at p < 0.05", nrow(x)),
         paste0(if (nrow(h)) paste(head(sprintf("%s %+.2f (p = %s; adjusted %+.2f, p = %s)", h$label, h$rho, fmt_p(h$p),
                                                   h$rho_adj, fmt_p(h$p_adj)), 12), collapse = "; ")
                else sprintf("Strongest: %s", paste(head(sprintf("%s %+.2f (p = %s)", x$label, x$rho, fmt_p(x$p)), 3), collapse = "; ")),
                sprintf(". Adjusted = partial Spearman for %s (a covariate is left out when it is the parameter itself). By chance alone about %.1f of %d reach p < 0.05.",
                        covs_label(cov_of(m)), 0.05 * nrow(x), nrow(x))))
}
answers <- answers_save(ans, cfg, "2_galanin_correlates", "answers.csv")

# ---- figures ------------------------------------------------------------------------------------------------------------------------------
F <- figs_new()
for (m in names(measures)) {
  r <- prot |> filter(measure == m, !is.na(p)) |> mutate(neuroendocrine = OlinkID %in% ne_found$OlinkID, sig = fdr < fdr_cut)
  fl <- if (any(r$fdr < fdr_cut)) max(r$p[r$fdr < fdr_cut]) else NA
  tag <- if (m == "Olink GAL") "GAL" else "ELISA"
  p <- ggplot(r, aes(rho, -log10(p))) +
    (if (!is.na(fl)) geom_hline(yintercept = -log10(fl), linetype = 2, colour = "grey50")) +
    geom_point(aes(colour = sig), size = 0.8, alpha = 0.7) +
    geom_point(data = r |> filter(neuroendocrine), shape = 21, size = 2.6, colour = "darkorange3", stroke = 0.9) +
    geom_text(data = bind_rows(head(r, 15), r |> filter(neuroendocrine, p < 0.05)) |> distinct(OlinkID, .keep_all = TRUE),
              aes(label = protein), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"),
                        labels = c(`FALSE` = "not significant", `TRUE` = sprintf("FDR < %g", fdr_cut)), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    labs(title = sprintf("Which proteins correlate with the %s? (LEIP serum)", m),
         subtitle = sprintf("%d proteins; dashed: FDR cutoff; orange circles: proteins released together with galanin (neuroendocrine vesicles)", nrow(r)),
         x = "Spearman rho", y = "-log10 p") + theme(legend.position = "bottom")
  fig(F, paste0("volcano_", tag), p, cfg, "2_galanin_correlates", sprintf("%s_vs_all_proteins.png", tag), width = 9, height = 7)
  tq <- head(r, 9)
  sc <- map(seq_len(nrow(tq)), \(i) tibble(panel = sprintf("%s: rho = %.2f", tq$protein[i], tq$rho[i]), x = d$Y[, tq$OlinkID[i]], y = measures[[m]])) |>
    bind_rows() |> mutate(panel = factor(panel, levels = unique(panel))) |> filter(!is.na(x), !is.na(y))
  p <- ggplot(sc, aes(x, y)) + geom_point(size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
    facet_wrap(~panel, scales = "free_x") + (if (m == "galanin ELISA") scale_y_continuous(trans = "log2")) +
    labs(title = sprintf("The %s vs its most strongly correlated proteins", m), x = "NPX of the other protein",
         y = if (m == "galanin ELISA") "galanin ELISA (pg/mL, log2 scale)" else "Olink GAL (NPX)")
  fig(F, paste0("top_", tag), p, cfg, "2_galanin_correlates", sprintf("%s_top_proteins.png", tag), width = 10, height = 8)
}
if (nrow(ne_tab)) {
  lev <- ne_tab |> group_by(gene) |> summarise(r = max(rho, na.rm = TRUE)) |> arrange(r) |> pull(gene)
  x <- ne_tab |> mutate(y = match(gene, lev) + if_else(measure == "Olink GAL", 0.15, -0.15), sig = coalesce(p < 0.05, FALSE),
                        gene = if_else(coalesce(measurable, TRUE), gene, paste0(gene, " (below LOD)")))
  p <- ggplot(x, aes(rho, y, colour = measure)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), alpha = 0.45) + geom_point(aes(shape = sig), size = 2.2) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    scale_colour_manual(values = c(`Olink GAL` = "firebrick", `galanin ELISA` = "steelblue"), name = NULL) +
    scale_y_continuous(breaks = seq_along(lev), labels = x$gene[match(lev, ne_tab$gene)], expand = expansion(add = 0.6)) +
    labs(title = "Galanin vs proteins released together with it (neuroendocrine dense-core vesicles)",
         subtitle = "pre-specified set; Spearman rho with 95% CI", x = "Spearman rho", y = NULL) + theme(legend.position = "bottom")
  fig(F, "neuroendocrine", p, cfg, "2_galanin_correlates", "neuroendocrine_proteins.png", width = 8, height = 2 + 0.3 * length(lev))
}
if (!is.null(gsea) && nrow(gsea)) {
  g <- gsea |> group_by(measure) |> slice_min(pval, n = 15, with_ties = FALSE) |> ungroup() |> mutate(sig = padj < fdr_cut)
  p <- ggplot(g, aes(NES, reorder(str_trunc(pathway, 55), NES), size = size, colour = sig)) + geom_point() +
    geom_vline(xintercept = 0, colour = "grey60") + facet_wrap(~measure, scales = "free_y") +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"), labels = c(`FALSE` = "padj >= 0.05", `TRUE` = "padj < 0.05"), name = NULL) +
    labs(title = "Gene sets among the proteins correlated with galanin (GSEA, top 15 by p)",
         subtitle = "NES > 0: enriched among proteins positively correlated with galanin", x = "normalised enrichment score", y = NULL)
  fig(F, "gene_sets", p, cfg, "2_galanin_correlates", "gene_sets.png", width = 14, height = 7)
}
clin_plot <- function(x, title_extra = "") {
  first <- if ("Olink GAL" %in% x$measure) "Olink GAL" else "galanin ELISA"
  lev <- c(x |> filter(measure == first) |> arrange(rho) |> pull(label), setdiff(unique(x$label), x$label[x$measure == first]))
  x <- x |> mutate(y = match(label, lev) + if_else(measure == "Olink GAL", 0.18, -0.18), sig = coalesce(p < 0.05, FALSE))
  ggplot(x, aes(rho, y, colour = measure)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), alpha = 0.45) + geom_point(aes(shape = sig), size = 2) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    scale_colour_manual(values = c(`Olink GAL` = "firebrick", `galanin ELISA` = "steelblue"), name = NULL) +
    scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
    labs(title = paste0("Galanin vs clinical parameters (LEIP)", title_extra), subtitle = "Spearman rho with 95% CI; filled = p < 0.05 (single tests)",
         x = "Spearman rho", y = NULL) + theme(legend.position = "bottom")
}
if (nrow(clin)) {
  fig(F, "clinical_all", clin_plot(clin), cfg, "2_galanin_correlates", "galanin_vs_clinical_all.png", width = 9, height = 2.5 + 0.2 * n_distinct(clin$label))
  sel <- clin |> group_by(parameter) |> filter(any(coalesce(p < 0.05, FALSE)) | any(key)) |> ungroup()
  fig(F, "clinical", clin_plot(sel, ": key parameters and p < 0.05"), cfg, "2_galanin_correlates", "galanin_vs_clinical.png",
            width = 9, height = 2.5 + 0.22 * n_distinct(sel$label))
}

# ---- tables -------------------------------------------------------------------------------------------------------------------------------
save_csv(prot, cfg, "2_galanin_correlates", "galanin_vs_all_proteins.csv")
save_csv(clin, cfg, "2_galanin_correlates", "galanin_vs_clinical.csv")
if (nrow(ne_tab)) save_csv(ne_tab, cfg, "2_galanin_correlates", "neuroendocrine_proteins.csv")
if (!is.null(gsea)) save_csv(gsea, cfg, "2_galanin_correlates", "gene_sets.csv")
figs_save(F, cfg, "2_galanin_correlates", "figures.rds")
writexl::write_xlsx(list(answers = answers, proteins = prot, proteins_both_measures = wide, neuroendocrine = ne_tab,
                         gene_sets = gsea %||% tibble(), clinical = clin) |> keep(\(x) is.data.frame(x) && ncol(x) > 0),
                    out_path(cfg, "2_galanin_correlates", "galanin_correlates.xlsx"))
for (i in seq_len(nrow(answers))) msg("%s: %s", answers$item[i], answers$verdict[i])
