# 03 - Aim 2: Olink galanin on its own. The galanin ELISA is ignored here (its specificity is not
# established); only the Olink assay GAL is used.
#   2a clinical parameters: Olink GAL vs every analysed clinical variable - Spearman; partial Spearman
#      adjusted for age, sex and Olink plate; women and men separately; leaving out one person at a
#      time; BH-FDR over the parameters. "Robust" (exploratory criteria): p < 0.05 unadjusted and
#      adjusted in the same direction, the same direction in women and in men, and still p < 0.05
#      whichever person is left out. Technical factors: Olink plate, flagged samples.
#   2b other Olink proteins: Olink GAL vs every measurable protein (Spearman; adjusted; BH-FDR);
#      proteins stored and released together with galanin (pre-specified, config: neuroendocrine_proteins,
#      tested even if below LOD; FDR within the set); gene-set enrichment (GSEA) of the correlation
#      ranking; the main axes of the serum proteome (principal components) and how the proteins most
#      correlated with GAL correlate with each other (heatmap)
# Out: output/2_olink_galanin/

source("R/utils.R")
source("R/leip.R")
cfg <- load_config()
set.seed(cfg$seed %||% 1)
d <- leip_load(cfg)
S <- d$S
fdr_cut <- cfg$fdr %||% 0.05
out <- "2_olink_galanin"
clear_outputs(cfg, out)
old <- file.path(cfg$paths$output, "2_galanin_correlates")           # output of the earlier version of this step
if (dir.exists(old)) { unlink(old, recursive = TRUE); msg("Removed %s (output of the earlier version of this step)", old) }
Q <- "2 Olink galanin on its own (ELISA ignored)"
ans <- answers_new()
if (!d$has_gal) {
  answer(ans, Q, "Olink galanin", "not possible", "Galanin (GAL) is not in the Olink data.")
  answers_save(ans, cfg, out, "answers.csv")
  quit(save = "no")
}
g <- d$npx
has_sex <- "sex_male" %in% names(S) && n_distinct(na.omit(S$sex_male)) == 2
sexlab <- if (has_sex) if_else(S$sex_male == 1, "men", "women") else rep("all", nrow(S))
chance_max <- \(k, n) r_crit(log(2) / k, n)          # typical strongest |rho| among k correlations by chance alone (n samples)
kw_p <- \(v, grp) { ok <- !is.na(v) & !is.na(grp); if (n_distinct(grp[ok]) < 2) NA_real_ else kruskal.test(v[ok], factor(grp[ok]))$p.value }

# ---- 2a: clinical parameters ---------------------------------------------------------------------------------------------
params <- setdiff(d$pinfo$parameter[d$pinfo$analysed], "galanin_elisa")
loo <- \(x, y) {                                     # Spearman rho, leaving out each person in turn
  k <- which(!is.na(x) & !is.na(y))
  map_dbl(k, \(i) { j <- setdiff(k, i); suppressWarnings(cor(x[j], y[j], method = "spearman")) })
}
in_sex <- \(x, sx) spearman1(x[S$sex_male %in% sx], g[S$sex_male %in% sx], min_n = 8)
clin <- map(params, \(pm) {
  x <- S[[pm]]
  u <- spearman1(x, g, min_n = 10)
  a <- spearman1(x, g, { z <- S[setdiff(d$cov, pm)]; if (ncol(z)) z }, min_n = 10)
  w <- if (has_sex && pm != "sex_male") in_sex(x, 0); mn <- if (has_sex && pm != "sex_male") in_sex(x, 1)
  lr <- loo(x, g)
  tibble(parameter = pm, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
         rho_adj = a$rho, ci_low_adj = a$ci_low, ci_high_adj = a$ci_high, p_adj = a$p,
         rho_women = w$rho %||% NA_real_, p_women = w$p %||% NA_real_, rho_men = mn$rho %||% NA_real_, p_men = mn$p %||% NA_real_,
         loo_min = min(lr), loo_max = max(lr),
         npx_difference = if (d$pinfo$type[d$pinfo$parameter == pm] == "binary") unname(median_diff(x, cbind(g))) else NA_real_)
}) |> bind_rows() |>
  mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH"),
         weakest_loo = if_else(rho > 0, loo_min, -loo_max),          # the weakest correlation when one person is left out
         robust = coalesce(p < 0.05 & p_adj < 0.05 & sign(rho_adj) == sign(rho) &
                             (is.na(rho_women) | sign(rho_women) == sign(rho)) & (is.na(rho_men) | sign(rho_men) == sign(rho)) &
                             weakest_loo >= r_crit(0.05, n - 1), FALSE)) |>
  left_join(d$pinfo |> select(parameter, label, type, key), by = "parameter") |>
  relocate(parameter, label, type, key) |> relocate(fdr, .after = p) |> relocate(fdr_adj, .after = p_adj) |> arrange(p)

# technical factors: Olink plate and flagged samples
gz <- zscore(g)
flags <- list(`Olink plate` = S$plate, `Olink QC warning` = col_or(S, "SampleQC") == "WARN", `Olink outlier` = col_or(S, "qc_outlier"),
              `low serum volume` = col_or(S, "low_serum") == 1, lipaemic = col_or(S, "lipamisch") == 1)
tech <- imap(flags, \(f, nm) {
  if (nm == "Olink plate") {
    ok <- !is.na(f) & !is.na(g); m <- tapply(g[ok], f[ok], median); m <- m[str_order(names(m), numeric = TRUE)]
    return(tibble(factor = nm, n_flagged = NA_integer_, detail = paste(sprintf("%s: %.2f", names(m), m), collapse = "; "), p = kw_p(g, f)))
  }
  f <- coalesce(as.logical(f), FALSE) & !is.na(g)
  if (!any(f)) return(tibble(factor = nm, n_flagged = 0L, detail = "no flagged sample", p = NA_real_))
  tibble(factor = nm, n_flagged = sum(f), detail = paste(sprintf("%s: GAL z = %+.1f", S$SubjectID[f], gz[f]), collapse = "; "),
         p = if (sum(f) >= 3 && sum(!f) >= 3) wilcox.test(g[f], g[!f])$p.value else NA_real_)
}) |> bind_rows()

# ---- 2b: other Olink proteins --------------------------------------------------------------------------------------------
Mp <- d$M[, setdiff(colnames(d$M), d$gal), drop = FALSE]
u <- spearman_vs(g, Mp); a <- spearman_vs(g, Mp, S[d$cov])
prot <- tibble(OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p, rho_adj = a$rho, p_adj = a$p) |>
  mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH")) |>
  left_join(d$det |> select(OlinkID, Assay, protein, frac_above_lod), by = "OlinkID") |>
  relocate(OlinkID, Assay, protein) |> relocate(fdr, .after = p) |> arrange(p)
n_prot <- sum(!is.na(prot$p))

# pre-specified: proteins stored and released together with galanin
ne <- tibble(gene = unlist(cfg$neuroendocrine_proteins)) |> mutate(OlinkID = map_chr(gene, \(x) find_assay(d$det, x)[1]))
ne_found <- ne |> filter(!is.na(OlinkID), OlinkID != d$gal)
ne_missing <- ne$gene[is.na(ne$OlinkID)]
ne_tab <- if (nrow(ne_found)) {
  Yn <- d$Y[, ne_found$OlinkID, drop = FALSE]
  un <- spearman_vs(g, Yn, min_n = 5); an <- spearman_vs(g, Yn, S[d$cov], min_n = 5)
  tibble(gene = ne_found$gene, OlinkID = un$id, n = un$n, rho = un$rho, ci_low = un$ci_low, ci_high = un$ci_high, p = un$p,
         rho_adj = an$rho, p_adj = an$p) |>
    mutate(fdr_in_set = p.adjust(p, "BH")) |> left_join(d$det |> select(OlinkID, frac_above_lod, measurable), by = "OlinkID")
} else tibble()

# gene sets
gsea <- NULL
if (isTRUE(cfg$gsea$run %||% TRUE) && ncol(Mp) >= 50)
  gsea <- tryCatch(rank_gsea(prot$Assay, prot$rho, cfg$gsea$collections %||% c("H", "C2:CP:REACTOME", "C5:GO:BP"),
                             cfg$gsea$min_size %||% 10, cfg$gsea$max_size %||% 500),
                   error = \(e) { msg("Gene-set enrichment skipped: %s", conditionMessage(e)); NULL })

# main axes of the serum proteome: principal components of all measurable proteins (missing values: protein median).
# For each: does GAL follow it, and what else does it reflect (plate, clinical variables, the proteins that load on it)?
X <- apply(Mp, 2, \(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v })
X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
pcs <- tibble(); pc_scores <- NULL
if (ncol(X) >= 10 && nrow(X) >= 10) {
  pca <- prcomp(X, center = TRUE, scale. = TRUE)
  k <- min(8, ncol(pca$x)); ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pc_scores <- pca$x[, seq_len(k), drop = FALSE]
  pcs <- map(seq_len(k), \(i) {
    s <- pc_scores[, i]; r <- spearman1(s, g)
    cl <- spearman_vs(s, S[params], min_n = 10) |> filter(!is.na(rho)) |> slice_max(abs(rho), n = 1, with_ties = FALSE)
    ld <- pca$rotation[, i] * pca$sdev[i]                               # correlation of each protein with the component
    top <- head(order(-abs(ld)), 6)
    tibble(component = paste0("PC", i), variance_pct = ve[i], rho_GAL = r$rho, ci_low = r$ci_low, ci_high = r$ci_high, p_GAL = r$p,
           plate_p = kw_p(s, S$plate),
           strongest_clinical = if (nrow(cl)) sprintf("%s (%+.2f)", d$pinfo$label[match(cl$id, d$pinfo$parameter)], cl$rho) else "",
           top_proteins = paste(sprintf("%s (%+.2f)", d$det$protein[match(colnames(X)[top], d$det$OlinkID)], ld[top]), collapse = ", "))
  }) |> bind_rows()
}

# ---- answers ---------------------------------------------------------------------------------------------------------------
h <- clin |> filter(p < 0.05); rb <- clin |> filter(robust)
answer(ans, Q, "2a Which clinical parameters go with Olink galanin?",
       sprintf("%d of %d parameters at p < 0.05 (%d at FDR < %g); about %.1f expected by chance", nrow(h), sum(!is.na(clin$p)),
               sum(clin$fdr < fdr_cut, na.rm = TRUE), fdr_cut, 0.05 * sum(!is.na(clin$p))),
       paste0(if (nrow(h)) paste(head(sprintf("%s %+.2f (p = %s; adjusted %+.2f, p = %s)", h$label, h$rho, fmt_p(h$p), h$rho_adj, fmt_p(h$p_adj)), 15),
                                 collapse = "; ") else "None",
              sprintf(". Adjusted = partial Spearman for %s (a covariate is left out when it is the parameter itself). With %d parameters, the strongest correlation expected by chance alone is about |rho| = %.2f.",
                      covs_label(d$cov), nrow(clin), chance_max(nrow(clin), median(clin$n, na.rm = TRUE)))))
answer(ans, Q, "2a Robust associations (exploratory criteria)",
       if (nrow(rb)) sprintf("%d: %s", nrow(rb), paste(rb$label, collapse = ", ")) else "none",
       paste0("Robust = p < 0.05 unadjusted and adjusted, in the same direction; the same direction in women and in men; still p < 0.05 whichever person is left out. ",
              if (nrow(rb)) paste0(paste(sprintf("%s: rho %+.2f, adjusted %+.2f, women %+.2f, men %+.2f, leaving one out %+.2f to %+.2f",
                                                 rb$label, rb$rho, rb$rho_adj, rb$rho_women, rb$rho_men, rb$loo_min, rb$loo_max), collapse = "; "), ".")
              else "No parameter meets all criteria."))
pl <- tech |> filter(factor == "Olink plate"); fl <- tech |> filter(factor != "Olink plate", n_flagged > 0)
answer(ans, Q, "2a Technical factors: does Olink galanin depend on the plate or on flagged samples?",
       if (isTRUE(pl$p < 0.05)) sprintf("plate effect (p = %s): the adjusted results account for it", fmt_p(pl$p)) else "no clear plate effect",
       paste0(if (nrow(pl)) sprintf("Median NPX per plate: %s (p = %s). ", pl$detail, fmt_p(pl$p)) else "",
              if (nrow(fl)) paste0(paste(sprintf("%s: %s", fl$factor, fl$detail), collapse = "; "), ". ") else "No flagged samples. ",
              sprintf("GAL is above LOD in %.0f%% of the samples.", 100 * d$det$frac_above_lod[d$det$OlinkID == d$gal])))
s <- prot |> filter(fdr < fdr_cut); pos <- prot |> filter(rho > 0); neg <- prot |> filter(rho < 0)
answer(ans, Q, "2b Which Olink proteins go with Olink galanin?",
       sprintf("%d of %d proteins at FDR < %g (%d positive, %d negative); adjusted for %s: %d", nrow(s), n_prot, fdr_cut,
               sum(s$rho > 0), sum(s$rho < 0), covs_label(d$cov), sum(prot$fdr_adj < fdr_cut, na.rm = TRUE)),
       sprintf("Strongest positive: %s. Strongest negative: %s. By chance alone about %d proteins reach p < 0.05, and the strongest of %d correlations is expected around |rho| = %.2f.",
               top_str(pos$protein, pos$rho, pos$p, 10), top_str(neg$protein, neg$rho, neg$p, 5), round(0.05 * n_prot), n_prot,
               chance_max(n_prot, median(prot$n, na.rm = TRUE))))
if (nrow(ne_tab)) {
  x <- ne_tab |> arrange(p); hn <- x |> filter(p < 0.05)
  answer(ans, Q, "2b Proteins stored and released together with galanin (neuroendocrine vesicles; pre-specified)",
         paste0(if (nrow(hn)) sprintf("%d of %d at p < 0.05 (%d at FDR < %g within the set)", nrow(hn), nrow(x), sum(x$fdr_in_set < fdr_cut), fdr_cut)
                else sprintf("none of %d at p < 0.05", nrow(x)),
                sprintf("; adjusted for %s: %d at p < 0.05", covs_label(d$cov), sum(x$p_adj < 0.05, na.rm = TRUE))),
         paste0(paste(sprintf("%s %+.2f (p = %s; adjusted %+.2f, p = %s)", x$gene, x$rho, fmt_p(x$p), x$rho_adj, fmt_p(x$p_adj)), collapse = "; "),
                if (length(ne_missing)) sprintf(". Not on the Olink panel: %s.", paste(ne_missing, collapse = ", ")) else "."))
}
if (!is.null(gsea)) {
  gs <- gsea |> filter(padj < fdr_cut)
  answer(ans, Q, "2b Gene sets (GSEA of the ranking by correlation with Olink galanin)",
         if (nrow(gs)) sprintf("%d gene sets at padj < %g", nrow(gs), fdr_cut) else sprintf("no gene set at padj < %g", fdr_cut),
         paste0(if (nrow(gs)) paste(head(sprintf("%s (NES %+.1f)", gs$pathway, gs$NES), 8), collapse = "; ") else "",
                if (nrow(gs)) ". " else "", "Top by p: ",
                paste(head(sprintf("%s (NES %+.1f, padj %s)", gsea$pathway, gsea$NES, fmt_p(gsea$padj)), 5), collapse = "; "), "."))
}
if (nrow(pcs)) {
  b <- pcs |> slice_min(p_GAL, n = 1, with_ties = FALSE)
  answer(ans, Q, "2b Does Olink galanin follow a main axis of variation of the serum proteome?",
         if (isTRUE(b$p_GAL < 0.05 / nrow(pcs))) sprintf("yes: %s (%.0f%% of the variance), rho %+.2f", b$component, b$variance_pct, b$rho_GAL)
         else if (isTRUE(b$p_GAL < 0.05)) sprintf("weakly: %s (%.0f%% of the variance), rho %+.2f (p = %s; %d components tested)", b$component,
                                                  b$variance_pct, b$rho_GAL, fmt_p(b$p_GAL), nrow(pcs))
         else sprintf("no: |rho| with each of the first %d components below %.2f", nrow(pcs), max(abs(pcs$rho_GAL), na.rm = TRUE) + 0.005),
         paste(sprintf("%s (%.0f%% of the variance): GAL rho %+.2f (p = %s); plate p = %s; strongest clinical link %s; proteins %s",
                       pcs$component, pcs$variance_pct, pcs$rho_GAL, fmt_p(pcs$p_GAL), fmt_p(pcs$plate_p), pcs$strongest_clinical,
                       pcs$top_proteins), collapse = ". "))
}
answers <- answers_save(ans, cfg, out, "answers.csv")

# ---- figures ---------------------------------------------------------------------------------------------------------------
F <- figs_new()
# 2a: forest plot - parameters with p < 0.05 (either analysis) and the key parameters
sel <- clin |> filter(coalesce(p < 0.05 | p_adj < 0.05, FALSE) | key)
if (nrow(sel)) {
  lev <- sel |> arrange(rho) |> pull(label)
  x <- bind_rows(sel |> transmute(label, analysis = "Spearman", rho, ci_low, ci_high, p),
                 sel |> transmute(label, analysis = paste("adjusted for", covs_label(d$cov)), rho = rho_adj, ci_low = ci_low_adj, ci_high = ci_high_adj, p = p_adj)) |>
    mutate(analysis = factor(analysis, unique(analysis)), y = match(label, lev) + if_else(analysis == "Spearman", 0.17, -0.17),
           sig = coalesce(p < 0.05, FALSE))
  p1 <- ggplot(x, aes(rho, y, colour = analysis)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), alpha = 0.45) + geom_point(aes(shape = sig), size = 2.1) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    scale_colour_manual(values = setNames(c("black", "#8e44ad"), levels(x$analysis)), name = NULL) +
    scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
    labs(title = "Olink galanin vs clinical parameters (LEIP; the ELISA is not used)",
         subtitle = sprintf("parameters with p < 0.05 in either analysis, and the key parameters; Spearman rho with 95%% CI; %d parameters tested", nrow(clin)),
         x = "rho with Olink GAL", y = NULL) + theme(legend.position = "bottom")
  fig(F, "clinical", p1, cfg, out, "olink_galanin_vs_clinical.png", width = 9, height = 2.5 + 0.25 * length(lev))
}
# 2a: scatter plots of the strongest associations
tq <- clin |> filter(!is.na(p)) |> head(9)
if (nrow(tq)) {
  sc <- map(seq_len(nrow(tq)), \(i) tibble(panel = sprintf("%s: rho %+.2f (adjusted %+.2f)%s", tq$label[i], tq$rho[i], tq$rho_adj[i], if (tq$robust[i]) ", robust" else ""),
                                           x = S[[tq$parameter[i]]], y = g, sex = sexlab)) |>
    bind_rows() |> mutate(panel = factor(panel, unique(panel))) |> filter(!is.na(x), !is.na(y))
  p2 <- ggplot(sc, aes(x, y)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30", linewidth = 0.6) +
    geom_point(aes(colour = sex), size = 1.6) + facet_wrap(~panel, scales = "free_x") +
    scale_colour_manual(values = c(women = "#c0392b", men = "#2471a3", all = "grey30")) +
    labs(title = "Olink galanin vs its most strongly associated clinical parameters", x = NULL, y = "Olink GAL (NPX)", colour = NULL)
  fig(F, "clinical_scatter", p2, cfg, out, "olink_galanin_top_clinical.png", width = 11, height = 8.5)
}
# 2b: volcano and the strongest proteins
r <- prot |> filter(!is.na(p)) |> mutate(neuroendocrine = OlinkID %in% ne_found$OlinkID, sig = fdr < fdr_cut)
fl <- if (any(r$fdr < fdr_cut)) max(r$p[r$fdr < fdr_cut]) else NA
p3 <- ggplot(r, aes(rho, -log10(p))) +
  (if (!is.na(fl)) geom_hline(yintercept = -log10(fl), linetype = 2, colour = "grey50")) +
  geom_point(aes(colour = sig), size = 0.8, alpha = 0.7) +
  geom_point(data = r |> filter(neuroendocrine), shape = 21, size = 2.6, colour = "darkorange3", stroke = 0.9) +
  geom_text(data = bind_rows(head(r, 15), r |> filter(neuroendocrine, p < 0.05)) |> distinct(OlinkID, .keep_all = TRUE),
            aes(label = protein), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"),
                      labels = c(`FALSE` = "not significant", `TRUE` = sprintf("FDR < %g", fdr_cut)), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(title = "Which Olink proteins go with Olink galanin? (LEIP serum)",
       subtitle = sprintf("%d proteins; orange circles: proteins released together with galanin (neuroendocrine vesicles)%s", nrow(r),
                          if (!is.na(fl)) "; dashed: FDR cutoff" else "; no protein passes the FDR"),
       x = "Spearman rho with Olink GAL", y = "-log10 p") + theme(legend.position = "bottom")
fig(F, "volcano", p3, cfg, out, "olink_galanin_vs_all_proteins.png", width = 9, height = 7)
tp <- head(r, 9)
if (nrow(tp)) {
  sc <- map(seq_len(nrow(tp)), \(i) tibble(panel = sprintf("%s: rho %+.2f", tp$protein[i], tp$rho[i]), x = d$Y[, tp$OlinkID[i]], y = g, sex = sexlab)) |>
    bind_rows() |> mutate(panel = factor(panel, unique(panel))) |> filter(!is.na(x), !is.na(y))
  p4 <- ggplot(sc, aes(x, y)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
    geom_point(aes(colour = sex), size = 1.4) + facet_wrap(~panel, scales = "free_x") +
    scale_colour_manual(values = c(women = "#c0392b", men = "#2471a3", all = "grey30")) +
    labs(title = "Olink galanin vs its most strongly correlated proteins", x = "NPX of the other protein", y = "Olink GAL (NPX)", colour = NULL)
  fig(F, "top_proteins", p4, cfg, out, "olink_galanin_top_proteins.png", width = 10, height = 8)
}
if (nrow(ne_tab)) {
  x <- ne_tab |> arrange(rho) |> mutate(y = row_number(), sig = coalesce(p < 0.05, FALSE),
                                        name = if_else(coalesce(measurable, TRUE), gene, paste0(gene, " (below LOD)")))
  p5 <- ggplot(x, aes(rho, y)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = "firebrick", alpha = 0.45) +
    geom_point(aes(shape = sig), colour = "firebrick", size = 2.2) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    scale_y_continuous(breaks = x$y, labels = x$name, expand = expansion(add = 0.6)) +
    labs(title = "Olink galanin vs proteins released together with it (neuroendocrine dense-core vesicles)",
         subtitle = "pre-specified set; Spearman rho with 95% CI", x = "rho with Olink GAL", y = NULL) + theme(legend.position = "bottom")
  fig(F, "neuroendocrine", p5, cfg, out, "neuroendocrine_proteins.png", width = 8, height = 2 + 0.3 * nrow(x))
}
if (!is.null(gsea) && nrow(gsea)) {
  x <- gsea |> slice_min(pval, n = 20, with_ties = FALSE) |> mutate(sig = padj < fdr_cut, name = str_trunc(pathway, 60))
  p6 <- ggplot(x, aes(NES, reorder(name, NES), size = size, colour = sig)) + geom_point() + geom_vline(xintercept = 0, colour = "grey60") +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"), labels = c(`FALSE` = "padj >= 0.05", `TRUE` = "padj < 0.05"), name = NULL) +
    labs(title = "Gene sets among the proteins correlated with Olink galanin (GSEA, top 20 by p)",
         subtitle = "NES > 0: enriched among proteins positively correlated with GAL", x = "normalised enrichment score", y = NULL)
  fig(F, "gene_sets", p6, cfg, out, "gene_sets.png", width = 11, height = 7)
}
if (nrow(pcs)) {
  x <- pcs |> mutate(name = factor(sprintf("%s (%.0f%%): %s", component, variance_pct, strongest_clinical),
                                   rev(sprintf("%s (%.0f%%): %s", component, variance_pct, strongest_clinical))),
                     sig = coalesce(p_GAL < 0.05, FALSE))
  p7 <- ggplot(x, aes(rho_GAL, name)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high, shape = sig), colour = "firebrick") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    labs(title = "Does Olink galanin follow a main axis of the serum proteome?",
         subtitle = "principal components of all measurable proteins: % of variance and the clinical parameter each is most correlated with",
         x = "Spearman rho of the component with Olink GAL", y = NULL) + theme(legend.position = "bottom")
  fig(F, "proteome_axes", p7, cfg, out, "proteome_axes.png", width = 10, height = 5.5)
}
top_ids <- head(prot$OlinkID[!is.na(prot$p)], 20)
if (length(top_ids) >= 3) {
  H <- cbind(g, d$Y[, top_ids, drop = FALSE])
  colnames(H) <- make.unique(c("GAL (Olink galanin)", d$det$protein[match(top_ids, d$det$OlinkID)]))
  R <- suppressWarnings(cor(H, method = "spearman", use = "pairwise.complete.obs"))
  ord <- colnames(R)[hclust(as.dist(1 - R))$order]
  hm <- as_tibble(R, rownames = "a") |> pivot_longer(-a, names_to = "b", values_to = "rho") |>
    mutate(a = factor(a, ord), b = factor(b, ord))
  p8 <- ggplot(hm, aes(a, b, fill = rho)) + geom_tile() + geom_text(aes(label = sprintf("%.1f", rho)), size = 2.1) +
    scale_fill_gradient2(low = "steelblue", high = "firebrick", limits = c(-1, 1)) +
    labs(title = "Olink galanin and its 20 most correlated proteins: do they form groups?",
         subtitle = "Spearman rho between all pairs, clustered; a block of proteins correlated with each other points to one shared process",
         x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  fig(F, "heatmap", p8, cfg, out, "top_proteins_heatmap.png", width = 10, height = 9)
}
figs_save(F, cfg, out, "figures.rds")

# ---- tables ----------------------------------------------------------------------------------------------------------------
save_csv(clin, cfg, out, "olink_galanin_vs_clinical.csv")
save_csv(tech, cfg, out, "technical_factors.csv")
save_csv(prot, cfg, out, "olink_galanin_vs_proteins.csv")
if (nrow(ne_tab)) save_csv(ne_tab, cfg, out, "neuroendocrine_proteins.csv")
if (!is.null(gsea)) save_csv(gsea, cfg, out, "gene_sets.csv")
if (nrow(pcs)) save_csv(pcs, cfg, out, "proteome_axes.csv")
writexl::write_xlsx(list(answers = answers, clinical = clin, technical_factors = tech, proteins = prot, neuroendocrine = ne_tab,
                         gene_sets = gsea %||% tibble(), proteome_axes = pcs) |> keep(\(x) is.data.frame(x) && ncol(x) > 0),
                    out_path(cfg, out, "olink_galanin.xlsx"))
for (i in seq_len(nrow(answers))) msg("%s: %s", answers$item[i], answers$verdict[i])
