# 05 - Supplementary: every measurable Olink protein vs every clinical variable of the LEIP file.
# Context for the galanin results (what else relates to HDL, sex, BMI ... in these sera) and a check
# of the data: associations known from population proteomics must show up (config: expected_associations).
#   Spearman (primary) and partial Spearman adjusted for age, sex and Olink plate; BH-FDR over the
#   proteins of each variable.
# Out: output/4_clinical_screen/

source("R/utils.R")
source("R/leip.R")
cfg <- load_config()
d <- leip_load(cfg)
S <- d$S; M <- d$M
fdr_cut <- cfg$fdr %||% 0.05
clear_outputs(cfg, "4_clinical_screen")
pinfo <- d$pinfo
params <- pinfo$parameter[pinfo$analysed]
ptype <- setNames(pinfo$type, pinfo$parameter)
ans <- answers_new()
Q <- "S Supplementary: proteins vs clinical parameters"

assoc <- map(params, \(pm) {
  x <- S[[pm]]
  u <- spearman_vs(x, M, min_n = 10)
  a <- spearman_vs(x, M, { z <- S[setdiff(d$cov, pm)]; if (ncol(z)) z }, min_n = 10)
  tibble(parameter = pm, OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
         median_diff = if (ptype[[pm]] == "binary") unname(median_diff(x, M)) else NA_real_,
         n_adj = a$n, rho_adj = a$rho, p_adj = a$p)
}) |> bind_rows() |>
  group_by(parameter) |> mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH")) |> ungroup() |>
  mutate(significant = coalesce(fdr < fdr_cut, FALSE), significant_adj = coalesce(fdr_adj < fdr_cut, FALSE)) |>
  left_join(d$det |> select(OlinkID, Assay, protein), by = "OlinkID") |>
  left_join(pinfo |> select(parameter, label, type, key), by = "parameter") |>
  relocate(parameter, label, type, key, OlinkID, Assay, protein) |>
  arrange(match(parameter, params), p)
nsig <- assoc |>
  group_by(parameter, label, type, key) |>
  summarise(n_samples = max(n, na.rm = TRUE), proteins_tested = sum(!is.na(p)),
            # positive/negative first: summarise() sees a summary as soon as it is created
            positive = sum(significant & rho > 0), negative = sum(significant & rho < 0), significant = sum(significant),
            positive_adj = sum(significant_adj & rho_adj > 0), negative_adj = sum(significant_adj & rho_adj < 0),
            significant_adjusted = sum(significant_adj),
            p_below_05 = sum(p < 0.05, na.rm = TRUE), expected_by_chance = round(0.05 * sum(!is.na(p))),
            top_proteins = top_str(protein, rho, p), .groups = "drop") |>
  relocate(significant, .before = positive) |> relocate(significant_adjusted, .before = positive_adj) |>
  arrange(desc(significant), desc(significant_adjusted))

sanity <- map(cfg$expected_associations %||% list(), \(e) {
  row <- tibble(protein = e$protein, parameter = e$parameter, label = param_label(e$parameter), expected = e$direction)
  oid <- find_assay(d$det, e$protein)[1]
  if (is.na(oid)) return(row |> mutate(status = "protein not measured"))
  if (!e$parameter %in% names(S) || sum(!is.na(S[[e$parameter]])) < 10) return(row |> mutate(status = "parameter not available"))
  r <- spearman1(S[[e$parameter]], d$Y[, oid], min_n = 10)
  ok <- isTRUE(r$p < 0.05 && sign(r$rho) == if (e$direction == "negative") -1 else 1)
  row |> mutate(frac_above_lod = d$det$frac_above_lod[d$det$OlinkID == oid], n = r$n, rho = r$rho, p = r$p,
                status = if (ok) "recovered" else "not recovered")
}) |> bind_rows()
if (nrow(sanity)) sanity <- ensure(sanity, c("frac_above_lod", "n", "rho", "p")) |> relocate(status, .after = last_col())

hit <- nsig |> filter(significant > 0)
answer(ans, Q, "Which clinical parameters show up in the LEIP serum proteome?",
       sprintf("%d of %d parameters with proteins at FDR < %g (adjusted for %s: %d)", nrow(hit), nrow(nsig), fdr_cut, covs_label(d$cov),
               sum(nsig$significant_adjusted > 0)),
       if (nrow(hit)) paste(head(sprintf("%s: %d at FDR < %g, strongest %s", hit$label, hit$significant, fdr_cut, hit$top_proteins), 10), collapse = "; ")
       else "No protein passes the FDR for any parameter.")
if (nrow(sanity)) {
  tested <- sanity |> filter(status %in% c("recovered", "not recovered"))
  answer(ans, Q, "Sanity check: do associations known from population studies show up?",
         sprintf("%d of %d recovered (p < 0.05, expected direction)%s", sum(tested$status == "recovered"), nrow(tested),
                 if (nrow(tested) < nrow(sanity)) sprintf("; %d not testable", nrow(sanity) - nrow(tested)) else ""),
         paste(sprintf("%s ~ %s: %s", sanity$protein, sanity$label,
                       if_else(is.na(sanity$rho), sanity$status, sprintf("rho %+.2f, p = %s", sanity$rho, fmt_p(sanity$p)))), collapse = "; "))
}
answers <- answers_save(ans, cfg, "4_clinical_screen", "answers.csv")

# ---- figures ----------------------------------------------------------------------------------------------------------------------------
F <- figs_new()
nb_par <- nsig |> filter(significant + significant_adjusted > 0 | key)
if (nrow(nb_par)) {
  nb <- nb_par |>
    select(label, `unadjusted|positive` = positive, `unadjusted|negative` = negative, `adjusted|positive` = positive_adj, `adjusted|negative` = negative_adj) |>
    pivot_longer(-label, names_to = c("analysis", "direction"), names_sep = "\\|", values_to = "n") |>
    mutate(n = if_else(direction == "negative", -n, n),
           label = factor(label, levels = nb_par |> arrange(significant, significant_adjusted) |> pull(label)),
           analysis = factor(analysis, c("unadjusted", "adjusted"), c("Spearman", paste("adjusted for", covs_label(d$cov)))))
  p <- ggplot(nb, aes(n, label, fill = direction)) + geom_col() + geom_vline(xintercept = 0, colour = "grey40") +
    facet_wrap(~analysis) + scale_fill_manual(values = c(positive = "firebrick", negative = "steelblue")) +
    labs(title = sprintf("LEIP: proteins correlated with each clinical parameter (FDR < %g)", fdr_cut),
         subtitle = sprintf("%d proteins per parameter; key parameters and all parameters with hits", ncol(M)),
         x = "proteins (negative < 0 < positive)", y = NULL, fill = NULL)
  fig(F, "nsig", p, cfg, "4_clinical_screen", "n_significant_per_parameter.png", width = 11, height = 2 + 0.2 * nrow(nb_par))
}
hk <- assoc |> filter(key, !is.na(p))
if (nrow(hk)) {
  top_hm <- hk |> group_by(OlinkID) |> summarise(best = min(p)) |> slice_min(best, n = 50, with_ties = FALSE) |> pull(OlinkID)
  x <- hk |> filter(OlinkID %in% top_hm)
  m <- x |> select(protein, label, rho) |> pivot_wider(names_from = label, values_from = rho)
  mm <- as.matrix(m[, -1]); mm[is.na(mm)] <- 0
  ord <- if (nrow(mm) > 2) m$protein[hclust(dist(mm))$order] else m$protein
  p <- ggplot(x |> mutate(protein = factor(protein, levels = ord), label = factor(label, levels = unique(pinfo$label[pinfo$key]))),
              aes(label, protein, fill = rho)) + geom_tile() +
    geom_text(aes(label = if_else(significant, "*", "")), size = 4, vjust = 0.75) +
    scale_fill_gradient2(low = "steelblue", high = "firebrick", limits = c(-1, 1)) +
    labs(title = "LEIP: correlation of proteins with the key clinical parameters",
         subtitle = sprintf("the %d proteins with the strongest association (Spearman rho; * FDR < %g)", length(top_hm), fdr_cut), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  fig(F, "heatmap", p, cfg, "4_clinical_screen", "heatmap_key_parameters.png", width = 10, height = 3 + 0.2 * length(top_hm))
  lab_d <- hk |> group_by(label) |> slice_min(p, n = 3, with_ties = FALSE) |> ungroup()
  p <- ggplot(hk, aes(rho, -log10(p))) + geom_point(aes(colour = significant), size = 0.6, alpha = 0.6) +
    geom_text(data = lab_d, aes(label = protein), size = 2.6, vjust = -0.6, check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"), labels = c(`FALSE` = "not significant", `TRUE` = sprintf("FDR < %g", fdr_cut)), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    facet_wrap(~label) + labs(title = "LEIP: every protein vs each key clinical parameter", x = "Spearman rho", y = "-log10 p") +
    theme(legend.position = "bottom")
  fig(F, "volcano", p, cfg, "4_clinical_screen", "volcano_key_parameters.png", width = 13, height = 10)
}

save_csv(assoc, cfg, "4_clinical_screen", "associations_all.csv.gz")
save_csv(nsig, cfg, "4_clinical_screen", "n_significant_per_parameter.csv")
if (nrow(sanity)) save_csv(sanity, cfg, "4_clinical_screen", "sanity_checks.csv")
figs_save(F, cfg, "4_clinical_screen", "figures.rds")
writexl::write_xlsx(list(answers = answers, parameters = pinfo, n_significant = nsig, significant = assoc |> filter(significant | significant_adj),
                         top10_per_parameter = assoc |> group_by(parameter) |> slice_min(p, n = 10, with_ties = FALSE) |> ungroup(),
                         sanity_checks = sanity) |> keep(\(x) is.data.frame(x) && ncol(x) > 0),
                    out_path(cfg, "4_clinical_screen", "clinical_screen.xlsx"))
for (i in seq_len(nrow(answers))) msg("%s: %s", answers$item[i], answers$verdict[i])
