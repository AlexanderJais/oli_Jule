# 16 - Which proteins correlate with TNFRSF9 (4-1BB / CD137) in dISF?
# Settings: config.yml -> tnfrsf9_correlation (anchor protein, target proteins, main targets IL33 / IL4).
#   0  Detectability: % of dISF samples above LOD for TNFRSF9 and the targets, per site and group.
#      Values below LOD are used as measured (Olink recommendation, no substitution); every targeted
#      correlation is repeated on the sample pairs where both proteins are above LOD.
#   1  Targeted correlations (Spearman rho, p, n): TNFRSF9 vs each target in the lesion site and in
#      non-lesional / healthy skin, (a) pooled, (b) within each group (AD, healthy; relapse, non-relapse;
#      lesional, ex-lesional), (c) per visit. A correlation seen only when pooling, not within groups
#      (healthy, relapse, non-relapse, lesional, ex-lesional), reflects group or skin-state differences rather than
#      a protein-protein relationship (flagged).
#   2  Proteome-wide: TNFRSF9 vs every dISF protein, adjusted for skin state / group, visit and plate:
#      partial Spearman correlation (primary ranking) and a mixed model (dream: protein ~ TNFRSF9 + state
#      + visit + plate + (1|subject)). Ranked lists with FDR; IL33 / IL4 and the other targets flagged.
#   3  Longitudinal: within-patient relationship (repeated-measures correlation on subject-centred values)
#      and visit-to-visit changes (delta-delta Spearman), AD patients, per site.
#   4  Cross-compartment (MicroAD only): dISF TNFRSF9 vs serum TNFRSF9 at the same visit; TNFRSF9 vs the
#      targets in MicroAD serum.
# CPUO samples are not used. Out: output/tnfrsf9_correlation/  (TNFRSF9_correlation.xlsx + figures)

source("R/utils.R")
source("R/models.R")
source("R/correlate.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clear_outputs(cfg, "tnfrsf9_correlation")
tc      <- cfg$tnfrsf9_correlation %||% list()
anchor  <- tc$anchor %||% "TNFRSF9"
targets <- unlist(tc$targets %||% c("IL33", "IL4", "CSF2", "IL6", "IL18", "CXCL8", "IL1RL1", "KIT", "KITLG", "TPSAB1", "TPSB2", "FCER1A"))
mains   <- unlist(tc$main_targets %||% c("IL33", "IL4"))
min_n   <- tc$min_n %||% 5
fdr     <- cfg$stats$fdr
out     <- \(...) file.path("tnfrsf9_correlation", ...)

# ---- proteins ----------------------------------------------------------------------------------------------
a_oid <- find_assay(clean, anchor)
if (is.na(a_oid)) stop(anchor, " is not in the NPX data.")
t_oid <- setNames(vapply(targets, \(g) find_assay(clean, g), ""), targets)
not_on_panel <- names(t_oid)[is.na(t_oid)]
shared_assay <- names(t_oid)[!is.na(t_oid) & duplicated(t_oid) | !is.na(t_oid) & duplicated(t_oid, fromLast = TRUE)]
t_oid <- t_oid[!is.na(t_oid)]
msg("Anchor %s (%s); targets on the panel: %s; not on the panel: %s", anchor, a_oid, paste(names(t_oid), collapse = ", "),
    if (length(not_on_panel)) paste(not_on_panel, collapse = ", ") else "-")
assay_name <- clean |> distinct(OlinkID, Assay) |> (\(d) setNames(d$Assay, d$OlinkID))()

# ---- dISF samples: site and group --------------------------------------------------------------------------
isf <- meta |>
  filter(matrix == "ISF", group %in% c("AD", "HC")) |>
  mutate(site_lab = case_when(group == "AD" & site == "L" ~ "lesion site",
                              group == "AD" & site == "NL" ~ "non-lesional skin",
                              group == "HC" ~ "healthy skin"),
         stratum = if_else(site_lab == "lesion site", "lesion site (AD)", "non-lesional / healthy skin"),
         grp = recode(group, HC = "healthy"),
         rel = if_else(relapse %in% c("relapse", "non-relapse"), relapse, NA_character_),
         visit_f = factor(visit, levels = paste0("V", 1:12)) |> droplevels(),
         state_g = if_else(group == "HC", "healthy", state)) |>
  filter(!is.na(site_lab), SampleID %in% unique(clean$SampleID[clean$matrix == "ISF"]))
vals <- clean |> filter(matrix == "ISF", OlinkID %in% c(a_oid, t_oid)) |> select(SampleID, OlinkID, value, below_lod)
val_of <- \(oid, col = "value") { v <- vals |> filter(OlinkID == oid); setNames(v[[col]], v$SampleID) }
isf$anchor <- val_of(a_oid)[isf$SampleID]
isf$anchor_blod <- val_of(a_oid, "below_lod")[isf$SampleID]

# ---- 0. detectability ----------------------------------------------------------------------------------------
lod_tab <- clean |> filter(matrix == "ISF", OlinkID %in% c(a_oid, t_oid)) |>
  inner_join(isf |> select(SampleID, site_lab, grp, state_g), by = "SampleID") |>
  mutate(protein = assay_name[OlinkID])
det_by <- \(...) lod_tab |> group_by(protein, ...) |>
  summarise(n = n(), pct_above_LOD = round(100 * mean(!below_lod, na.rm = TRUE), 1), median_NPX = median(value, na.rm = TRUE),
            median_LOD = median(LOD, na.rm = TRUE), .groups = "drop")
lod_summary <- bind_rows(
  det_by() |> mutate(site = "all dISF", group = "all"),
  det_by(site_lab, grp) |> rename(site = site_lab, group = grp),
  det_by(site_lab, state_g) |> filter(site_lab == "lesion site") |> transmute(protein, site = site_lab, group = paste("AD", state_g), n, pct_above_LOD, median_NPX, median_LOD)) |>
  relocate(protein, site, group) |> arrange(match(protein, assay_name[c(a_oid, t_oid)]), site, group)
overall_det <- lod_summary |> filter(site == "all dISF") |> select(protein, pct_above_LOD)
det_notes <- c(
  "Below-LOD handling: values below LOD are used as measured (no substitution, Olink recommendation). Each targeted correlation is also given for the sample pairs where both proteins are above LOD (rho_both_above_LOD).",
  map_chr(seq_len(nrow(overall_det)), \(i) {
    d <- overall_det[i, ]
    if (is.na(d$pct_above_LOD)) return(sprintf("%s: no LOD available for dISF - detectability unknown.", d$protein))
    sprintf("%s: %.0f%% of dISF samples above LOD%s", d$protein, d$pct_above_LOD,
            if (d$pct_above_LOD < 50) " - LARGELY BELOW LOD: correlations with this protein mostly reflect background noise and must be interpreted with great caution." else
              if (d$pct_above_LOD < 80) " - partly below LOD: check rho_both_above_LOD." else ".")
  }),
  if (length(not_on_panel)) sprintf("Not measured on the Olink Explore HT panel: %s.", paste(not_on_panel, collapse = ", ")),
  if (length(shared_assay)) sprintf("Measured by one combined assay: %s.", paste(shared_assay, collapse = ", ")))
walk(det_notes, \(x) msg("%s", x))

# ---- 1. targeted correlations -----------------------------------------------------------------------------------
subsets <- function(d) {
  s <- list(`pooled (all groups)` = rep(TRUE, nrow(d)),
            `AD` = d$grp == "AD", `healthy` = d$grp == "healthy",
            `AD relapse` = d$rel %in% "relapse", `AD non-relapse` = d$rel %in% "non-relapse",
            `AD lesional` = d$state_g %in% "lesional", `AD ex-lesional` = d$state_g %in% "ex-lesional")
  for (v in levels(d$visit_f)) s[[paste("visit", v)]] <- d$visit_f %in% v
  s
}
targeted <- map(names(t_oid), \(tg) {
  d0 <- isf |> mutate(target = val_of(t_oid[[tg]])[SampleID], target_blod = val_of(t_oid[[tg]], "below_lod")[SampleID])
  map(c("lesion site (AD)", "non-lesional / healthy skin", "all dISF"), \(st) {
    d <- if (st == "all dISF") d0 else d0 |> filter(stratum == st)
    imap(subsets(d), \(sel, nm) {
      dd <- d[sel, ]
      if (!nrow(dd)) return(NULL)
      above <- !dd$anchor_blod %in% TRUE & !dd$target_blod %in% TRUE
      r <- spearman_n(dd$anchor, dd$target, min_n)
      ra <- spearman_n(dd$anchor[above], dd$target[above], min_n)
      tibble(target = tg, site = st, subset = nm,
             subset_type = case_when(nm == "pooled (all groups)" ~ "a pooled", str_detect(nm, "^visit") ~ "c per visit", TRUE ~ "b within group"),
             n = r$n, rho = r$rho, p = r$p, n_subjects = n_distinct(dd$SubjectID[!is.na(dd$anchor) & !is.na(dd$target)]),
             n_both_above_LOD = ra$n, rho_both_above_LOD = ra$rho, p_both_above_LOD = ra$p)
    }) |> bind_rows()
  }) |> bind_rows()
}) |> bind_rows() |> filter(n > 0)
# pooled-only flag: pooled p < 0.05 but no within-group correlation with p < 0.05 in the same direction
flag <- targeted |> group_by(target, site) |>
  summarise(pooled_p = p[subset == "pooled (all groups)"][1], pooled_rho = rho[subset == "pooled (all groups)"][1],
            # within groups = subsets smaller than the pooled set (in the lesion site, "AD" equals the pooled set)
            within_hit = any(subset_type == "b within group" & n < n[subset == "pooled (all groups)"][1] &
                               p < 0.05 & sign(rho) == sign(pooled_rho[1]), na.rm = TRUE), .groups = "drop") |>
  mutate(interpretation = case_when(is.na(pooled_p) ~ "not estimable",
                                    pooled_p < 0.05 & !within_hit ~ "pooled only - likely driven by group differences, not by a protein-protein relationship",
                                    pooled_p < 0.05 ~ "pooled and within groups - consistent relationship",
                                    TRUE ~ "no pooled correlation (p >= 0.05)"))
targeted <- targeted |> left_join(flag |> select(target, site, interpretation), by = c("target", "site")) |>
  arrange(match(target, names(t_oid)), site, subset_type)
save_csv(targeted, cfg, out("targeted_correlations.csv"))
print(flag |> filter(target %in% mains) |> as.data.frame())

# ---- 2. proteome-wide ---------------------------------------------------------------------------------------------
pw_strata <- list(`all dISF` = isf, `lesion site (AD)` = isf |> filter(stratum == "lesion site (AD)"),
                  `non-lesional / healthy skin` = isf |> filter(stratum == "non-lesional / healthy skin"))
det_isf <- read_csv(file.path(cfg$paths$output, "qc", "assay_detection.csv"), show_col_types = FALSE) |>
  filter(matrix == "ISF") |> select(OlinkID, pct_above_LOD_dISF = frac_detected_all) |> mutate(pct_above_LOD_dISF = round(100 * pct_above_LOD_dISF, 1))
proteome <- imap(pw_strata, \(d, st) {
  d <- d |> filter(!is.na(anchor), SampleID %in% colnames(wide$ISF)) |> as.data.frame()
  if (nrow(d) < 8) return(NULL)
  e <- wide$ISF[setdiff(rownames(wide$ISF), a_oid), d$SampleID, drop = FALSE]
  covars <- c("state_g", "visit_f", "plate")
  msg("Proteome-wide (%s): %d samples, %d proteins", st, nrow(d), nrow(e))
  ps <- partial_spearman(d$anchor, e, d, covars, min_n = 8)
  # tryCatch: a rank-deficient design (e.g. plate identical to visit in a stratum) must not stop the step;
  # the partial Spearman ranking (primary) is then still reported
  mm <- tryCatch(fit_contrasts(e, d |> mutate(anchor_z = anchor),
                      as.formula(paste("~ anchor_z +", paste(covars[vapply(covars, \(v) n_distinct(d[[v]]) > 1, logical(1))], collapse = " + "),
                                       "+ (1 | SubjectID)")),
                      c(TNFRSF9 = "anchor_z"), paste("proteome-wide", st), cfg$stats$min_group_n),
                 error = \(e) { msg("  mixed model (%s) not fitted: %s", st, conditionMessage(e)); NULL })
  ps |> mutate(fdr_partial = p.adjust(p_partial, "BH")) |>
    left_join(if (is.null(mm)) tibble(OlinkID = character()) else
                mm |> transmute(OlinkID, slope_mixed = logFC, t_mixed = t, p_mixed = P.Value, fdr_mixed = adj.P.Val, method_mixed = method),
              by = "OlinkID") |>
    mutate(site = st, Assay = assay_name[OlinkID], .before = 1)
}) |> bind_rows() |>
  left_join(det_isf, by = "OlinkID") |>
  group_by(site) |> arrange(p_partial, .by_group = TRUE) |> mutate(rank = row_number(), of = n()) |> ungroup() |>
  mutate(target = OlinkID %in% t_oid, main_target = OlinkID %in% t_oid[intersect(mains, names(t_oid))],
         significant = coalesce(fdr_partial < fdr, FALSE)) |>
  relocate(site, rank, of, Assay, OlinkID, rho_partial, p_partial, fdr_partial, n)
save_csv(proteome, cfg, out("proteome_wide_correlation.csv"))
where_targets <- proteome |> filter(target) |> select(site, Assay, rank, of, rho_partial, p_partial, fdr_partial, slope_mixed, fdr_mixed, pct_above_LOD_dISF)
msg("Proteome-wide: %s", paste(proteome |> group_by(site) |> summarise(t = sprintf("%s %d of %d FDR < %g", site[1], sum(significant), n(), fdr)) |> pull(t), collapse = "; "))
print(where_targets |> filter(Assay %in% assay_name[t_oid[intersect(mains, names(t_oid))]]) |> as.data.frame())

# ---- 3. longitudinal (within-patient) ---------------------------------------------------------------------------------
longit <- map(names(t_oid), \(tg) {
  d0 <- isf |> filter(grp == "AD") |> mutate(target = val_of(t_oid[[tg]])[SampleID])
  map(c("lesion site", "non-lesional skin"), \(st) {
    d <- d0 |> filter(site_lab == st) |> arrange(SubjectID, visit_num)
    rm <- rmcorr_n(d$anchor, d$target, d$SubjectID, min_n)
    dd <- d |> group_by(SubjectID) |> mutate(d_anchor = anchor - lag(anchor), d_target = target - lag(target)) |> ungroup()
    dl <- spearman_n(dd$d_anchor, dd$d_target, min_n)
    tibble(target = tg, site = st, n_samples = nrow(d), subjects = n_distinct(d$SubjectID),
           r_within_subject = rm$r_within, p_within_subject = rm$p_within, n_within = rm$n_rm,
           delta_rho = dl$rho, delta_p = dl$p, n_visit_changes = dl$n)
  }) |> bind_rows()
}) |> bind_rows()
save_csv(longit, cfg, out("longitudinal_within_patient.csv"))

# ---- 4. cross-compartment (MicroAD only; RELAD / RELAD2 / LEIP not used) -------------------------------------------------
ser_vals <- clean |> filter(matrix == "Serum", cohort == "MicroAD", group %in% c("AD", "HC"), OlinkID %in% c(a_oid, t_oid)) |>
  select(SampleID, SubjectID, visit, group, relapse, OlinkID, value, below_lod)
ser_anchor <- ser_vals |> filter(OlinkID == a_oid) |> select(SubjectID, visit, serum_anchor = value)
cross <- isf |> filter(cohort == "MicroAD") |> inner_join(ser_anchor, by = c("SubjectID", "visit")) |>
  (\(d) map(c("lesion site", "non-lesional skin", "healthy skin", "all"), \(st) {
    x <- if (st == "all") d else d |> filter(site_lab == st)
    bind_cols(tibble(comparison = sprintf("dISF %s vs serum %s, same subject and visit", anchor, anchor), site = st),
              spearman_n(x$anchor, x$serum_anchor, min_n), rmcorr_n(x$anchor, x$serum_anchor, x$SubjectID, min_n))
  }) |> bind_rows())()
ser_w <- ser_vals |> mutate(protein = assay_name[OlinkID]) |> select(SampleID, SubjectID, group, relapse, protein, value) |>
  pivot_wider(names_from = protein, values_from = value)
ser_det <- ser_vals |> mutate(protein = assay_name[OlinkID]) |> group_by(protein) |>
  summarise(serum_pct_above_LOD = round(100 * mean(!below_lod, na.rm = TRUE), 1), .groups = "drop")
serum_targeted <- map(names(t_oid), \(tg) {
  tn <- assay_name[[t_oid[[tg]]]]; an <- assay_name[[a_oid]]
  if (!all(c(tn, an) %in% names(ser_w))) return(NULL)
  map(c("pooled", "AD", "healthy"), \(g) {
    x <- if (g == "pooled") ser_w else ser_w |> filter(recode(group, HC = "healthy") == g)
    bind_cols(tibble(target = tg, group = g), spearman_n(x[[an]], x[[tn]], min_n), rmcorr_n(x[[an]], x[[tn]], x$SubjectID, min_n))
  }) |> bind_rows()
}) |> bind_rows() |> mutate(protein = unname(assay_name[t_oid[target]])) |>
  left_join(ser_det |> rename(target_serum_pct_above_LOD = serum_pct_above_LOD), by = "protein") |> select(-protein)
save_csv(cross, cfg, out("cross_compartment.csv")); save_csv(serum_targeted, cfg, out("serum_correlations.csv"))

# ---- figures -------------------------------------------------------------------------------------------------------------
grp_cols <- c(`AD relapse` = "firebrick", `AD non-relapse` = "steelblue", `AD (relapse unknown)` = "grey40", healthy = "darkgreen")
scatter_data <- \(tg) isf |> mutate(target = val_of(t_oid[[tg]])[SampleID], target_blod = val_of(t_oid[[tg]], "below_lod")[SampleID],
                                    colour = case_when(grp == "healthy" ~ "healthy", rel %in% "relapse" ~ "AD relapse",
                                                       rel %in% "non-relapse" ~ "AD non-relapse", TRUE ~ "AD (relapse unknown)"),
                                    lod = if_else(anchor_blod %in% TRUE | target_blod %in% TRUE, "below LOD (either protein)", "both above LOD"),
                                    site_lab = factor(site_lab, levels = c("lesion site", "non-lesional skin", "healthy skin"))) |>
  filter(!is.na(anchor), !is.na(target))
lab_r <- \(tg) targeted |> filter(target == tg, subset == "pooled (all groups)", site != "all dISF") |>
  mutate(site_lab = if_else(site == "lesion site (AD)", "lesion site", "non-lesional skin"),
         label = sprintf("pooled rho = %.2f, p = %s, n = %d", rho, format.pval(p, digits = 2), n))
report_plots <- list()
for (tg in intersect(c(mains, names(t_oid)), names(t_oid))) {
  d <- scatter_data(tg)
  d$panel <- if_else(d$site_lab == "lesion site", "lesion site (AD)", "non-lesional / healthy skin")
  lr <- lab_r(tg) |> mutate(panel = site)
  p <- ggplot(d, aes(anchor, target)) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, colour = "grey50", linewidth = 0.4, linetype = 2) +
    geom_point(aes(colour = colour, shape = lod), size = 1.8, alpha = 0.85) +
    geom_text(data = lr, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5, size = 3, inherit.aes = FALSE) +
    scale_colour_manual(values = grp_cols) + scale_shape_manual(values = c(`both above LOD` = 16, `below LOD (either protein)` = 1)) +
    facet_wrap(~panel, scales = "free") +
    labs(title = sprintf("%s vs %s in dISF", anchor, assay_name[[t_oid[[tg]]]]),
         subtitle = sprintf("%% above LOD: %s %.0f%%, %s %.0f%%; open symbols = below LOD; dashed = pooled trend (see targeted_correlations.csv for within-group rho)",
                            anchor, overall_det$pct_above_LOD[overall_det$protein == assay_name[[a_oid]]],
                            assay_name[[t_oid[[tg]]]], overall_det$pct_above_LOD[overall_det$protein == assay_name[[t_oid[[tg]]]]]),
         x = sprintf("%s (NPX)", anchor), y = sprintf("%s (NPX)", assay_name[[t_oid[[tg]]]]), colour = NULL, shape = NULL) +
    theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
  save_plot(p, cfg, out("scatter", sprintf("%s_vs_%s.png", anchor, tg)), width = 10, height = 5.5)
  if (tg %in% mains) report_plots[[paste("scatter", tg)]] <- p
}
hm <- targeted |> filter(subset_type != "c per visit", site != "all dISF") |>
  mutate(subset = factor(subset, levels = unique(subset)), target = factor(target, levels = rev(names(t_oid))))
p <- ggplot(hm, aes(subset, target, fill = rho)) + geom_tile() +
  geom_text(aes(label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")), size = 3.5) +
  scale_fill_gradient2(low = "steelblue", high = "firebrick", limits = c(-1, 1)) + facet_wrap(~site) +
  labs(title = sprintf("%s vs target proteins in dISF: Spearman rho pooled and within groups", anchor),
       subtitle = "* p < 0.05, ** < 0.01, *** < 0.001 (unadjusted); per-visit results in the workbook", x = NULL, y = NULL, fill = "rho") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p, cfg, out("targeted_heatmap.png"), width = 11, height = 2.5 + 0.35 * length(t_oid))
report_plots$heatmap <- p
pw <- proteome |> filter(!is.na(rho_partial))
p <- ggplot(pw, aes(rho_partial, -log10(p_partial))) +
  geom_point(aes(colour = significant), size = 0.7, alpha = 0.6) +
  geom_point(data = \(x) filter(x, target), colour = "black", size = 1.6) +
  geom_text(data = \(x) filter(x, target), aes(label = Assay), size = 2.8, vjust = -0.7) +
  geom_text(data = \(x) x |> filter(significant, !target) |> group_by(site) |> slice_min(p_partial, n = 8), aes(label = Assay),
            size = 2.4, vjust = 1.4, colour = "grey30", check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"), labels = c(`FALSE` = "n.s.", `TRUE` = sprintf("FDR < %g", fdr))) +
  facet_wrap(~site) +
  labs(title = sprintf("Proteome-wide correlation with %s in dISF (adjusted for skin state, visit, plate)", anchor),
       subtitle = "black = target proteins", x = "partial Spearman rho", y = "-log10 p", colour = NULL) + theme(legend.position = "bottom")
save_plot(p, cfg, out("proteome_wide_correlation.png"), width = 13, height = 5.5)
report_plots$proteome <- p
saveRDS(list(plots = report_plots, notes = det_notes, flag = flag, where = where_targets, longit = longit, cross = cross,
             serum = serum_targeted, lod = lod_summary, anchor = anchor, mains = mains, not_on_panel = not_on_panel,
             top = proteome |> filter(significant) |> group_by(site) |> slice_min(p_partial, n = 15) |> ungroup()),
        out_path(cfg, out("summary.rds")))

# ---- workbook --------------------------------------------------------------------------------------------------------------
readme <- tibble(item = c("question", "anchor", "targets", "not on panel", "samples", "detectability and below-LOD handling",
                          "targeted (1)", "proteome-wide (2)", "longitudinal (3)", "cross-compartment (5)"),
                 note = c("Which measured proteins correlate with TNFRSF9 (4-1BB / CD137) in dISF, and do the target proteins?",
                          sprintf("%s (%s)", anchor, a_oid), paste(names(t_oid), collapse = ", "),
                          if (length(not_on_panel)) paste(not_on_panel, collapse = ", ") else "-",
                          sprintf("dISF of AD patients (lesion site, non-lesional skin) and healthy volunteers; CPUO excluded. %d samples.", nrow(isf)),
                          paste(det_notes, collapse = " "),
                          "Spearman rho, p and n: pooled, within groups (AD, healthy, relapse, non-relapse, lesional, ex-lesional) and per visit, per site. 'interpretation' flags correlations seen only when pooling. Pooled correlations use repeated samples of the same patients (not independent): use the within-group and longitudinal results to judge them.",
                          sprintf("Partial Spearman correlation adjusted for skin state / group, visit and plate (ranked, BH FDR per site); mixed model protein ~ %s + state + visit + plate + (1|subject) as a second method (slope = NPX change per NPX of %s). Proteins passing the dISF detection filter.", anchor, anchor),
                          "AD patients: repeated-measures correlation (subject-centred values) and Spearman correlation of visit-to-visit changes (delta-delta).",
                          "MicroAD only (RELAD, RELAD2 and LEIP excluded): dISF vs serum TNFRSF9 at the same subject and visit; TNFRSF9 vs targets in MicroAD serum."))
writexl::write_xlsx(list(README = readme, LOD_summary = lod_summary, targeted_correlations = targeted, pooled_vs_within = flag,
                         targets_in_ranking = where_targets, proteome_wide = proteome, longitudinal = longit,
                         cross_compartment = cross, serum_correlations = serum_targeted),
                    out_path(cfg, out("TNFRSF9_correlation.xlsx")))
msg("TNFRSF9 correlation analysis: %s", file.path(cfg$paths$output, "tnfrsf9_correlation"))
