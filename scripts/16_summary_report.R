# 16 - Executive summary PDF of all findings (step 16; step 17 only exports data)
# Collects the key numbers, tables and figures of steps 01-15 into output/Executive_summary.pdf,
# plus Executive_summary_tables.xlsx with the full list behind every shortened list in the PDF.
# Every section is optional: if an earlier step did not run, its page says so instead of failing.

source("R/utils.R")
source("R/models.R")
source("R/report.R")
cfg  <- load_config()
fdr  <- cfg$stats$fdr
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean_path <- file.path(cfg$paths$output, "data", "npx_clean.rds")
clean <- if (file.exists(clean_path)) readRDS(clean_path) else NULL

# study numbers used in the text (computed, so the text stays right when the data change)
rel_pat <- meta |> filter(cohort == "MicroAD", group == "AD", relapse %in% c("relapse", "non-relapse")) |> distinct(SubjectID, relapse)
n_rel <- sum(rel_pat$relapse == "relapse"); n_non <- sum(rel_pat$relapse == "non-relapse")
rel_txt <- sprintf("%d relapsing vs %d non-relapsing patients in MicroAD", n_rel, n_non)
rel_plates <- meta |> filter(matrix == "ISF", SubjectID %in% rel_pat$SubjectID[rel_pat$relapse == "relapse"]) |> distinct(plate) |> pull(plate)
per_visit_n <- meta |> filter(matrix == "ISF", cohort == "MicroAD", group == "AD", !is.na(visit_num)) |>
  distinct(visit_num, SubjectID) |> count(visit_num) |> pull(n)

# ---- gather results -------------------------------------------------------------------------------------
sq      <- out_csv(cfg, "qc", "sample_qc.csv")
lod     <- out_csv(cfg, "qc", "lod.csv")
det     <- out_csv(cfg, "qc", "assay_detection.csv")
cv      <- out_csv(cfg, "qc", "sample_control_cv.csv")
flags   <- out_csv(cfg, "metadata", "data_flags.csv")
isf     <- out_csv(cfg, "models", "ISF_results.csv")
isf_sum <- out_csv(cfg, "models", "ISF_summary.csv")
ser     <- out_csv(cfg, "models", "Serum_results.csv")
ser_sum <- out_csv(cfg, "models", "Serum_summary.csv")
delta   <- out_csv(cfg, "models", "ISF_relapse_delta_results.csv")
agree   <- out_csv(cfg, "models", "Serum_AD_vs_controls_agreement.csv")
vis     <- out_csv(cfg, "models", "ISF_by_visit_results.csv")
vis_sum <- out_csv(cfg, "models", "ISF_by_visit_summary.csv")
cons    <- out_csv(cfg, "visit_course", "consistency_across_visits.csv")
lstate  <- out_csv(cfg, "visit_course", "lesion_site_state_per_visit.csv")
corr    <- out_csv(cfg, "isf_serum", "isf_serum_correlation.csv")
gsea    <- out_csv(cfg, "enrichment", "gsea_results.csv")
leip    <- out_csv(cfg, "leip_reference", "leip_reference_summary.csv")
leip_a  <- out_csv(cfg, "leip_reference", "leip_clinical_associations.csv")
prof    <- out_csv(cfg, "isf_profile", "isf_detection_profile.csv")
det_mx  <- out_csv(cfg, "matrix_comparison", "detection_by_matrix.csv")
enr     <- out_csv(cfg, "matrix_comparison", "relative_enrichment.csv")
conc_s  <- out_csv(cfg, "matrix_comparison", "disease_signal_concordance_summary.csv")
traj_s  <- out_csv(cfg, "trajectories", "trajectory_summary.csv")
fov     <- out_csv(cfg, "focus", "focus_overview.csv")
kq_ans  <- out_csv(cfg, "key_questions", "answers.csv")
kq_q1   <- out_csv(cfg, "key_questions", "Q1_tests.csv")
kq_auc  <- out_csv(cfg, "key_questions", "Q4_auc.csv")
ov15    <- out_csv(cfg, "serum_vs_disf", "overlap_summary.csv")

n_sig <- \(r, mdl, ct) if (is.null(r)) NA else sum(r$significant[r$model == mdl & r$contrast == ct], na.rm = TRUE)
sig_names <- \(r, mdl, ct, dir, n = 10) {
  if (is.null(r)) return("n/a")
  d <- r |> filter(model == mdl, contrast == ct, significant, if (dir == "up") logFC > 0 else logFC < 0) |> arrange(P.Value)
  if (!nrow(d)) "none" else paste(head(d$Assay, n), collapse = ", ")
}
na_txt <- \(x) if (is.null(x) || length(x) == 0 || all(is.na(x))) "n/a" else x

out_file <- out_path(cfg, "Executive_summary.pdf")
pdf_open(out_file)

# ---- 1. title & data ------------------------------------------------------------------------------------------
section("Overview", {
  n_mx <- meta |> count(matrix, cohort, group)
  qc_line <- if (!is.null(sq)) sprintf("%d samples measured; %d failed Olink QC, %d with warnings, %d flagged as outliers (flagged, not removed).",
                                       nrow(sq), sum(sq$SampleQC == "FAIL"), sum(sq$SampleQC == "WARN"), sum(sq$outlier, na.rm = TRUE)) else "n/a"
  lod_line <- if (!is.null(lod)) paste(names(table(lod$LOD_source)), table(lod$LOD_source), sep = ": ", collapse = "; ") else "n/a"
  kept <- if (!is.null(det)) det |> group_by(matrix) |> summarise(n = n(), kept = sum(keep)) else NULL
  page_text("O-MicroAD - Olink Explore HT: executive summary",
    c("## Study material",
      sprintf("%d samples: %s.", nrow(meta), paste(sprintf("%s %d", names(table(meta$matrix)), table(meta$matrix)), collapse = ", ")),
      paste(sprintf("%s %s %s: %d", n_mx$matrix, n_mx$cohort, n_mx$group, n_mx$n), collapse = "; "),
      "## Data and quality control",
      qc_line,
      paste("LOD source per assay -", lod_line),
      if (!is.null(kept)) paste(sprintf("%s: %d of %d proteins pass the detection filter (>= %d%% of samples above LOD in at least one group)",
                                        kept$matrix, kept$kept, kept$n, round(100 * cfg$qc$min_detect_frac)), collapse = "; ") else "n/a",
      if (!is.null(cv)) sprintf("Sample controls: median inter-plate CV %.1f%%, intra-plate CV %.1f%%.",
                                100 * median(cv$inter_cv, na.rm = TRUE), 100 * median(cv$intra_cv, na.rm = TRUE)) else "n/a",
      sprintf("Metadata issues flagged (not corrected): %s.", if (is.null(flags) || !nrow(flags)) "none" else
        paste(sprintf("%s (%d)", names(table(flags$issue)), table(flags$issue)), collapse = "; ")),
      "## How to read this summary",
      sprintf("Proteome-wide results use the Benjamini-Hochberg false discovery rate (FDR < %g) within each comparison. Effects are differences in PC-normalised NPX (log2 scale). Focus proteins (e.g. CD137) are pre-specified and judged by their own p-value. Relapse analyses are exploratory (%s).", fdr, rel_txt),
      paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M"), "from", normalizePath(cfg$paths$output))),
    subtitle = "dermal interstitial fluid (dISF) and serum - automatically generated from the pipeline results")
})

# ---- 1b. key questions (step 15) ------------------------------------------------------------------------------------
section("Key questions", {
  if (is.null(kq_ans)) stop("step 15 results not found")
  items <- character()
  for (q in unique(kq_ans$question)) {
    a <- kq_ans |> filter(question == q)
    items <- c(items, paste("##", q), sprintf("%s: %s.", a$item, if_else(str_detect(a$verdict, "^dISF"), a$verdict,
                                                        paste0(str_to_upper(substr(a$verdict, 1, 1)), substring(a$verdict, 2)))))
  }
  page_text("Key questions - answers", c(items,
    "## Note", paste0("Single pre-specified tests (p < 0.05, not corrected across questions). Relapse results are exploratory: ", rel_txt, "; RELAD/RELAD2 serum provides larger groups. Evidence for every answer: next pages and key_questions/key_questions.xlsx.")),
    subtitle = "automatically derived from the tests; please interpret with the evidence tables", size = 9.5)
  page_table("Key questions - evidence", kq_ans |> transmute(question = str_extract(question, "^Q[0-9]"), item, verdict = str_trunc(verdict, 60)),
             rows_per_page = 30, note = "Full evidence text: key_questions/answers.csv")
  defn <- kq_ans |> filter(str_detect(item, "mast cell score - definition|markers used"))
  if (nrow(defn)) page_text("How the mast cell score is computed", c(
    paste(defn$verdict, defn$evidence, sep = ": "),
    "Why the score and single markers can disagree: the score averages the markers. If one marker (e.g. a tryptase) differs between groups but the others do not, or move in the opposite direction, the single marker can be significant while the score is not. The evidence pages that follow show each marker and the score separately."))
  ev_path <- file.path(cfg$paths$output, "key_questions", "evidence_plots.rds")
  if (file.exists(ev_path)) for (pl in readRDS(ev_path)) page_plot(pl)
  if (!is.null(kq_q1)) {
    ents <- unique(kq_q1$entity)
    p <- ggplot(kq_q1, aes(label, factor(entity, levels = rev(ents)), fill = estimate)) + geom_tile() +
      geom_text(aes(label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")), size = 4) +
      scale_fill_gradient2(low = "steelblue", high = "firebrick") + facet_grid(~type, scales = "free_x", space = "free_x") +
      labs(title = "Q1  Mast cell markers: elevated in AD, or associated with relapse?",
           subtitle = "effect (log2 / score units); * p < 0.05, ** < 0.01, *** < 0.001", x = NULL, y = NULL) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    page_plot(p)
  }
  if (!is.null(kq_auc)) {
    p <- ggplot(kq_auc |> filter(!is.na(AUC)), aes(AUC, predictor, colour = if_else(str_detect(predictor, "^dISF"), "dISF", "serum"))) +
      geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
      geom_pointrange(aes(xmin = ci_low, xmax = ci_high)) + facet_wrap(~entity) + coord_cartesian(xlim = c(0, 1)) +
      scale_colour_manual(values = c(dISF = "firebrick", serum = "steelblue")) +
      labs(title = "Q4/Q5  Do values BEFORE the relapse separate relapsers from non-relapsers? (dISF vs serum)",
           subtitle = sprintf("AUC with 95%% bootstrap CI; 0.5 = no separation; MicroAD n = %d vs %d patients (CI unreliable), RELAD/RELAD2 larger", n_rel, n_non),
           x = "AUC", y = NULL, colour = NULL)
    page_plot(p)
  }
})

# ---- 2. key findings ---------------------------------------------------------------------------------------------
section("Key findings", {
  cs <- if (!is.null(cons)) cons |> filter(contrast == "Lsite_vs_HC") else NULL
  n_vis <- if (!is.null(cs) && nrow(cs)) max(cs$n_visits) else NA
  items <- c(
    "## Aim 1 - dISF proteome",
    if (!is.null(prof)) sprintf("%d proteins are detectable in dISF (>= %d%% of samples above LOD); %d robustly (>= 90%%). %d are detectable only in lesional AD skin: %s.",
                                sum(prof$frac_detected >= cfg$qc$min_detect_frac, na.rm = TRUE), round(100 * cfg$qc$min_detect_frac),
                                sum(prof$frac_detected >= 0.9, na.rm = TRUE), sum(prof$lesion_restricted, na.rm = TRUE),
                                paste(top_names(prof |> filter(lesion_restricted)), "... (full list: Executive_summary_tables.xlsx, sheet lesion_restricted)")) else "dISF profile: n/a",
    sprintf("Lesional vs non-lesional AD skin (all visits): %s proteins (up: %s ...; down: %s ...; full list: sheet dISF_L_vs_NL).",
            na_txt(n_sig(isf, "states_all_visits", "AD_L_vs_NL")), sig_names(isf, "states_all_visits", "AD_L_vs_NL", "up", 8),
            sig_names(isf, "states_all_visits", "AD_L_vs_NL", "down", 5)),
    sprintf("Ex-lesional vs non-lesional (residual signature after clearing): %s proteins (sheet dISF_xL_vs_NL); lesional vs ex-lesional: %s (sheet dISF_L_vs_xL).",
            na_txt(n_sig(isf, "states_all_visits", "AD_xL_vs_NL")), na_txt(n_sig(isf, "states_all_visits", "AD_L_vs_xL"))),
    sprintf("Lesional AD skin vs healthy skin: %s proteins (sheet dISF_L_vs_HC); non-lesional AD vs healthy skin: %s (sheet dISF_NL_vs_HC).",
            na_txt(n_sig(isf, "states_all_visits", "AD_L_vs_HC")), na_txt(n_sig(isf, "states_all_visits", "AD_NL_vs_HC"))),
    "## Visit by visit (lesion site vs healthy skin)",
    if (!is.null(vis_sum)) paste(vis_sum |> filter(contrast == "Lsite_vs_HC") |>
                                   transmute(t = sprintf("%s: %d (n = %d patients + controls)", model, n_sig, n_subjects)) |> pull(t), collapse = "; ") else "n/a",
    if (!is.null(cs) && nrow(cs)) sprintf("Regulated at all %d visits, same direction: %d proteins at FDR < %g (%s), %d at p < 0.05 at every visit.",
                                          n_vis, sum(cs$all_visits_fdr), fdr, top_names(cs |> filter(all_visits_fdr), 12), sum(cs$all_visits_nominal)) else "n/a",
    "## Aim 2 - dISF vs blood",
    if (!is.null(det_mx)) sprintf("Detected in both matrices: %d; dISF only: %d; serum only: %d.",
                                  sum(det_mx$detected_in == "both"), sum(det_mx$detected_in == "dISF only"), sum(det_mx$detected_in == "serum only")) else "n/a",
    if (!is.null(enr)) paste(map_chr(unique(enr$model), \(md) sprintf("Relatively enriched in dISF vs serum (>= 2-fold vs the typical protein), %s: %d, e.g. %s.",
                               str_remove(md, "relative enrichment, "), sum(enr$direction == "enriched in dISF" & enr$model == md),
                               top_names(enr |> filter(direction == "enriched in dISF", model == md) |> arrange(desc(rel_log2_isf_vs_serum)), 6))),
                             collapse = " ") else "n/a",
    if (!is.null(corr)) sprintf("dISF and serum levels move together within patients over visits for %d (lesion site) and %d (non-lesional skin) proteins; between-patient correlation: %d / %d.",
                                sum(corr$fdr_within < fdr & corr$site == "L", na.rm = TRUE), sum(corr$fdr_within < fdr & corr$site == "NL", na.rm = TRUE),
                                sum(corr$fdr_between < fdr & corr$site == "L", na.rm = TRUE), sum(corr$fdr_between < fdr & corr$site == "NL", na.rm = TRUE)) else "n/a",
    sprintf("Serum AD vs in-study healthy controls: %s proteins; vs LEIP biobank: %s; significant against both (same direction): %s (%s).",
            na_txt(n_sig(ser, "AD_vs_HC_in_study", "AD_vs_HC")), na_txt(n_sig(ser, "AD_vs_Biobank", "AD_vs_Biobank")),
            if (is.null(agree)) "n/a" else sum(agree$agree_both_controls, na.rm = TRUE),
            if (is.null(agree)) "n/a" else top_names(agree |> filter(agree_both_controls), 10)),
    "## Aim 3 - disease course and relapse (exploratory)",
    sprintf("Relapse, cleared skin (xL - NL): %s proteins; serum MicroAD relapse: %s; RELAD/RELAD2 serum relapse: %s.",
            if (is.null(delta)) "n/a" else sum(delta$significant), na_txt(n_sig(ser, "MicroAD_relapse", "relapse_vs_non")),
            na_txt(n_sig(ser, "RELAD_relapse", "relapse_vs_non"))),
    if (!is.null(traj_s)) paste(sprintf("%s: %d", traj_s$model, traj_s$n_sig), collapse = "; ") else "Trajectories: n/a",
    "## LEIP population reference",
    if (!is.null(leip)) sprintf("%d dISF-serum correlated proteins checked in LEIP; %d associate with clinical parameters (FDR < %g); %d differ between LEIP and in-study healthy controls (possible pre-analytical effect).",
                                nrow(leip), if (is.null(leip_a)) 0 else n_distinct(leip_a$Assay[leip_a$significant]), fdr,
                                sum(leip$source_shift, na.rm = TRUE)) else "n/a")
  # focus proteins (one line each) and serum vs dISF overlap
  if (!is.null(fov)) {
    items <- c(items, "## Focus proteins (single-protein tests, p-values uncorrected)")
    for (pr in unique(fov$label)) {
      d <- fov |> filter(label == pr)
      g <- \(cmp) { r <- d |> filter(comparison == cmp); if (!nrow(r)) "n/a" else sprintf("%+.2f (p = %.2g)", r$estimate[1], r$p[1]) }
      items <- c(items, sprintf("%s: dISF lesional vs non-lesional %s; ex-lesional vs non-lesional %s; non-lesional vs healthy %s; serum AD vs healthy %s.",
                                pr, g("dISF: states_all_visits AD_L_vs_NL"), g("dISF: states_all_visits AD_xL_vs_NL"),
                                g("dISF: states_all_visits AD_NL_vs_HC"), g("serum: AD_vs_HC_in_study AD_vs_HC")))
    }
    missing_fp <- setdiff(toupper(unlist(cfg$focus_proteins)), toupper(unique(fov$protein)))
    if (length(missing_fp)) items <- c(items, sprintf("Not measured in this dataset: %s.", paste(missing_fp, collapse = ", ")))
  }
  if (!is.null(ov15)) {
    o <- ov15 |> filter(tier == "nominal", visit == "all visits")
    items <- c(items, "## Serum vs dISF (same question in both matrices, all visits pooled, p < 0.05)",
               sprintf("%s, %s: dISF %d, serum %d, both %d (same direction %d); dISF only %d (+%d not measurable in serum); serum only %d.",
                       o$question, o$isf_site, o$dISF_significant, o$serum_significant, o$both_same + o$both_opposite, o$both_same,
                       o$dISF_only, o$dISF_only_not_in_serum, o$serum_only))
  }
  page_text("Key findings", items, subtitle = sprintf("automatically extracted; FDR < %g unless stated", fdr), size = 9.5)
})

# ---- 3. overview table of all comparisons -----------------------------------------------------------------------
section("All comparisons", {
  tab <- bind_rows(isf_sum |> mutate(matrix = "dISF"), vis_sum |> mutate(matrix = "dISF per visit"),
                   ser_sum |> mutate(matrix = "serum")) |>
    select(matrix, model, contrast, samples = n_samples, subjects = n_subjects, proteins = n_assays, significant = n_sig, up = n_up, down = n_down)
  page_table("All comparisons: number of significant proteins", tab, rows_per_page = 28,
             note = sprintf("FDR < %g within each comparison. Model details: models/*_summary.csv", fdr))
})

# ---- 4. visit course ----------------------------------------------------------------------------------------------
section("Visit course", {
  if (is.null(vis)) stop("step 13 results not found")
  cnt <- vis |> group_by(model, contrast) |> summarise(up = sum(significant & logFC > 0), down = sum(significant & logFC < 0), .groups = "drop") |>
    pivot_longer(c(up, down), names_to = "dir", values_to = "n")
  p <- ggplot(cnt, aes(model, if_else(dir == "up", n, -n), fill = dir)) + geom_col() + geom_hline(yintercept = 0) +
    facet_wrap(~contrast) + scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
    labs(title = "dISF per visit: significant proteins", subtitle = "Lsite = tracked lesion (lesional at V1 / relapse, ex-lesional after clearing)",
         x = NULL, y = "proteins (down < 0 < up)", fill = NULL)
  page_plot(p)
  if (!is.null(lstate)) page_table("State of the tracked lesion site per visit", lstate)
  for (ct in c("Lsite_vs_HC", "Lsite_vs_NL")) {
    ct_lab <- c(Lsite_vs_HC = "lesion site vs healthy skin", Lsite_vs_NL = "lesion site vs non-lesional skin")[[ct]]
    sheet <- c(Lsite_vs_HC = "all_visits_Lsite_vs_HC", Lsite_vs_NL = "all_visits_Lsite_vs_NL")[[ct]]
    cs <- cons |> filter(contrast == ct)
    sel <- cs |> filter(all_visits_fdr); tier <- sprintf("FDR < %g at every visit", fdr)
    if (!nrow(sel)) { sel <- cs |> filter(all_visits_nominal); tier <- "p < 0.05 at every visit" }
    if (!nrow(sel)) { page_text(sprintf("Proteins regulated at all visits (%s)", ct_lab), "No protein was significant at every visit."); next }
    top <- sel |> slice_min(max_p, n = 16, with_ties = FALSE)
    d <- vis |> filter(contrast == ct, OlinkID %in% top$OlinkID) |>
      mutate(se = abs(logFC / t), Assay = factor(Assay, levels = unique(top$Assay)))
    p <- ggplot(d, aes(model, logFC, group = 1)) + geom_hline(yintercept = 0, colour = "grey60") +
      geom_errorbar(aes(ymin = logFC - 1.96 * se, ymax = logFC + 1.96 * se), width = 0.2, colour = "grey50") +
      geom_line(colour = "grey40") + geom_point(aes(colour = significant), size = 1.8) +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), labels = c(`FALSE` = "n.s.", `TRUE` = "FDR sig.")) +
      facet_wrap(~Assay, scales = "free_y") +
      labs(title = sprintf("Time course: %s, proteins regulated at all visits (%s)", ct_lab, tier),
           subtitle = sprintf("top %d of %d; full list and per-visit data: Executive_summary_tables.xlsx, sheet %s", nrow(top), nrow(sel), sheet),
           x = NULL, y = "difference in NPX (log2, 95% CI)", colour = NULL)
    page_plot(p)
    page_table(sprintf("Proteins regulated at all visits (%s, %s)", ct_lab, tier),
               sel |> slice_min(max_p, n = 60, with_ties = FALSE) |>
                 select(Assay, direction, mean_logFC, visits = n_visits, `FDR-sig visits` = n_fdr, max_p),
               note = sprintf("Full list: Executive_summary_tables.xlsx, sheet %s", sheet))
  }
})

# ---- 5. volcano plots, all visits ---------------------------------------------------------------------------------
section("dISF volcano plots", {
  if (is.null(isf)) stop("step 04 results not found")
  page_plot(volcano(isf |> filter(model == "states_all_visits",
                                  contrast %in% c("AD_L_vs_NL", "AD_xL_vs_NL", "AD_L_vs_HC", "AD_NL_vs_HC")),
                    "dISF, all visits: skin states (dashed line = FDR cutoff of each panel)", fdr))
})

# ---- 6. pathways ---------------------------------------------------------------------------------------------------
section("Pathways", {
  if (is.null(gsea)) stop("step 07 results not found")
  g <- gsea |> filter(model == "states_all_visits", contrast == "AD_L_vs_NL", padj < fdr) |>
    arrange(padj) |> head(25) |> transmute(pathway = str_trunc(pathway, 60), NES = NES, padj = padj, size = size)
  page_table("Gene sets enriched in lesional vs non-lesional dISF (GSEA)", g,
             note = "positive NES = higher in lesional skin; all contrasts: enrichment/gsea_results.csv")
})

# ---- 7. dISF vs serum ------------------------------------------------------------------------------------------------
section("dISF vs serum", {
  if (!is.null(det_mx)) {
    p <- ggplot(det_mx |> count(detected_in), aes(reorder(detected_in, -n), n)) + geom_col(fill = "grey45") +
      geom_text(aes(label = n), vjust = -0.3) +
      labs(title = "Where are proteins detectable? (MicroAD samples)", x = NULL, y = "proteins")
    page_plot(p)
  }
  if (!is.null(enr)) for (md in unique(enr$model))
    page_table(sprintf("Proteins most enriched in dISF relative to serum - %s", str_remove(md, "relative enrichment, ")),
               enr |> filter(model == md, direction == "enriched in dISF") |>
                 arrange(desc(rel_log2_isf_vs_serum)) |> head(25) |>
                 select(Assay, rel_log2_isf_vs_serum, adj.P.Val, n_samples, n_subjects),
               note = "relative log2 ratio centred on the typical protein; candidates for local (skin) production. Full lists: Executive_summary_tables.xlsx")
  if (!is.null(corr)) page_table("Proteins whose dISF and serum levels move together within patients",
                                 corr |> filter(fdr_within < fdr) |> arrange(fdr_within) |> head(25) |>
                                   select(Assay, site, n_pairs, r_within, fdr_within, r_between, fdr_between))
  if (!is.null(conc_s)) page_table("Do skin disease signals appear in blood?", conc_s)
})

# ---- 8. serum --------------------------------------------------------------------------------------------------------
section("Serum", {
  if (!is.null(agree)) page_table("Serum: AD vs controls, significant against in-study AND biobank controls",
                                  agree |> filter(agree_both_controls) |> select(Assay, logFC_AD_vs_HC, adj.P.Val_AD_vs_HC,
                                                                                  logFC_AD_vs_Biobank, adj.P.Val_AD_vs_Biobank))
})

# ---- 9. serum vs dISF per visit (step 14) ---------------------------------------------------------------------------
section("Serum vs dISF", {
  if (is.null(ov15)) stop("step 14 results not found")
  cat_cols <- c(`both, same direction` = "purple3", `both, opposite direction` = "orange3", `dISF only` = "firebrick",
                `dISF only - not measurable in serum` = "darkred", `serum only` = "steelblue")
  for (tr in c("nominal", "FDR")) {
    b <- ov15 |> filter(tier == tr) |>
      mutate(visit = factor(visit, levels = c(paste0("V", 1:12), "all visits"))) |>
      select(question, visit, isf_site, `both, same direction` = both_same, `both, opposite direction` = both_opposite,
             `dISF only` = dISF_only, `dISF only - not measurable in serum` = dISF_only_not_in_serum, `serum only` = serum_only) |>
      pivot_longer(-c(question, visit, isf_site), names_to = "cat", values_to = "n") |>
      mutate(cat = factor(cat, levels = names(cat_cols)))
    p <- ggplot(b, aes(visit, n, fill = cat)) + geom_col() + facet_grid(question ~ isf_site, scales = "free_y") +
      scale_fill_manual(values = cat_cols) +
      labs(title = sprintf("Serum vs dISF: overlap of significant proteins per visit (%s)",
                           if (tr == "FDR") sprintf("FDR < %g", fdr) else "p < 0.05, exploratory"),
           subtitle = "red = information only dISF provides; purple/orange = seen in both; blue = only serum",
           x = NULL, y = "proteins", fill = NULL) + theme(legend.position = "bottom")
    page_plot(p)
  }
  page_table("Serum vs dISF overlap (p < 0.05)", ov15 |> filter(tier == "nominal") |>
               select(question, visit, isf_site, dISF = dISF_significant, serum = serum_significant, both_same, both_opposite,
                      dISF_only, `dISF only, not in serum` = dISF_only_not_in_serum, serum_only),
             note = "Venn diagrams, volcano plots and protein lists: serum_vs_disf/", rows_per_page = 32)
})

# ---- 10. focus proteins -------------------------------------------------------------------------------------------------
section("Focus proteins", {
  if (is.null(fov)) stop("step 12 results not found")
  hm <- fov |> mutate(comparison = factor(comparison, levels = unique(comparison)))
  p <- ggplot(hm, aes(comparison, label, fill = estimate)) + geom_tile() +
    geom_text(aes(label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")), size = 4) +
    scale_fill_gradient2(low = "steelblue", high = "firebrick") +
    labs(title = "Focus proteins: effect (log2) per comparison", subtitle = "single-protein tests: * p < 0.05, ** < 0.01, *** < 0.001",
         x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  page_plot(p)
  if (!is.null(clean)) {
    fps <- unique(fov$protein)
    d <- clean |> filter(Assay %in% fps, matrix == "ISF") |>
      mutate(cond = case_when(group == "HC" ~ "healthy skin", group == "AD" & state == "non-lesional" ~ "AD non-lesional",
                              group == "AD" & state == "ex-lesional" ~ "AD ex-lesional", group == "AD" & state == "lesional" ~ "AD lesional",
                              group == "CPUO" ~ paste("CPUO", state)),
             cond = factor(cond, levels = c("healthy skin", "AD non-lesional", "AD ex-lesional", "AD lesional", "CPUO non-lesional", "CPUO lesional")),
             protein = fov$label[match(Assay, fov$protein)])
    p <- ggplot(d |> filter(!is.na(cond)), aes(cond, value)) + geom_boxplot(outlier.shape = NA, fill = "grey92") +
      geom_jitter(width = 0.15, size = 0.5) + facet_wrap(~protein, scales = "free_y") +
      labs(title = "Focus proteins in dISF by skin state", x = NULL, y = "NPX") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))
    page_plot(p)
  }
  page_table("Focus proteins: key single-protein tests",
             fov |> select(protein, comparison, estimate, ci_low, ci_high, p, proteome_wide_FDR), rows_per_page = 32,
             note = "p = single-protein test (primary for pre-specified proteins). Per-protein reports and figures: focus/<protein>/")
})

# ---- 11. methods & caveats ---------------------------------------------------------------------------------------------
section("Methods", {
  page_text("Methods, caveats and where to find everything", c(
    "## Methods",
    sprintf("NPX: %s (plate-control normalised; plates were not matrix-randomised). LOD: Olink fixed LOD file where available (per-sample for count-based assays), else Olink's negative-control method.", cfg$npx_column),
    "Per-protein models: limma/dream with empirical Bayes moderation; subject as random effect for repeated samples, plate as fixed effect in dISF; cohort and plate in serum. One-term mixed models: limma + duplicateCorrelation.",
    "dISF vs serum: repeated-measures correlation (within patients) and Spearman on patient means (between patients); relative enrichment = centred paired log2 difference.",
    "## Caveats",
    if (length(rel_plates) == 1) sprintf("All %d MicroAD relapsers' dISF samples are on %s (relapse and plate cannot be separated): dISF relapse results are exploratory.", n_rel, rel_plates)
    else sprintf("Relapse analyses compare %s: exploratory.", rel_txt),
    sprintf("Per-visit and V1 comparisons have %s AD patients per visit: absence of significance is not absence of an effect (strict FDR cutoff with few samples).",
            if (length(per_visit_n)) paste(unique(range(per_visit_n)), collapse = "-") else "few"),
    "LEIP biobank serum differs pre-analytically; AD vs biobank differences are only trusted when they agree with the in-study controls.",
    "Sex is available for MicroAD (manifest) and LEIP, age only for LEIP; the models are not adjusted for them.",
    "## Output folders (all under the output directory)",
    "qc/ - QC tables and plots | models/ - all model results and volcano plots (models/volcano/) | visit_course/ - per-visit analysis and time courses",
    "isf_profile/ - dISF proteome | matrix_comparison/ - dISF vs serum | serum_vs_disf/ - overlap per visit (Venn) | isf_serum/ - correlations",
    "trajectories/ - disease course | leip_reference/ | enrichment/ | focus/ - focus proteins (overview + one folder per protein)"),
    size = 10)
})

grDevices::dev.off()

# ---- companion workbook: the full lists behind the summary ---------------------------------------------------------------
sig_list <- \(r, mdl, ct) if (is.null(r)) NULL else r |> filter(model == mdl, contrast == ct, significant) |>
  transmute(Assay, OlinkID, direction = if_else(logFC > 0, "up", "down"), logFC, P.Value, FDR = adj.P.Val, n_samples, n_subjects) |>
  arrange(P.Value)
per_visit_wide <- \(ct) if (is.null(vis)) NULL else vis |> filter(contrast == ct) |>
  select(Assay, OlinkID, visit = model, logFC, FDR = adj.P.Val) |>
  pivot_wider(names_from = visit, values_from = c(logFC, FDR), names_glue = "{visit}_{.value}")
all_vis <- \(ct) if (is.null(cons)) NULL else cons |> filter(contrast == ct, all_visits_fdr | all_visits_nominal) |>
  transmute(Assay, OlinkID, tier = if_else(all_visits_fdr, "FDR at every visit", "p < 0.05 at every visit"), direction, mean_logFC,
            visits = n_visits, fdr_sig_visits = n_fdr, max_p) |>
  left_join(per_visit_wide(ct) |> select(-Assay), by = "OlinkID") |> arrange(desc(tier == "FDR at every visit"), max_p)
enr_sheets <- if (is.null(enr)) list() else
  map(split(enr, enr$model), \(d) d |> filter(direction != "not different") |>
        transmute(Assay, OlinkID, direction, rel_log2_isf_vs_serum, P.Value, FDR = adj.P.Val, n_samples, n_subjects) |>
        arrange(desc(rel_log2_isf_vs_serum))) |>
  setNames(paste0("dISF_vs_serum_", c(`relative enrichment, AD lesional skin` = "lesional", `relative enrichment, AD ex-lesional skin` = "ex_lesional",
                                      `relative enrichment, AD non-lesional skin` = "non_lesional", `relative enrichment, healthy skin` = "healthy")[names(split(enr, enr$model))]))
tabs <- c(list(
  lesion_restricted = if (!is.null(prof)) prof |> filter(lesion_restricted) |> select(-any_of("lesion_restricted")) else NULL,
  dISF_L_vs_NL = sig_list(isf, "states_all_visits", "AD_L_vs_NL"),
  dISF_xL_vs_NL = sig_list(isf, "states_all_visits", "AD_xL_vs_NL"),
  dISF_L_vs_xL = sig_list(isf, "states_all_visits", "AD_L_vs_xL"),
  dISF_L_vs_HC = sig_list(isf, "states_all_visits", "AD_L_vs_HC"),
  dISF_NL_vs_HC = sig_list(isf, "states_all_visits", "AD_NL_vs_HC"),
  dISF_V1_L_vs_NL = sig_list(isf, "baseline_V1", "AD_L_vs_NL"),
  dISF_V1_L_vs_HC = sig_list(isf, "baseline_V1", "AD_L_vs_HC"),
  all_visits_Lsite_vs_HC = all_vis("Lsite_vs_HC"),
  all_visits_Lsite_vs_NL = all_vis("Lsite_vs_NL"),
  all_visits_NL_vs_HC = all_vis("NL_vs_HC"),
  per_visit_Lsite_vs_NL_all = per_visit_wide("Lsite_vs_NL"),
  per_visit_Lsite_vs_HC_all = per_visit_wide("Lsite_vs_HC")),
  enr_sheets,
  list(serum_AD_vs_HC = sig_list(ser, "AD_vs_HC_in_study", "AD_vs_HC"),
       serum_AD_both_controls = if (!is.null(agree)) agree |> filter(agree_both_controls) else NULL,
       key_questions = kq_ans, focus_proteins = fov)) |> compact()
tabs <- tabs[vapply(tabs, nrow, 1L) > 0]
index <- tibble(sheet = names(tabs), rows = vapply(tabs, nrow, 1L),
                content = c(lesion_restricted = "proteins detectable only in lesional AD dISF (step 09)",
                            dISF_L_vs_NL = sprintf("significant: lesional vs non-lesional AD skin, all visits (FDR < %g)", fdr),
                            dISF_xL_vs_NL = "significant: ex-lesional vs non-lesional, all visits",
                            dISF_L_vs_xL = "significant: lesional vs ex-lesional, all visits",
                            dISF_L_vs_HC = "significant: lesional AD skin vs healthy skin, all visits",
                            dISF_NL_vs_HC = "significant: non-lesional AD skin vs healthy skin, all visits",
                            dISF_V1_L_vs_NL = "significant: lesional vs non-lesional at V1",
                            dISF_V1_L_vs_HC = "significant: lesional vs healthy at V1",
                            all_visits_Lsite_vs_HC = "regulated at every visit, lesion site vs healthy skin, with per-visit logFC/FDR",
                            all_visits_Lsite_vs_NL = "regulated at every visit, lesion site vs non-lesional skin, with per-visit logFC/FDR",
                            all_visits_NL_vs_HC = "regulated at every visit, non-lesional vs healthy skin",
                            per_visit_Lsite_vs_NL_all = "ALL proteins: lesion site vs non-lesional per visit (logFC, FDR)",
                            per_visit_Lsite_vs_HC_all = "ALL proteins: lesion site vs healthy skin per visit (logFC, FDR)",
                            dISF_vs_serum_lesional = "relatively enriched in dISF (or serum) - AD lesional skin vs serum",
                            dISF_vs_serum_ex_lesional = "relatively enriched in dISF (or serum) - AD ex-lesional skin vs serum",
                            dISF_vs_serum_non_lesional = "relatively enriched in dISF (or serum) - AD non-lesional skin vs serum",
                            dISF_vs_serum_healthy = "relatively enriched in dISF (or serum) - healthy skin vs serum",
                            serum_AD_vs_HC = "significant: serum AD vs in-study healthy controls",
                            serum_AD_both_controls = "serum AD vs controls: significant against in-study AND biobank controls",
                            key_questions = "answers and evidence for the key questions",
                            focus_proteins = "focus proteins: key single-protein tests")[names(tabs)])
writexl::write_xlsx(c(list(index = index), tabs), out_path(cfg, "Executive_summary_tables.xlsx"))
msg("Full lists: %s", out_path(cfg, "Executive_summary_tables.xlsx"))
msg("Executive summary: %s", out_file)
