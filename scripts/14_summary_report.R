# 14 - Executive summary PDF of all findings (runs last)
# Collects the key numbers, tables and figures of steps 01-13 into output/Executive_summary.pdf.
# Every section is optional: if an earlier step did not run, its page says so instead of failing.

source("R/utils.R")
source("R/models.R")
source("R/report.R")
cfg  <- load_config()
fdr  <- cfg$stats$fdr
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean_path <- file.path(cfg$paths$output, "data", "npx_clean.rds")
clean <- if (file.exists(clean_path)) readRDS(clean_path) else NULL

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
      sprintf("Proteome-wide results use the Benjamini-Hochberg false discovery rate (FDR < %g) within each comparison. Effects are differences in PC-normalised NPX (log2 scale). Focus proteins (e.g. CD137) are pre-specified and judged by their own p-value. Relapse analyses are exploratory (4 relapsers in MicroAD).", fdr),
      paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M"), "from", normalizePath(cfg$paths$output))),
    subtitle = "dermal interstitial fluid (dISF) and serum - automatically generated from the pipeline results")
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
                                top_names(prof |> filter(lesion_restricted))) else "dISF profile: n/a",
    sprintf("Lesional vs non-lesional AD skin (all visits): %s proteins (up: %s ...; down: %s ...).",
            na_txt(n_sig(isf, "states_all_visits", "AD_L_vs_NL")), sig_names(isf, "states_all_visits", "AD_L_vs_NL", "up", 8),
            sig_names(isf, "states_all_visits", "AD_L_vs_NL", "down", 5)),
    sprintf("Ex-lesional vs non-lesional (residual signature after clearing): %s proteins; lesional vs ex-lesional: %s.",
            na_txt(n_sig(isf, "states_all_visits", "AD_xL_vs_NL")), na_txt(n_sig(isf, "states_all_visits", "AD_L_vs_xL"))),
    sprintf("Lesional AD skin vs healthy skin: %s proteins; non-lesional AD vs healthy skin: %s.",
            na_txt(n_sig(isf, "states_all_visits", "AD_L_vs_HC")), na_txt(n_sig(isf, "states_all_visits", "AD_NL_vs_HC"))),
    "## Visit by visit (lesion site vs healthy skin)",
    if (!is.null(vis_sum)) paste(vis_sum |> filter(contrast == "Lsite_vs_HC") |>
                                   transmute(t = sprintf("%s: %d (n = %d patients + controls)", model, n_sig, n_subjects)) |> pull(t), collapse = "; ") else "n/a",
    if (!is.null(cs) && nrow(cs)) sprintf("Regulated at all %d visits, same direction: %d proteins at FDR < %g (%s), %d at p < 0.05 at every visit.",
                                          n_vis, sum(cs$all_visits_fdr), fdr, top_names(cs |> filter(all_visits_fdr), 12), sum(cs$all_visits_nominal)) else "n/a",
    "## Aim 2 - dISF vs blood",
    if (!is.null(det_mx)) sprintf("Detected in both matrices: %d; dISF only: %d; serum only: %d.",
                                  sum(det_mx$detected_in == "both"), sum(det_mx$detected_in == "dISF only"), sum(det_mx$detected_in == "serum only")) else "n/a",
    if (!is.null(enr)) sprintf("Relatively enriched in dISF (>= 2-fold vs the typical protein, non-lesional skin): %d, e.g. %s.",
                               sum(enr$direction == "enriched in dISF" & str_detect(enr$model, "non-lesional")),
                               top_names(enr |> filter(direction == "enriched in dISF", str_detect(model, "non-lesional")) |> arrange(desc(rel_log2_isf_vs_serum)), 8)) else "n/a",
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
  # focus proteins
  for (fp in unlist(cfg$focus_proteins %||% character())) {
    f <- list.files(file.path(cfg$paths$output, "focus"), pattern = "_report\\.xlsx$", recursive = TRUE, full.names = TRUE)
    f <- f[str_detect(basename(f), fixed(fp))]
    if (!length(f)) next
    t <- readxl::read_excel(f[1], sheet = "prespecified_tests")
    key <- t |> filter(paste(model, contrast) %in% c("states_all_visits AD_L_vs_NL", "states_all_visits AD_xL_vs_NL",
                                                    "states_all_visits AD_NL_vs_HC", "AD_vs_HC_in_study AD_vs_HC",
                                                    "relapse_delta_xL_minus_NL relapse_vs_non"))
    items <- c(items, sprintf("## Focus protein %s%s", fp, if (fp == "TNFRSF9") " (CD137 / 4-1BB)" else ""),
               sprintf("%s %s: %+.2f log2 (95%% CI %.2f to %.2f), p = %.2g", key$matrix, key$contrast, key$estimate,
                       key$ci_low, key$ci_high, key$p))
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
  cs <- cons |> filter(contrast == "Lsite_vs_HC")
  sel <- cs |> filter(all_visits_fdr); tier <- sprintf("FDR < %g at every visit", fdr)
  if (!nrow(sel)) { sel <- cs |> filter(all_visits_nominal); tier <- "p < 0.05 at every visit" }
  if (nrow(sel)) {
    top <- sel |> slice_min(max_p, n = 16, with_ties = FALSE)
    d <- vis |> filter(contrast == "Lsite_vs_HC", OlinkID %in% top$OlinkID) |>
      mutate(se = abs(logFC / t), Assay = factor(Assay, levels = top$Assay))
    p <- ggplot(d, aes(model, logFC, group = 1)) + geom_hline(yintercept = 0, colour = "grey60") +
      geom_errorbar(aes(ymin = logFC - 1.96 * se, ymax = logFC + 1.96 * se), width = 0.2, colour = "grey50") +
      geom_line(colour = "grey40") + geom_point(aes(colour = significant), size = 1.8) +
      scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), labels = c(`FALSE` = "n.s.", `TRUE` = "FDR sig.")) +
      facet_wrap(~Assay, scales = "free_y") +
      labs(title = sprintf("Time course: lesion site vs healthy skin, proteins regulated at all visits (%s)", tier),
           subtitle = sprintf("top %d of %d; full list: visit_course/consistency_across_visits.csv", nrow(top), nrow(sel)),
           x = NULL, y = "difference in NPX (log2, 95% CI)", colour = NULL)
    page_plot(p)
    page_table(sprintf("Proteins regulated at all visits (lesion site vs healthy skin, %s)", tier),
               sel |> slice_min(max_p, n = 60, with_ties = FALSE) |>
                 select(Assay, direction, mean_logFC, visits = n_visits, `FDR-sig visits` = n_fdr, max_p))
  } else page_text("Proteins regulated at all visits", "No protein was significant at every visit.")
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
  if (!is.null(enr)) page_table("Proteins most enriched in dISF relative to serum (non-lesional / healthy skin)",
                                enr |> filter(str_detect(model, "non-lesional"), direction == "enriched in dISF") |>
                                  arrange(desc(rel_log2_isf_vs_serum)) |> head(25) |>
                                  select(Assay, rel_log2_isf_vs_serum, adj.P.Val),
                                note = "relative log2 ratio centred on the typical protein; candidates for local (skin) production")
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

# ---- 9. focus proteins -------------------------------------------------------------------------------------------------
section("Focus proteins", {
  for (fp in unlist(cfg$focus_proteins %||% character())) {
    f <- list.files(file.path(cfg$paths$output, "focus"), pattern = "_report\\.xlsx$", recursive = TRUE, full.names = TRUE)
    f <- f[str_detect(basename(f), fixed(fp))]
    if (!length(f)) next
    lab <- if (fp == "TNFRSF9") "TNFRSF9 (CD137 / 4-1BB)" else fp
    if (!is.null(clean)) {
      d <- clean |> filter(Assay == fp | OlinkID == fp, matrix == "ISF") |>
        mutate(cond = case_when(group == "HC" ~ "healthy skin", group == "AD" & state == "non-lesional" ~ "AD non-lesional",
                                group == "AD" & state == "ex-lesional" ~ "AD ex-lesional", group == "AD" & state == "lesional" ~ "AD lesional",
                                group == "CPUO" ~ paste("CPUO", state)),
               cond = factor(cond, levels = c("healthy skin", "AD non-lesional", "AD ex-lesional", "AD lesional", "CPUO non-lesional", "CPUO lesional")))
      if (nrow(d)) {
        p <- ggplot(d |> filter(!is.na(cond)), aes(cond, value)) + geom_boxplot(outlier.shape = NA, fill = "grey92") +
          geom_jitter(width = 0.15, size = 1) + labs(title = sprintf("%s in dISF by skin state", lab), x = NULL, y = "NPX") +
          theme(axis.text.x = element_text(angle = 20, hjust = 1))
        page_plot(p)
      }
    }
    t <- readxl::read_excel(f[1], sheet = "prespecified_tests")
    page_table(sprintf("%s: pre-specified tests", lab),
               t |> select(matrix, model, contrast, estimate, ci_low, ci_high, p, any_of("proteome_wide_FDR")),
               note = "p = single-protein test (primary for a pre-specified protein); proteome_wide_FDR for comparison. Full report: focus/")
  }
})

# ---- 10. methods & caveats ---------------------------------------------------------------------------------------------
section("Methods", {
  page_text("Methods, caveats and where to find everything", c(
    "## Methods",
    sprintf("NPX: %s (plate-control normalised; plates were not matrix-randomised). LOD: Olink fixed LOD file where available (per-sample for count-based assays), else Olink's negative-control method.", cfg$npx_column),
    "Per-protein models: limma/dream with empirical Bayes moderation; subject as random effect for repeated samples, plate as fixed effect in dISF; cohort and plate in serum. One-term mixed models: limma + duplicateCorrelation.",
    "dISF vs serum: repeated-measures correlation (within patients) and Spearman on patient means (between patients); relative enrichment = centred paired log2 difference.",
    "## Caveats",
    "All four MicroAD relapsers' dISF samples are on plate 1: dISF relapse results are exploratory.",
    "Per-visit and V1 comparisons have 6-11 patients per visit: absence of significance is not absence of an effect (strict FDR cutoff with few samples).",
    "LEIP biobank serum differs pre-analytically; AD vs biobank differences are only trusted when they agree with the in-study controls.",
    "Age and sex are only available for LEIP; serum models are not adjusted for them.",
    "## Output folders (all under the output directory)",
    "qc/ - QC tables and plots | models/ - all model results and volcano plots (models/volcano/) | visit_course/ - per-visit analysis and time courses",
    "isf_profile/ - dISF proteome | matrix_comparison/ - dISF vs serum | isf_serum/ - correlations | trajectories/ - disease course | leip_reference/ | enrichment/ | focus/ - focus proteins"),
    size = 10)
})

grDevices::dev.off()
msg("Executive summary: %s", out_file)
