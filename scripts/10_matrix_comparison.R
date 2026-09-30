# 10 - Aim 2: compare the dISF and blood (serum) proteome, MicroAD
#   a) detection: which proteins are measurable in dISF only, serum only, or both
#   b) relative enrichment in dISF: paired ISF - serum NPX difference per subject-visit, centred on
#      the median difference of all proteins in that pair. Both matrices are normalised to the same
#      plate control, so the difference is a relative log2 ratio; centring removes the overall
#      dilution of ISF. Proteins far above the typical protein are candidates for local (skin)
#      production. Only proteins detected in both matrices are tested; 'enriched' needs FDR < fdr
#      and a relative difference of at least stats$min_rel_log2 (default 1 = 2-fold).
#   c) disease signals: are skin effects (dISF models, step 04) mirrored in serum (step 05)?
#      Serum side: MicroAD only (AD_vs_HC_MicroAD) - RELAD / RELAD2 serum is not used for dISF comparisons.
#      And, visit by visit, does the lesional-minus-non-lesional ISF difference track serum?
# Correlation of ISF and serum levels over time is in step 06.
# Out: output/matrix_comparison/*

source("R/utils.R")
source("R/models.R")
cfg   <- load_config()
clear_outputs(cfg, "matrix_comparison")
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
min_f <- cfg$qc$min_detect_frac
min_rel <- cfg$stats$min_rel_log2 %||% 1   # minimum relative difference (log2) to call a protein enriched
assay_map <- clean |> distinct(OlinkID, Assay)

micro <- clean |> filter(cohort == "MicroAD")

# ---- a) detection ------------------------------------------------------------------------------
det <- micro |>
  group_by(OlinkID, Assay, matrix) |>
  summarise(frac = mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = matrix, values_from = frac, names_prefix = "detected_") |>
  mutate(detected_in = case_when(detected_ISF >= min_f & detected_Serum >= min_f ~ "both",
                                 detected_ISF >= min_f ~ "dISF only",
                                 detected_Serum >= min_f ~ "serum only",
                                 TRUE ~ "neither"))
save_csv(det, cfg, "matrix_comparison", "detection_by_matrix.csv")
msg("Detection (MicroAD): %s", paste(names(table(det$detected_in)), table(det$detected_in), collapse = ", "))
p <- ggplot(det, aes(detected_Serum, detected_ISF, colour = detected_in)) + geom_point(size = 0.8, alpha = 0.7) +
  geom_hline(yintercept = min_f, linetype = 2) + geom_vline(xintercept = min_f, linetype = 2) +
  labs(title = "Share of samples above LOD: dISF vs serum", x = "serum", y = "dISF", colour = NULL)
save_plot(p, cfg, "matrix_comparison", "detection_isf_vs_serum.png", width = 7, height = 6)

# ---- b) relative enrichment in dISF ----------------------------------------------------------------
both <- det |> filter(detected_in == "both") |> pull(OlinkID) |>
  intersect(rownames(wide$ISF)) |> intersect(rownames(wide$Serum))
mk <- meta |> filter(cohort == "MicroAD")
serum_ids <- mk |> filter(matrix == "Serum", SampleID %in% colnames(wide$Serum)) |> select(SubjectID, visit, serum_id = SampleID)
pairs <- mk |> filter(matrix == "ISF", SampleID %in% colnames(wide$ISF)) |>
  transmute(SubjectID, visit, isf_id = SampleID, group, plate,
            site = case_when(site == "L" ~ "lesional site", TRUE ~ "non-lesional / healthy skin")) |>
  inner_join(serum_ids, by = c("SubjectID", "visit")) |>
  mutate(SampleID = paste(isf_id, serum_id, sep = "_"), ones = 1)

enrich <- NULL
if (length(both) >= 10 && nrow(pairs) >= 4) {
  delta <- wide$ISF[both, pairs$isf_id, drop = FALSE] - wide$Serum[both, pairs$serum_id, drop = FALSE]
  delta <- sweep(delta, 2, apply(delta, 2, median, na.rm = TRUE))   # centre each pair
  colnames(delta) <- pairs$SampleID
  enrich <- map(unique(pairs$site), \(st) {
    pp <- pairs |> filter(site == st)
    fit_contrasts(delta, pp, ~ 0 + ones + (1 | SubjectID), c(isf_vs_serum = "ones"),
                  paste("relative enrichment,", st), cfg$stats$min_group_n)
  }) |> bind_rows()
}
if (!is.null(enrich) && nrow(enrich)) {
  enrich <- enrich |>
    annotate_results(assay_map, cfg$stats$fdr) |>
    rename(rel_log2_isf_vs_serum = logFC) |>
    mutate(direction = case_when(significant & rel_log2_isf_vs_serum >= min_rel ~ "enriched in dISF",
                                 significant & rel_log2_isf_vs_serum <= -min_rel ~ "enriched in serum",
                                 TRUE ~ "not different"))
  save_csv(enrich, cfg, "matrix_comparison", "relative_enrichment.csv")
  print(count(enrich, model, direction))
  p <- enrich |> group_by(model) |> mutate(rank = rank(-rel_log2_isf_vs_serum)) |> ungroup() |>
    ggplot(aes(rank, rel_log2_isf_vs_serum, colour = direction)) + geom_point(size = 0.7) +
    geom_text(data = \(d) d |> group_by(model) |> slice_max(rel_log2_isf_vs_serum, n = 10),
              aes(label = Assay), size = 2.5, hjust = -0.2, colour = "black", check_overlap = TRUE) +
    facet_wrap(~model) +
    scale_colour_manual(values = c(`enriched in dISF` = "firebrick", `enriched in serum` = "steelblue", `not different` = "grey60")) +
    labs(title = "Relative dISF / serum level (centred log2 ratio, detected in both)", x = "protein rank", y = "relative log2 ratio")
  save_plot(p, cfg, "matrix_comparison", "relative_enrichment.png", width = 11, height = 5)
} else {
  enrich <- NULL
  msg("Too few proteins detected in both matrices or too few matched pairs - enrichment skipped.")
}

# ---- c) disease signals: skin vs blood -----------------------------------------------------------------
isf_res   <- file.path(cfg$paths$output, "models", "ISF_results.csv")
serum_res <- file.path(cfg$paths$output, "models", "Serum_results.csv")
if (file.exists(isf_res) && file.exists(serum_res)) {
  ir <- read_csv(isf_res, show_col_types = FALSE)
  sr <- read_csv(serum_res, show_col_types = FALSE)
  pick <- \(r, mdl, ct, nm) r |> filter(model == mdl, contrast == ct) |>
    transmute(OlinkID, "{nm}_logFC" := logFC, "{nm}_fdr" := adj.P.Val)
  comparisons <- list(
    c(isf_model = "states_all_visits", isf_ct = "AD_L_vs_NL", serum_model = "MicroAD_active_vs_cleared",
      serum_ct = "active_vs_cleared", label = "lesion activity: skin (L vs NL) vs blood (active vs cleared visits)"),
    c(isf_model = "states_all_visits", isf_ct = "AD_L_vs_HC", serum_model = "AD_vs_HC_MicroAD",
      serum_ct = "AD_vs_HC", label = "disease: lesional skin vs healthy skin, AD vs healthy serum"),
    c(isf_model = "states_all_visits", isf_ct = "AD_NL_vs_HC", serum_model = "AD_vs_HC_MicroAD",
      serum_ct = "AD_vs_HC", label = "systemic: non-lesional AD skin vs healthy skin, AD vs healthy serum")
  )
  conc <- map(comparisons, \(cp) {
    a <- pick(ir, cp[["isf_model"]], cp[["isf_ct"]], "isf"); b <- pick(sr, cp[["serum_model"]], cp[["serum_ct"]], "serum")
    if (!nrow(a) || !nrow(b)) return(NULL)
    inner_join(a, b, by = "OlinkID") |> mutate(comparison = cp[["label"]]) |>
      mutate(category = case_when(isf_fdr < cfg$stats$fdr & serum_fdr < cfg$stats$fdr &
                                    sign(isf_logFC) == sign(serum_logFC) ~ "both, same direction",
                                  isf_fdr < cfg$stats$fdr & serum_fdr < cfg$stats$fdr ~ "both, opposite",
                                  isf_fdr < cfg$stats$fdr ~ "dISF only",
                                  serum_fdr < cfg$stats$fdr ~ "serum only", TRUE ~ "neither"))
  }) |> bind_rows()
  if (nrow(conc)) conc <- conc |> left_join(assay_map, by = "OlinkID") |> relocate(comparison, Assay, .after = OlinkID)
  if (nrow(conc)) {
    save_csv(conc, cfg, "matrix_comparison", "disease_signal_concordance.csv")
    cs <- conc |> group_by(comparison) |>
      summarise(n = n(), spearman_logFC = cor(isf_logFC, serum_logFC, method = "spearman", use = "complete.obs"),
                both_same = sum(category == "both, same direction"), isf_only = sum(category == "dISF only"),
                serum_only = sum(category == "serum only"), .groups = "drop")
    save_csv(cs, cfg, "matrix_comparison", "disease_signal_concordance_summary.csv")
    print(as.data.frame(cs))
    p <- ggplot(conc, aes(isf_logFC, serum_logFC, colour = category)) +
      geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
      geom_point(size = 0.8, alpha = 0.7) +
      geom_text(data = conc |> filter(category == "both, same direction") |> group_by(comparison) |>
                  slice_min(isf_fdr, n = 8), aes(label = Assay), size = 2.5, vjust = -0.6, colour = "black") +
      facet_wrap(~str_wrap(comparison, 45), scales = "free") +
      labs(title = "Disease effects in dISF vs serum", x = "dISF effect (log2)", y = "serum effect (log2)", colour = NULL) +
      theme(legend.position = "bottom")
    save_plot(p, cfg, "matrix_comparison", "disease_signal_concordance.png", width = 13, height = 5.5)
  }
} else msg("Model results from steps 04/05 not found - concordance skipped.")

# visit-level: does the lesional signal in skin (L - NL) track serum at the same visit?
lp <- mk |> filter(matrix == "ISF", group == "AD", SampleID %in% colnames(wide$ISF)) |>
  select(SubjectID, visit, site, SampleID) |>
  pivot_wider(names_from = site, values_from = SampleID) |>
  inner_join(serum_ids, by = c("SubjectID", "visit")) |> filter(!is.na(L), !is.na(NL))
shared <- intersect(rownames(wide$ISF), rownames(wide$Serum))
if (nrow(lp) >= 8 && length(shared)) {
  vt <- map(shared, \(a) {
    x <- wide$ISF[a, lp$L] - wide$ISF[a, lp$NL]; y <- wide$Serum[a, lp$serum_id]
    ok <- !is.na(x) & !is.na(y) & lp$SubjectID %in% names(which(table(lp$SubjectID) >= 2))
    if (sum(ok) < 8) return(NULL)
    r <- suppressWarnings(rmcorr::rmcorr(participant = subj, measure1 = x, measure2 = y,
                                         dataset = data.frame(subj = factor(lp$SubjectID[ok]), x = x[ok], y = y[ok])))
    tibble(OlinkID = a, n_visits = sum(ok), r_within = r$r, p = r$p)
  }) |> bind_rows()
  if (nrow(vt)) vt <- vt |> mutate(fdr = p.adjust(p, "BH")) |> left_join(assay_map, by = "OlinkID") |> arrange(p)
  if (nrow(vt)) {
    save_csv(vt, cfg, "matrix_comparison", "lesion_signal_vs_serum_by_visit.csv")
    msg("%d proteins: serum level tracks the lesional-minus-non-lesional skin difference within patients (FDR < %.2f)",
        sum(vt$fdr < cfg$stats$fdr), cfg$stats$fdr)
  } else vt <- NULL
}

writexl::write_xlsx(Filter(Negate(is.null), list(
  detection = det, relative_enrichment = enrich,
  disease_concordance = if (exists("conc") && nrow(conc)) conc else NULL,
  lesion_signal_vs_serum = if (exists("vt")) vt else NULL)),
  out_path(cfg, "matrix_comparison", "matrix_comparison.xlsx"))
