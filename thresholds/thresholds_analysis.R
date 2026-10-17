# Thresholds (part of run_all.R; shortcut: run_thresholds.R): which thresholds remove data or decide results, and are they
# too strict for these data? Each one is measured on the data:
#   1 samples: Olink SampleQC FAIL (removed), WARN and outliers (flagged only), low volume (flagged only)
#   2 proteins: the detection (LOD) filter at other values - assays kept per matrix; the proteins of interest
#   3 results: five comparisons re-run with other detection filters (proteins tested, p < 0.05, FDR); do the proteins the
#     filter removes carry more p < 0.05 than chance?; gene sets (GSEA) and the signature gene set at other filters
#   4 group sizes per visit against stats$min_group_n and visit_min_subjects; the smallest possible exact p
#   5 relapse target criteria at V1: where the candidates are lost (p, exact p, LOD in the relapsers)
#   6 inventory: every threshold in the code - value, where, what it does
# In:  output of steps 01-02 (data/npx_clean.rds, qc/sample_qc.csv); MSigDB / KEGG gene sets
# Out: output/thresholds/ - Thresholds_summary.pdf, thresholds.xlsx, figures as PNG (16:9, 300 dpi, for PowerPoint) and PDF

source("R/utils.R")
source("R/qc.R")
source("R/qc_overview.R")
source("R/visit3.R")
source("R/relapse_v1.R")
source("R/relapse_followup.R")
source("R/thresholds.R")
cfg <- load_config()
th <- cfg$thresholds %||% list()
cur <- cfg$qc$min_detect_frac %||% 0.5
fr <- sort(unique(c(unlist(th$detect_fracs %||% c(0, 0.1, 0.25, 0.5, 0.75, 0.9)), cur)))
fr_gsea <- sort(unique(c(unlist(th$gsea_detect_fracs %||% c(0.25, 0.75)), cur)))
fdr <- cfg$stats$fdr %||% 0.05; fdr_lv <- sort(unique(c(unlist(th$fdr_levels %||% c(0.05, 0.1, 0.2)), fdr)))
k_cur <- cfg$qc$outlier_sd %||% 3; sds <- sort(unique(c(unlist(th$outlier_sds %||% c(2, 2.5, 3, 4)), k_cur)))
min_n <- cfg$stats$min_group_n %||% 3; vmin <- cfg$stats$visit_min_subjects %||% 5
gmin <- cfg$enrichment$min_size %||% 10; gmax <- cfg$enrichment$max_size %||% 500
sig_sets <- unlist(cfg$visits_v1_v4$signature_sets %||% "HALLMARK_MTORC1_SIGNALING")
colours <- cfg$qc_overview$colours %||% list(above = "#69005F", below = "#FF506E")
out <- "thresholds"
clear_outputs(cfg, out)
set.seed(1)
fmt_p <- \(x) ifelse(is.na(x), "-", formatC(x, digits = 2, format = "g"))
sig3 <- \(x) x |> mutate(across(where(is.double), \(v) signif(v, 3)))
mx_name <- c(ISF = "dISF", Serum = "serum")
manifest_txt <- sprintf("sample manifest: %s", basename(cfg$paths$manifest))

clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
info <- clean |> distinct(OlinkID, Assay, UniProt)

# =============================================================================================================================
# 1  Samples: what sample QC removes or only flags
# =============================================================================================================================
sq_path <- file.path(cfg$paths$output, "qc", "sample_qc.csv")
sq <- if (file.exists(sq_path)) read_csv(sq_path, show_col_types = FALSE) else NULL
beyond <- \(x, ok, k) abs(x - mean(x[ok], na.rm = TRUE)) > k * sd(x[ok], na.rm = TRUE)      # the step 02 outlier rule
lv <- unlist(cfg$qc$low_volume_ul %||% list(ISF = 20, Serum = 40))
vol <- meta |> filter(matrix %in% names(lv)) |> group_by(matrix) |>
  summarise(with_volume = sum(!is.na(volume_ul)), low_volume_flagged = sum(volume_ul < lv[matrix], na.rm = TRUE), .groups = "drop")
samples_qc <- if (is.null(sq)) tibble(matrix = character()) else sq |> group_by(matrix) |>
  summarise(samples = n(), removed_FAIL = sum(SampleQC == "FAIL"), WARN_kept = sum(SampleQC == "WARN"),
            outliers_flagged = sum(coalesce(outlier, FALSE) & SampleQC != "FAIL"), .groups = "drop")
samples_qc <- samples_qc |> full_join(vol, by = "matrix") |> mutate(matrix = mx_name[matrix])
outliers_by_sd <- if (is.null(sq)) tibble() else map(sds, \(k) sq |> group_by(matrix) |>
  summarise(outlier_sd = k, outliers = sum((beyond(median_npx, SampleQC != "FAIL", k) | beyond(iqr_npx, SampleQC != "FAIL", k)) & SampleQC != "FAIL", na.rm = TRUE),
            .groups = "drop")) |> bind_rows() |> mutate(matrix = mx_name[matrix], current = outlier_sd == k_cur)
removed_samples <- if (is.null(sq)) tibble() else sq |> filter(SampleQC == "FAIL") |> select(SampleID, matrix) |>
  left_join(meta |> select(SampleID, SubjectID, cohort, group, visit, site), by = "SampleID") |> mutate(matrix = mx_name[matrix])

# =============================================================================================================================
# 2  Proteins: the detection (LOD) filter at other values
# =============================================================================================================================
det <- assay_detection(clean, cur)       # per assay and matrix: largest share above LOD in a detection group, share of all samples
det_tbl <- map(fr, \(f) det |> group_by(matrix) |>
  summarise(detection_threshold = f, assays = n(), assays_without_LOD = sum(!lod_available),
            kept_in_one_group = sum(!lod_available | max_group_frac >= f), kept_of_all_samples = sum(!lod_available | frac_detected_all >= f),
            .groups = "drop")) |> bind_rows() |>
  mutate(matrix = mx_name[matrix], current = detection_threshold == cur,
         pct_kept_in_one_group = round(100 * kept_in_one_group / assays, 1), pct_kept_of_all_samples = round(100 * kept_of_all_samples / assays, 1)) |>
  arrange(matrix, detection_threshold)
det_groups <- clean |> distinct(matrix, det_group, SampleID) |> count(matrix, det_group, name = "samples") |>
  mutate(matrix = mx_name[matrix], samples_needed_above_LOD = ceiling(cur * samples))

# proteins of interest (the protein lists in config.yml): share above LOD per detection group, kept by step 02
rf <- cfg$relapse_followup %||% list()
poi_lists <- list(`focus proteins` = cfg$focus_proteins, `key questions` = c(cfg$key_questions$proteins, cfg$key_questions$mast_cell_markers),
                  `leukocytes / T cells` = c(rf$leukocyte_markers, rf$leukocyte_support, rf$t_cell_markers, rf$t_cell_activation),
                  cytokines = rf$cytokines, `tracked (V1-V4)` = cfg$visits_v1_v4$track_proteins,
                  `TNFRSF9 targets` = c(cfg$tnfrsf9_correlation$anchor, cfg$tnfrsf9_correlation$targets), extra = th$extra_proteins)
poi <- imap(poi_lists, \(g, nm) tibble(list = nm, gene = unlist(g))) |> bind_rows() |> filter(!is.na(gene)) |> distinct(gene, .keep_all = TRUE)
poi_ids <- assay_ids(info, poi$gene)
poi_missing <- setdiff(poi$gene, poi_ids$gene)
poi_by_group <- clean |> filter(OlinkID %in% poi_ids$OlinkID) |> group_by(matrix, OlinkID, det_group) |>
  summarise(pct_above_LOD = 100 * mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  left_join(poi_ids, by = "OlinkID") |> left_join(poi, by = "gene") |>
  left_join(det |> select(matrix, OlinkID, keep), by = c("matrix", "OlinkID")) |> mutate(matrix = mx_name[matrix])

# =============================================================================================================================
# 3  Results: five comparisons re-run with other detection filters
# =============================================================================================================================
ad_isf <- clean |> filter(cohort == "MicroAD", group == "AD", matrix == "ISF")
healthy <- \(mx) clean |> filter(cohort == "MicroAD", group == "HC", matrix == mx) |> distinct(SampleID, SubjectID)
pairs_at <- function(v) {
  x <- ad_isf |> filter(site %in% c("L", "NL"), visit_num == v) |> distinct(SampleID, SubjectID, site)
  x |> filter(SubjectID %in% (x |> count(SubjectID) |> filter(n == 2) |> pull(SubjectID))) |>
    mutate(grp = factor(if_else(site == "L", "lesion site", "non-lesional"), c("non-lesional", "lesion site")))
}
two_groups <- \(a, b, la, lb) bind_rows(a |> mutate(grp = la), b |> mutate(grp = lb)) |> mutate(grp = factor(grp, c(lb, la)))
cmps <- list(
  relapse_V1 = list(label = "Relapse vs non-relapse, dISF lesion site V1", short = "relapse V1", paired = FALSE,
                    s = ad_isf |> filter(site == "L", visit_num == 1, relapse %in% c("relapse", "non-relapse")) |> distinct(SampleID, SubjectID, relapse) |>
                      mutate(grp = factor(relapse, c("non-relapse", "relapse")))),
  lesional_V1 = list(label = "Lesion site vs non-lesional, dISF V1 (paired)", short = "lesional V1", paired = TRUE, s = pairs_at(1)),
  lesional_V3 = list(label = "Lesion site vs non-lesional, dISF V3 (paired)", short = "lesional V3", paired = TRUE, s = pairs_at(3)),
  AD_vs_healthy_dISF = list(label = "AD lesion site V1 vs healthy skin, dISF", short = "AD vs healthy dISF", paired = FALSE,
                            s = two_groups(ad_isf |> filter(site == "L", visit_num == 1) |> distinct(SampleID, SubjectID), healthy("ISF"), "AD lesion site", "healthy")),
  AD_vs_healthy_serum = list(label = "AD V1 vs healthy, serum (MicroAD)", short = "AD vs healthy serum", paired = FALSE,
                             s = two_groups(clean |> filter(cohort == "MicroAD", group == "AD", matrix == "Serum", visit_num == 1) |> distinct(SampleID, SubjectID),
                                            healthy("Serum"), "AD", "healthy")))
count_res <- \(r) tibble(proteins_tested = nrow(r), p_below_0.01 = sum(r$p < 0.01, na.rm = TRUE), p_below_0.05 = sum(r$p < 0.05, na.rm = TRUE),
                         !!!set_names(map(fdr_lv, \(q) sum(r$FDR < q, na.rm = TRUE)), sprintf("FDR_below_%g", fdr_lv)))
runs <- imap(cmps, \(cm, nm) {
  n <- table(cm$s$grp); n_txt <- paste(rev(sprintf("%s %d", names(n), as.integer(n))), collapse = " vs ")
  if (!nrow(cm$s) || min(n) < min_n) return(list(label = cm$label, short = cm$short, n_txt = n_txt, skipped = sprintf("not run: fewer than %d samples in a group (%s)", min_n, n_txt)))
  res <- map(set_names(fr, as.character(fr)), \(f) de_limma(clean, cm$s, paired = cm$paired, min_f = f, fdr = fdr))
  tab <- imap(res, \(r, f) count_res(r) |> mutate(detection_threshold = as.numeric(f), .before = 1)) |> bind_rows()
  if (nm == "relapse_V1") {               # how many p < 0.05 random labels give, at each filter (every way to choose the relapse group)
    tab$chance_median_p_below_0.05 <- map_dbl(res, \(r) {
      M <- patient_matrix(clean, cm$s, r$OlinkID); M <- M[r$OlinkID, sort(colnames(M)), drop = FALSE]
      rc <- relabel_check(M, cm$s$grp[match(colnames(M), cm$s$SubjectID)] == "relapse")
      median(colSums(rc$up_hits) + colSums(rc$down_hits))
    })
  }
  msg("%s (%s): %s proteins tested at detection filters %s", cm$label, n_txt, paste(tab$proteins_tested, collapse = " / "), paste(fr, collapse = " / "))
  list(label = cm$label, short = cm$short, n_txt = n_txt, paired = cm$paired, s = cm$s, res = res, tab = tab)
})
runs_ok <- keep(runs, \(x) is.null(x$skipped))
cmp_tbl <- imap(runs, \(x, nm) if (!is.null(x$skipped)) tibble(comparison = x$label, short = x$short, samples = x$n_txt, note = x$skipped) else
  x$tab |> mutate(comparison = x$label, short = x$short, samples = x$n_txt, current = detection_threshold == cur, .before = 1)) |> bind_rows()

# the proteins the filter removes (tested at the laxest filter, not at the current one): more p < 0.05 than chance?
f0 <- min(fr); f_strict <- fr[fr > cur][1]
removed <- imap(runs_ok, \(x, nm) {
  r0 <- x$res[[as.character(f0)]]; rc <- x$res[[as.character(cur)]]
  rm_ <- r0 |> filter(!OlinkID %in% rc$OlinkID); kp <- r0 |> filter(OlinkID %in% rc$OlinkID)
  lost <- if (is.na(f_strict)) tibble() else rc |> filter(p < 0.05, !OlinkID %in% x$res[[as.character(f_strict)]]$OlinkID)
  list(row = tibble(comparison = x$label, short = x$short, removed_by_filter = nrow(rm_), removed_p_below_0.05 = sum(rm_$p < 0.05, na.rm = TRUE),
                    removed_expected_by_chance = 0.05 * nrow(rm_), kept = nrow(kp), kept_p_below_0.05 = sum(kp$p < 0.05, na.rm = TRUE),
                    kept_expected_by_chance = 0.05 * nrow(kp),
                    FDR_hits_current = sum(rc$FDR < fdr, na.rm = TRUE), FDR_hits_without_filter = sum(r0$FDR < fdr, na.rm = TRUE),
                    p_below_0.05_lost_with_stricter_filter = if (is.na(f_strict)) NA_integer_ else nrow(lost)),
       hits = rm_ |> filter(p < 0.05) |> mutate(comparison = x$label, .before = 1) |> select(comparison, Assay, OlinkID, logFC, p, FDR, starts_with("pct_above_LOD_"), starts_with("mean_NPX_")))
})
removed_tbl <- map(removed, "row") |> bind_rows() |>
  mutate(removed_ratio_to_chance = round(removed_p_below_0.05 / pmax(removed_expected_by_chance, 1e-9), 2),
         kept_ratio_to_chance = round(kept_p_below_0.05 / pmax(kept_expected_by_chance, 1e-9), 2))
removed_hits <- map(removed, "hits") |> bind_rows()

# gene sets: GSEA at three detection filters, and the signature gene set(s)
sets <- v3_gene_sets(cfg)
gsea_cmp <- intersect(unlist(th$gsea_comparisons %||% c("relapse_V1", "lesional_V3")), names(runs_ok))
gsea_runs <- map(set_names(gsea_cmp), \(nm) map(set_names(fr_gsea, as.character(fr_gsea)), \(f) gsea_keyed(runs_ok[[nm]]$res[[as.character(f)]], sets, gmin, gmax)))
gsea_tbl <- imap(gsea_runs, \(gl, nm) imap(gl, \(g, f) tibble(comparison = runs_ok[[nm]]$label, detection_threshold = as.numeric(f), gene_sets_tested = nrow(g),
  FDR_below_0.05 = sum(g$FDR < 0.05), FDR_below_0.1 = sum(g$FDR < 0.1), FDR_below_0.25 = sum(g$FDR < 0.25),
  FDR_below_0.05_higher = sum(g$FDR < 0.05 & g$NES > 0), FDR_below_0.05_lower = sum(g$FDR < 0.05 & g$NES < 0))) |> bind_rows()) |> bind_rows() |>
  mutate(current = detection_threshold == cur)
sig_tbl <- imap(gsea_runs, \(gl, nm) imap(gl, \(g, f) {
  x <- if (nrow(g)) g |> filter(pathway %in% sig_sets) else tibble()
  tibble(comparison = runs_ok[[nm]]$label, detection_threshold = as.numeric(f), gene_set = sig_sets,
         NES = x$NES[match(sig_sets, x$pathway)], p = x$p[match(sig_sets, x$pathway)], FDR = x$FDR[match(sig_sets, x$pathway)],
         measured_proteins = x$measured_proteins_in_set[match(sig_sets, x$pathway)])
}) |> bind_rows()) |> bind_rows()
# gene sets with enough measured proteins at other minimum sizes (current filter)
set_sizes <- imap(runs_ok, \(x, nm) {
  genes <- unique(unlist(str_split(x$res[[as.character(cur)]]$Assay, "_"))); sz <- map_int(sets, \(g) sum(unique(g) %in% genes))
  tibble(comparison = x$label, min_size = sort(unique(c(5L, 10L, 15L, 20L, as.integer(gmin)))), max_size = gmax) |>
    mutate(gene_sets_eligible = map_int(min_size, \(m) sum(sz >= m & sz <= gmax)), current = min_size == gmin)
}) |> bind_rows()

# =============================================================================================================================
# 4  Group sizes per visit against the minimums
# =============================================================================================================================
visits <- sort(unique(ad_isf$visit_num))
grp_tbl <- map(visits, \(v) {
  L <- ad_isf |> filter(site == "L", visit_num == v) |> distinct(SubjectID, relapse, relapse_visit)
  at_rel <- L$relapse %in% "relapse" & coalesce(L$relapse_visit <= v, FALSE)
  nr <- sum(L$relapse %in% "relapse" & !at_rel); nn <- sum(L$relapse %in% "non-relapse")
  n_ad <- n_distinct(ad_isf$SubjectID[ad_isf$visit_num == v])
  lab <- if (nr && nn) choose(nr + nn, nr) else NA_real_
  tibble(visit = paste0("V", v), AD_patients_dISF = n_ad, pairs_lesion_site_non_lesional = n_distinct(pairs_at(v)$SubjectID),
         relapse = nr, non_relapse = nn, relapsers_at_relapse_visit_left_out = sum(at_rel),
         step13_visit_analysed = n_ad >= vmin, relapse_comparison_tested = min(nr, nn) >= min_n,
         labellings = lab, smallest_exact_p = if (is.na(lab)) NA_real_ else (if (nr == nn) 2 else 1) / lab)
}) |> bind_rows() |> mutate(exact_p_below_0.05_possible = smallest_exact_p < 0.05)

# =============================================================================================================================
# 5  Relapse target criteria at V1 (as in the V1 relapse screen): where the candidates are lost
# =============================================================================================================================
funnel <- if (!"relapse_V1" %in% names(runs_ok)) tibble() else {
  x <- runs_ok$relapse_V1; r <- x$res[[as.character(cur)]]
  M <- patient_matrix(clean, x$s, r$OlinkID); M <- M[r$OlinkID, sort(colnames(M)), drop = FALSE]
  rc <- relabel_check(M, x$s$grp[match(colnames(M), x$s$SubjectID)] == "relapse")
  r <- r |> left_join(rc$per_protein |> select(OlinkID, p_exact_up, p_exact_down), by = "OlinkID")
  relapse_funnel(r, cur)
}

# =============================================================================================================================
# 6  Inventory of the thresholds in the code, and the verdict per threshold
# =============================================================================================================================
inventory <- threshold_inventory(cfg)
verdict <- threshold_verdicts(cfg, samples_qc, det_tbl, det, removed_tbl, cmp_tbl, grp_tbl, set_sizes, funnel, cur, fdr, f0, poi_by_group, poi_missing, det_groups)

# =============================================================================================================================
# Excel file
# =============================================================================================================================
writexl::write_xlsx(list(
  README = tibble(item = c("question", "verdict", "samples", "detection", "comparisons", "removed_by_filter", "gene_sets", "group_sizes", "relapse_funnel", "inventory", "data"),
    note = c("Which thresholds remove data or decide results, and are they too strict for these data?",
             "one row per threshold: current value, what it does, its effect on these data, verdict",
             sprintf("Olink SampleQC FAIL samples are removed; WARN samples and outliers (median or IQR of NPX beyond mean +/- %g SD) are only flagged; low volume (< %s uL) is only flagged. Sheet outliers_by_sd: outliers at other SD values", k_cur, paste(sprintf("%s %g", names(lv), lv), collapse = ", ")),
             sprintf("assays kept at each detection threshold (step 02 rule: share of samples above LOD in at least one detection group - sheet detection_groups; or share of all samples of the matrix). Current: %g", cur),
             sprintf("five comparisons re-run at each detection threshold (limma; a protein is tested if above LOD in >= that share of the samples of one of the two groups): proteins tested, p < 0.05, FDR levels; relapse V1 also the median of p < 0.05 proteins under all random labellings"),
             sprintf("proteins tested without the filter (threshold %g) but removed by the current one: how many have p < 0.05, against 5%% expected by chance; sheet removed_p_below_0.05 lists them with the share above LOD per group", f0),
             sprintf("GSEA at detection thresholds %s: gene sets with FDR below 0.05 / 0.1 / 0.25, the signature gene set(s) (%s); gene sets eligible at other minimum sizes (measured proteins %d-%d now)",
                     paste(fr_gsea, collapse = ", "), paste(sig_sets, collapse = ", "), gmin, gmax),
             sprintf("patients per visit: AD with dISF (step 13 analyses a visit with >= %d), relapse vs non-relapse at the lesion site (tested with >= %d per group; relapsers at their relapse visit left out), smallest possible exact p", vmin, min_n),
             sprintf("V1 lesion site, relapse vs non-relapse: proteins left after each criterion of the V1 relapse screen (p < 0.05, exact one-sided p < 0.05, above LOD in >= %.0f%% of the group), and with looser criteria", 100 * cur),
             "every threshold in the code: value, config key or hard-coded, where, what it does, type",
             manifest_txt)),
  verdict = verdict, inventory = inventory,
  samples = samples_qc, outliers_by_sd = outliers_by_sd, removed_samples = removed_samples,
  detection = det_tbl, detection_groups = det_groups,
  proteins_of_interest = poi_by_group |> select(list, gene, Assay, OlinkID, matrix, det_group, pct_above_LOD, kept_by_step02 = keep) |> sig3(),
  comparisons = cmp_tbl, removed_by_filter = removed_tbl, removed_p_below_0.05 = removed_hits |> sig3(),
  gene_sets = gsea_tbl, signature_gene_set = sig_tbl |> sig3(), gene_set_sizes = set_sizes,
  group_sizes = grp_tbl |> sig3(), relapse_funnel = funnel),
  out_path(cfg, out, "thresholds.xlsx"))

# =============================================================================================================================
# Figures: PNG (300 dpi, 16:9, large fonts - safe to paste into PowerPoint) and PDF, collected for the PDF summary
# =============================================================================================================================
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
figs <- list()
fig <- function(p, file, title, note = NULL, width = 12.5, height = 6.4) {
  if (is.null(p)) return(invisible(NULL))
  if (!is.null(note)) note <- str_wrap(note, 115)
  pt <- p + labs(title = title, subtitle = note)
  for (type in c("pdf", "png")) {
    open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
    print(set_family(pt, fi$family))
    close_device(fi)
  }
  figs[[length(figs) + 1]] <<- list(p = p, title = title, note = note)
  invisible(p)
}
dc <- det_tbl |> filter(current)
fig(plot_detection_curve(det_tbl, cur, colours), "1_assays_kept_by_detection_threshold", "Proteins kept at each detection (LOD) threshold",
    sprintf("current setting %g (dashed): %s. Step 02 rule: share above LOD in at least one detection group (solid); dashed: of all samples of the matrix", cur,
            paste(sprintf("%s %d of %d (%.0f%%)", dc$matrix, dc$kept_in_one_group, dc$assays, dc$pct_kept_in_one_group), collapse = ", ")))
if (nrow(cmp_tbl) && "proteins_tested" %in% names(cmp_tbl))
  fig(plot_comparisons(cmp_tbl |> filter(!is.na(proteins_tested)), cur, fdr, colours), "2_comparisons_by_detection_threshold",
      "Results of five comparisons at each detection threshold",
      sprintf("proteins tested, with p < 0.05 and with FDR < %.2g; relapse V1: grey = median under random labels (chance); dashed line = current setting %g", fdr, cur))
if (nrow(removed_tbl))
  fig(plot_removed(removed_tbl, colours), "3_removed_proteins_vs_chance", "Do the proteins removed by the detection filter carry signal?",
      sprintf("proteins with p < 0.05 among those the current filter (%g) removes and among those it keeps, against the 5%% expected by chance (no filter = %g)", cur, f0))
if (nrow(poi_by_group))
  fig(plot_poi(poi_by_group, cur, colours), "4_proteins_of_interest_above_LOD", "Proteins of interest: share of samples above LOD per group",
      sprintf("x = proteins not kept by step 02 (below %.0f%% in every group)%s", 100 * cur,
              if (length(poi_missing)) sprintf("; not on the panel: %s", paste(poi_missing, collapse = ", ")) else ""), height = 6.4)
fig(plot_group_sizes(grp_tbl, min_n, colours), "5_group_sizes_per_visit", "Relapse vs non-relapse at the lesion site: patients per visit",
    sprintf("relapsers at their relapse visit left out; dashed line = stats$min_group_n (%d); label = smallest possible exact p", min_n))
if (nrow(funnel)) fig(plot_funnel(funnel, colours), "6_relapse_criteria_V1", "V1 relapse screen: proteins left after each criterion",
                      "dISF lesion site V1, relapse vs non-relapse; the target list needs all criteria (top to bottom)")

# =============================================================================================================================
# PDF summary
# =============================================================================================================================
source("R/report.R")
pdf_open(out_path(cfg, out, "Thresholds_summary.pdf"))
section("Verdict", page_text("Thresholds: how strict are they for these data?", verdict_lines(verdict),
                             subtitle = sprintf("O-MicroAD, Olink Explore HT. Generated %s. %s. Details: thresholds.xlsx (sheet verdict)", format(Sys.Date()), manifest_txt),
                             size = 9.5, width = 150))
section("Comparisons", page_table("Five comparisons at each detection threshold", cmp_tbl |> select(any_of(c("comparison", "samples", "detection_threshold", "proteins_tested",
  "p_below_0.05", sprintf("FDR_below_%g", fdr_lv), "chance_median_p_below_0.05", "note"))) |> mutate(comparison = str_trunc(comparison, 45)),
  note = "detection_threshold: share of the samples of one of the two groups above LOD needed for a protein to be tested", rows_per_page = 32))
section("Removed", page_table("Proteins removed by the detection filter vs chance", removed_tbl |> transmute(comparison = str_trunc(comparison, 45), removed = removed_by_filter,
  `removed p<0.05` = removed_p_below_0.05, `expected` = round(removed_expected_by_chance, 1), kept, `kept p<0.05` = kept_p_below_0.05, `expected ` = round(kept_expected_by_chance, 1),
  `FDR hits now` = FDR_hits_current, `FDR hits no filter` = FDR_hits_without_filter),
  note = sprintf("no filter = detection threshold %g; expected = 5%% of the proteins (chance)", f0)))
section("Group sizes", page_table("Patients per visit and group", grp_tbl |> transmute(visit, AD = AD_patients_dISF, pairs = pairs_lesion_site_non_lesional, relapse, non_relapse,
  left_out = relapsers_at_relapse_visit_left_out, step13 = if_else(step13_visit_analysed, "yes", "no"), relapse_test = if_else(relapse_comparison_tested, "yes", "no"),
  min_exact_p = fmt_p(smallest_exact_p)),
  note = sprintf("step13: visit analysed (>= %d AD patients); relapse_test: both groups >= %d; min_exact_p: smallest exact p possible with these group sizes", vmin, min_n)))
for (f in figs) section(f$title, page_plot(f$p + labs(title = f$title, subtitle = f$note)))
section("Inventory", page_table("Every threshold in the code", inventory |> transmute(area, threshold = str_trunc(threshold, 40), value = str_trunc(value, 18),
  where = str_trunc(where, 34), type), rows_per_page = 30))
invisible(grDevices::dev.off())
msg("Thresholds: %s", file.path(cfg$paths$output, out))
