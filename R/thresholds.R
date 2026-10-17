# Helpers for the threshold analysis (thresholds/): the relapse criteria funnel, the verdict per threshold, the inventory of
# all thresholds in the code, and the figures (16:9, large fonts, for PowerPoint).

#' Proteins left after each criterion of the V1 relapse screen (relapse_v1/: p < 0.05, exact one-sided p < 0.05, above LOD in
#' >= min_f of the group with the higher values), and with looser criteria. r: de_limma() result with p_exact_up / p_exact_down.
relapse_funnel <- function(r, min_f) {
  lod_rel <- r$`pct_above_LOD_relapse` >= 100 * min_f; lod_non <- r$`pct_above_LOD_non-relapse` >= 100 * min_f
  lod_rel25 <- r$`pct_above_LOD_relapse` >= 25; lod_non25 <- r$`pct_above_LOD_non-relapse` >= 25
  up <- coalesce(r$logFC > 0, FALSE); dn <- coalesce(r$logFC < 0, FALSE)
  p05 <- coalesce(r$p < 0.05, FALSE); p10 <- coalesce(r$p < 0.1, FALSE)
  ex_up <- coalesce(r$p_exact_up < 0.05, FALSE); ex_dn <- coalesce(r$p_exact_down < 0.05, FALSE)
  lod_txt <- sprintf("+ above LOD in >= %.0f%% of the group", 100 * min_f)
  bind_rows(
    tibble(direction = "higher in relapsers", criterion = c("proteins tested", "p < 0.05", "+ exact p < 0.05", lod_txt),
           kind = "screen", proteins = c(nrow(r), sum(up & p05), sum(up & p05 & ex_up), sum(up & p05 & ex_up & coalesce(lod_rel, FALSE)))),
    tibble(direction = "lower in relapsers", criterion = c("proteins tested", "p < 0.05", "+ exact p < 0.05", lod_txt),
           kind = "screen", proteins = c(nrow(r), sum(dn & p05), sum(dn & p05 & ex_dn), sum(dn & p05 & ex_dn & coalesce(lod_non, FALSE)))),
    tibble(direction = "higher in relapsers", kind = "looser",
           criterion = c("looser: without the exact p", "looser: above LOD in >= 25% of the group", "looser: p < 0.1 (with exact p and LOD)"),
           proteins = c(sum(up & p05 & coalesce(lod_rel, FALSE)), sum(up & p05 & ex_up & coalesce(lod_rel25, FALSE)), sum(up & p10 & ex_up & coalesce(lod_rel, FALSE)))),
    tibble(direction = "lower in relapsers", kind = "looser",
           criterion = c("looser: without the exact p", "looser: above LOD in >= 25% of the group", "looser: p < 0.1 (with exact p and LOD)"),
           proteins = c(sum(dn & p05 & coalesce(lod_non, FALSE)), sum(dn & p05 & ex_dn & coalesce(lod_non25, FALSE)), sum(dn & p10 & ex_dn & coalesce(lod_non, FALSE)))))
}

#' One row per threshold: the setting, its effect on these data and a verdict (is it too strict?).
threshold_verdicts <- function(cfg, samples_qc, det_tbl, det, removed_tbl, cmp_tbl, grp_tbl, set_sizes, funnel, cur, fdr, f0, poi_by_group, poi_missing,
                               det_groups, removed_samples, id_mismatch, na_tbl) {
  num <- \(x) if (length(x) && !is.na(x[1])) x[1] else 0
  sq_of <- \(m, col) num(samples_qc[[col]][samples_qc$matrix == m])
  rows <- list()
  add <- function(threshold, setting, effect, verdict) rows[[length(rows) + 1]] <<- tibble(threshold, setting, effect, verdict)
  # samples
  part_fail <- if ("pct_values_FAIL" %in% names(removed_samples) && nrow(removed_samples)) median(removed_samples$pct_values_FAIL, na.rm = TRUE) else NA_real_
  add("Olink sample QC: FAIL", "a sample with a FAIL in any block is removed as a whole (always)",
      sprintf("%d of %d dISF and %d of %d serum samples removed%s", sq_of("dISF", "removed_FAIL"), sq_of("dISF", "samples"), sq_of("serum", "removed_FAIL"), sq_of("serum", "samples"),
              if (is.na(part_fail)) "" else sprintf("; median %.0f%% of their values failed (sheet removed_samples)", part_fail)),
      if (!is.na(part_fail) && part_fail < 50)
        "stricter than needed: most values of the removed samples passed QC - removing only the failed blocks would keep the rest (not changed here)"
      else if (max(sq_of("dISF", "removed_FAIL") / max(sq_of("dISF", "samples"), 1), sq_of("serum", "removed_FAIL") / max(sq_of("serum", "samples"), 1)) > 0.1)
        "many samples removed - check qc/sample_qc.csv" else "standard Olink QC - not limiting")
  n_npx_only <- sum(id_mismatch$problem %in% "in NPX data, not in manifest"); n_man_only <- sum(id_mismatch$problem %in% "in manifest, not in NPX data")
  add("Samples not in the manifest", "NPX samples without a row in the manifest are left out (step 02)",
      sprintf("%d samples of the NPX files are not in the manifest; %d manifest samples have no NPX data", n_npx_only, n_man_only),
      if (n_npx_only) "check qc/sample_id_mismatches.csv: these samples are in no analysis - add them to the manifest if they belong to the study" else "none - not limiting")
  add("Olink sample QC: WARN", sprintf("qc$drop_sample_qc_warn = %s", tolower(isTRUE(cfg$qc$drop_sample_qc_warn))),
      sprintf("%d dISF and %d serum WARN samples %s", sq_of("dISF", "WARN_kept"), sq_of("serum", "WARN_kept"), if (isTRUE(cfg$qc$drop_sample_qc_warn)) "removed" else "kept (flagged)"),
      if (isTRUE(cfg$qc$drop_sample_qc_warn)) "removes samples - set false to keep them" else "not limiting (kept)")
  add("Sample outliers", sprintf("qc$outlier_sd = %g, qc$drop_outliers = %s", cfg$qc$outlier_sd %||% 3, tolower(isTRUE(cfg$qc$drop_outliers))),
      sprintf("%d dISF and %d serum samples flagged", sq_of("dISF", "outliers_flagged"), sq_of("serum", "outliers_flagged")),
      if (isTRUE(cfg$qc$drop_outliers)) "removes the flagged samples" else "not limiting (flagged only; sheet outliers_by_sd)")
  add("Low sample volume", sprintf("qc$low_volume_ul = %s", paste(sprintf("%s %g", names(unlist(cfg$qc$low_volume_ul)), unlist(cfg$qc$low_volume_ul)), collapse = ", ")),
      if (sq_of("dISF", "with_volume") + sq_of("serum", "with_volume") == 0) "no sample volumes in the manifest" else
        sprintf("%d dISF and %d serum samples flagged", sq_of("dISF", "low_volume_flagged"), sq_of("serum", "low_volume_flagged")),
      "not limiting (flagged only)")
  # proteins
  dc <- det_tbl |> filter(current)
  kept_txt <- paste(sprintf("%s keeps %d of %d assays (%.0f%%)", dc$matrix, dc$kept_in_one_group, dc$assays, dc$pct_kept_in_one_group), collapse = "; ")
  strong <- removed_tbl |> filter(removed_p_below_0.05 >= 10, removed_ratio_to_chance > 1.5)
  rm_txt <- if (nrow(removed_tbl)) paste(sprintf("%s: %d of %d removed proteins with p < 0.05 (chance %.0f)", removed_tbl$short,
                                                 removed_tbl$removed_p_below_0.05, removed_tbl$removed_by_filter, removed_tbl$removed_expected_by_chance), collapse = "; ") else ""
  add("Detection (LOD) filter", sprintf("qc$min_detect_frac = %g (above LOD in >= %.0f%% of the samples of one group)", cur, 100 * cur),
      paste0(kept_txt, ". ", rm_txt),
      if (nrow(strong)) sprintf("possibly too strict for %s: the removed proteins carry more p < 0.05 than chance - see sheet removed_p_below_0.05 (values mostly below LOD: read with care; lower qc$min_detect_frac to include them)",
                                paste(strong$short, collapse = "; ")) else
        "not too strict: the removed proteins carry about as many p < 0.05 as chance (noise below LOD); keeping them only adds tests")
  small <- det_groups |> filter(samples < 10) |> arrange(samples)
  add("Detection groups (step 02 filter)", "a protein is kept if above LOD in enough samples of ANY one detection group",
      if (nrow(small)) paste(sprintf("%s %s: %d samples (%d above LOD suffice)", small$matrix, str_remove(small$det_group, ":$"), small$samples, small$samples_needed_above_LOD), collapse = "; ")
      else "all detection groups have >= 10 samples",
      if (nrow(small)) "lenient, not strict: small groups let a protein pass with very few samples above LOD; the dedicated analyses use their own filter on the two groups compared"
      else "not limiting")
  lost <- poi_by_group |> distinct(matrix, gene, keep) |> filter(!coalesce(keep, TRUE)) |> group_by(matrix) |>
    summarise(genes = paste(sort(gene), collapse = ", "), n = n(), .groups = "drop")
  add("Detection filter: proteins of interest", sprintf("the protein lists in config.yml (%d proteins)", n_distinct(poi_by_group$gene) + length(poi_missing)),
      paste0(if (nrow(lost)) paste(sprintf("not kept in %s (below %.0f%% above LOD in every group): %s", lost$matrix, 100 * cur, lost$genes), collapse = "; ")
             else "all kept in both matrices", if (length(poi_missing)) sprintf("; not on the panel: %s", paste(poi_missing, collapse = ", ")) else ""),
      if (nrow(lost)) sprintf("mostly below LOD in that matrix, so tests on them are noise-dominated; to test them anyway lower qc$min_detect_frac (share above LOD per group: sheet proteins_of_interest; values: export/)")
      else "not limiting")
  no_lod <- det |> group_by(matrix) |> summarise(n = sum(!lod_available), .groups = "drop")
  add("Assays without LOD", "kept by step 02 (cannot be judged); not tested by the dedicated analyses",
      paste(sprintf("%s %d", c(ISF = "dISF", Serum = "serum")[no_lod$matrix], no_lod$n), collapse = ", "),
      if (sum(no_lod$n)) "check qc/lod.csv: these assays are missing from the dedicated analyses" else "none - not limiting")
  add("Missing values in the main models", "a protein with > 20% missing values in a model's samples is not tested; the rest are imputed (hard-coded)",
      if (nrow(na_tbl)) paste(sprintf("%s: %d of %d kept proteins have > 20%% missing values over all samples", na_tbl$matrix, na_tbl$proteins_over_20pct_missing, na_tbl$proteins), collapse = "; ")
      else "npx_wide.rds not found",
      if (nrow(na_tbl) && sum(na_tbl$proteins_over_20pct_missing)) "removes these proteins from the main models (steps 04-17)" else "not limiting (hardly any missing values)")
  # significance
  fcol <- sprintf("FDR_below_%g", fdr)
  cc <- if (fcol %in% names(cmp_tbl)) cmp_tbl |> filter(current, !is.na(.data[[fcol]])) else tibble()
  rel <- cc |> filter(short == "relapse V1"); rel_n <- if (nrow(rel)) rel$samples[1] else "the relapse"
  add("False discovery rate", sprintf("stats$fdr = %g", fdr),
      if (nrow(cc)) paste(sprintf("%s: %d with FDR < %g (p < 0.05: %d%s)", cc$short, cc[[fcol]], fdr, cc$p_below_0.05,
                                  if ("chance_median_p_below_0.05" %in% names(cc)) ifelse(is.na(cc$chance_median_p_below_0.05), "", sprintf(", chance %.0f", cc$chance_median_p_below_0.05)) else ""),
                          collapse = "; ") else "no comparison run",
      if (nrow(rel) && rel[[fcol]][1] == 0) sprintf("too strict for relapse vs non-relapse (%s gives no FDR hits): the relapse screens use p < 0.05 with the exact relabelling check instead - compare p < 0.05 with chance", rel_n)
      else "adequate where groups are larger; relapse screens use p < 0.05 with the exact relabelling check")
  # group sizes
  nt <- grp_tbl |> filter(!relapse_comparison_tested)
  add("Minimum group size", sprintf("stats$min_group_n = %d", cfg$stats$min_group_n %||% 3),
      if (nrow(nt)) sprintf("relapse vs non-relapse not tested at %s", paste(sprintf("%s (%d vs %d)", nt$visit, nt$relapse, nt$non_relapse), collapse = ", ")) else "all visits tested",
      if (any(nt$visit %in% paste0("V", 1:4))) "limits V1-V4 - consider a smaller minimum (results then rest on very few patients)" else "not limiting for V1-V4")
  ns <- grp_tbl |> filter(!step13_visit_analysed)
  add("Visits in the visit course (step 13)", sprintf("stats$visit_min_subjects = %d", cfg$stats$visit_min_subjects %||% 5),
      if (nrow(ns)) sprintf("not analysed: %s", paste(sprintf("%s (%d AD patients)", ns$visit, ns$AD_patients_dISF), collapse = ", ")) else "all visits analysed",
      if (nrow(ns)) "few patients left at these visits - lowering the minimum would rest on very few patients" else "not limiting")
  ex <- grp_tbl |> filter(relapse_comparison_tested)
  add("Exact relabelling p < 0.05 (relapse)", "p < 0.05 over all ways to choose the relapse group",
      paste(sprintf("%s smallest possible p %s", ex$visit, formatC(ex$smallest_exact_p, digits = 2, format = "g")), collapse = "; "),
      if (any(!ex$exact_p_below_0.05_possible)) sprintf("cannot be reached at %s", paste(ex$visit[!ex$exact_p_below_0.05_possible], collapse = ", ")) else
        "reachable at every tested visit, but only by near-perfect separation where the groups are smallest")
  # gene sets
  ss <- set_sizes |> group_by(min_size) |> summarise(n = min(gene_sets_eligible), .groups = "drop")
  gm <- cfg$enrichment$min_size %||% 10
  add("Gene set size", sprintf("enrichment$min_size = %d, max_size = %d (measured proteins)", gm, cfg$enrichment$max_size %||% 500),
      paste(sprintf(">= %d proteins: %d sets", ss$min_size, ss$n), collapse = "; "),
      "standard - smaller sets give unstable enrichment scores")
  # relapse targets
  if (nrow(funnel)) {
    fu <- funnel |> filter(kind == "screen", direction == "higher in relapsers")
    add("Relapse target criteria (V1)", "p < 0.05 + exact p < 0.05 + above LOD in the relapsers",
        paste(sprintf("%s %d", fu$criterion, fu$proteins), collapse = " -> "),
        sprintf("intended: with %s the exact p protects against chance; sheet relapse_funnel shows looser variants", rel_n))
  }
  add("Firm verdict in the key questions (step 15)", ">= 10 patients per group (hard-coded)",
      sprintf("relapse groups here: %s", rel_n),
      "by design never reached with these groups: relapse verdicts stay 'possible' at most")
  bind_rows(rows) |>
    mutate(status = case_when(str_detect(verdict, "^by design") ~ "by design",
                              str_detect(verdict, "^(possibly too strict|too strict|stricter than needed|limits|cannot be reached|many samples removed|removes|check)") ~ "may limit these data",
                              TRUE ~ "not limiting"), .after = threshold)
}

#' The verdict as text for the first PDF page: what may limit the analysis first, then the rest.
verdict_lines <- function(verdict) {
  one <- \(v) sprintf("%s [%s]: %s -> %s", v$threshold, v$setting, v$effect, v$verdict)
  grp <- c("may limit these data", "by design", "not limiting")
  n <- table(factor(verdict$status, grp))
  c(sprintf("%d thresholds checked: %d may limit these data, %d limit by design, %d not limiting.", nrow(verdict), n[[1]], n[[2]], n[[3]]),
    unlist(map(grp, \(g) { v <- verdict[verdict$status == g, ]; if (nrow(v)) c(paste("##", str_to_sentence(g)), map_chr(seq_len(nrow(v)), \(i) one(v[i, ]))) })))
}

# ---- figures --------------------------------------------------------------------------------------------------------------
slide_theme <- \(base = 16) theme_bw(base_size = base) + theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title.position = "plot")
pct_axis <- \(x) paste0(round(100 * x), "%")

plot_detection_curve <- function(det_tbl, cur, colours) {
  d <- det_tbl |> select(matrix, detection_threshold, `in at least one group (step 02)` = kept_in_one_group, `of all samples of the matrix` = kept_of_all_samples) |>
    pivot_longer(-c(matrix, detection_threshold), names_to = "rule", values_to = "assays")
  ggplot(d, aes(detection_threshold, assays, colour = matrix, linetype = rule)) +
    geom_vline(xintercept = cur, linetype = 2, colour = "grey55") +
    geom_line(linewidth = 1) + geom_point(size = 2.6) +
    scale_colour_manual(values = c(dISF = colours$above, serum = colours$below), name = NULL) +
    scale_x_continuous(labels = pct_axis, breaks = sort(unique(d$detection_threshold))) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.12))) +
    labs(x = "share of samples above LOD needed", y = "assays kept", linetype = NULL) + slide_theme() +
    guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2))
}

plot_comparisons <- function(cmp_tbl, cur, fdr, colours) {
  fcol <- sprintf("FDR_below_%g", fdr); flab <- sprintf("FDR < %g", fdr); lv <- c("proteins tested", "p < 0.05", flab)
  lvc <- unique(str_wrap(cmp_tbl$comparison, 24))                     # in the order of the analysis
  d <- cmp_tbl |> transmute(comparison = factor(str_wrap(comparison, 24), lvc), detection_threshold, `proteins tested` = proteins_tested, `p < 0.05` = p_below_0.05,
                            !!flab := .data[[fcol]]) |>
    pivot_longer(-c(comparison, detection_threshold), names_to = "measure", values_to = "n") |> mutate(measure = factor(measure, lv))
  ch <- if ("chance_median_p_below_0.05" %in% names(cmp_tbl)) cmp_tbl |> filter(!is.na(chance_median_p_below_0.05)) |>
    transmute(comparison = factor(str_wrap(comparison, 24), lvc), detection_threshold, measure = factor("p < 0.05", lv), n = chance_median_p_below_0.05) else tibble()
  p <- ggplot(d, aes(detection_threshold, n)) + geom_vline(xintercept = cur, linetype = 2, colour = "grey60") +
    geom_line(colour = colours$above, linewidth = 0.9) + geom_point(colour = colours$above, size = 1.8)
  if (nrow(ch)) p <- p + geom_line(data = ch, colour = "grey55", linewidth = 0.9) + geom_point(data = ch, colour = "grey55", size = 1.8)
  p + facet_grid(measure ~ comparison, scales = "free_y", switch = "y") +
    scale_x_continuous(labels = pct_axis, breaks = c(0, 0.5, 1), limits = c(0, 1)) + scale_y_continuous(limits = c(0, NA)) +
    labs(x = "share of samples above LOD needed (in one of the two groups)", y = NULL) + slide_theme(12) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))
}

plot_removed <- function(removed_tbl, colours) {
  d <- removed_tbl |> transmute(comparison = factor(str_wrap(comparison, 22), unique(str_wrap(comparison, 22))),
                                `removed by the filter` = removed_ratio_to_chance, `kept by the filter` = kept_ratio_to_chance,
                                lab_r = sprintf("%d of %d", removed_p_below_0.05, removed_by_filter), lab_k = sprintf("%d of %d", kept_p_below_0.05, kept)) |>
    pivot_longer(c(`removed by the filter`, `kept by the filter`), names_to = "set", values_to = "ratio") |>
    mutate(label = if_else(set == "removed by the filter", lab_r, lab_k), set = factor(set, c("kept by the filter", "removed by the filter")))
  ggplot(d, aes(comparison, ratio, fill = set)) + geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(aes(label = label), position = position_dodge(width = 0.8), vjust = -0.4, size = 4.2) +
    scale_fill_manual(values = c(`kept by the filter` = colours$above, `removed by the filter` = colours$below), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "p < 0.05: observed / expected by chance\n(1 = chance, dashed)") + slide_theme(15)
}

plot_poi <- function(poi_by_group, cur, colours, max_rows = 40) {
  best <- poi_by_group |> group_by(gene) |> summarise(best = max(pct_above_LOD, na.rm = TRUE), .groups = "drop")
  at_risk <- best |> filter(best < 100 * min(cur + 0.25, 1)) |> arrange(best) |> slice_head(n = max_rows)   # near or below the threshold
  if (!nrow(at_risk)) return(NULL)
  d <- poi_by_group |> filter(gene %in% at_risk$gene) |>
    mutate(col = paste(matrix, str_remove(det_group, ":$"), sep = "\n"), gene = factor(gene, rev(at_risk$gene)),
           label = paste0(round(pct_above_LOD), if_else(coalesce(keep, TRUE), "", "x")))
  ggplot(d, aes(col, gene, fill = pct_above_LOD)) + geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label), size = 3.2) +
    scale_fill_gradient(low = "white", high = colours$above, limits = c(0, 100), name = "% above\nLOD") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 11) + theme(panel.grid = element_blank(), plot.title.position = "plot")
}

plot_group_sizes <- function(grp_tbl, min_n, colours) {
  d <- grp_tbl |> select(visit, relapse, `non-relapse` = non_relapse, smallest_exact_p) |>
    pivot_longer(c(relapse, `non-relapse`), names_to = "group", values_to = "patients") |> mutate(group = factor(group, c("relapse", "non-relapse")))
  lab <- grp_tbl |> mutate(y = pmax(relapse, non_relapse), label = if_else(is.na(smallest_exact_p), "not testable", paste0("min p ", formatC(smallest_exact_p, digits = 2, format = "g"))))
  ggplot(d, aes(visit, patients, fill = group)) + geom_hline(yintercept = min_n, linetype = 2, colour = "grey40") +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_text(data = lab, aes(visit, y, label = label), inherit.aes = FALSE, vjust = -0.6, size = 4.5) +
    scale_fill_manual(values = c(relapse = colours$above, `non-relapse` = colours$below), name = NULL) +
    scale_y_continuous(breaks = function(l) seq(0, ceiling(l[2])), expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "patients (lesion site, dISF)") + slide_theme()
}

plot_funnel <- function(funnel, colours) {
  d <- funnel |> filter(kind == "screen") |> mutate(criterion = factor(criterion, rev(unique(criterion))))
  ggplot(d, aes(proteins, criterion, fill = direction)) + geom_col(width = 0.7) +
    geom_text(aes(label = proteins), hjust = -0.2, size = 5) + facet_wrap(~direction) +
    scale_fill_manual(values = c(`higher in relapsers` = colours$above, `lower in relapsers` = colours$below), guide = "none") +
    scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300, 1000, 3000), expand = expansion(mult = c(0, 0.15))) +
    labs(x = "proteins (log scale)", y = NULL) + slide_theme()
}

#' Every threshold in the code that removes data or decides a result: value (from config.yml where it is set there), where it
#' is applied, what it does, and whether it removes data, decides significance / a verdict, or only flags.
threshold_inventory <- function(cfg) {
  val <- function(key) {                                   # current value of a config key such as "qc$min_detect_frac"
    x <- cfg; for (k in str_split(key, "\\$")[[1]]) x <- x[[k]]
    if (is.null(x)) return("(default)")
    v <- unlist(x); nm <- names(v); v <- if (is.logical(v)) tolower(v) else as.character(v)
    paste(if (!is.null(nm)) paste(nm, v) else v, collapse = ", ")
  }
  inv <- tribble(
    ~area, ~threshold, ~config_key, ~fixed_value, ~where, ~effect, ~type,
    # samples
    "samples", "Olink SampleQC FAIL", "", "FAIL in any block", "R/qc.R sample_qc(), scripts/02_import_qc.R", "a sample is FAIL if any of its rows (blocks) fails; the whole sample is then removed from all analyses", "removes data",
    "samples", "samples missing from the manifest", "", "not in the manifest", "scripts/02_import_qc.R", "NPX samples without a manifest row are dropped (listed in qc/sample_id_mismatches.csv)", "removes data",
    "samples", "relapse label spelling", "", "relapse / non-relapse / no relapse / dropout ...", "R/metadata.R harmonise_relapse()", "any other spelling of Relapse becomes missing: the patient leaves every relapse comparison", "removes data",
    "samples", "one serum sample per person (AD vs healthy)", "", "V1", "R/design.R", "MicroAD serum AD vs healthy and the cross-sectional serum models use V1 only", "removes data",
    "samples", "before-relapse samples", "", "visit < relapse visit (steps 14, 17)", "steps 14, 17; visits_v1_v4/; steps 04, 05, 11", "steps 14 and 17 use visits before the relapse; visits_v1_v4 drops relapsers from their relapse visit on; the relapse models of steps 04, 05, 11 have no such filter", "removes data",
    "samples", "Olink SampleQC WARN", "qc$drop_sample_qc_warn", "", "scripts/02_import_qc.R", "true = WARN samples removed; false = kept and flagged", "removes data / flag",
    "samples", "sample outlier rule", "qc$outlier_sd", "", "R/qc.R sample_qc()", "median or IQR of NPX beyond mean +/- k SD of the matrix (non-FAIL samples): flagged", "flag only",
    "samples", "remove outliers", "qc$drop_outliers", "", "scripts/02_import_qc.R", "true = flagged outliers removed", "removes data / flag",
    "samples", "low sample volume", "qc$low_volume_ul", "", "R/metadata.R flag_metadata()", "samples below the volume are flagged in metadata/data_flags.csv", "flag only",
    "samples", "dropouts and labels", "", "relapse label relapse / non-relapse", "relapse_v1/, relapse_followup/, visits_v1_v4/, step 15", "patients without a relapse / non-relapse label (dropouts) are left out of relapse comparisons", "removes data",
    "samples", "relapse visit in per-visit relapse tests", "visits_v1_v4$include_relapse_visit", "", "visits_v1_v4/", "false = a relapser's sample at the relapse visit (lesional again) is left out", "removes data",
    # proteins
    "proteins", "LOD source", "qc$lod_method", "", "R/qc.R add_lod()", "auto = Olink fixed LOD file, else LOD from the negative controls; decides below_lod for every value", "decides LOD",
    "proteins", "AssayQC WARN values", "qc$drop_assay_qc_warn", "", "scripts/02_import_qc.R", "true = values with AssayQC WARN set to missing; false = kept", "removes data / flag",
    "proteins", "AssayQC FAIL values", "", "never removed", "scripts/02_import_qc.R", "values with AssayQC FAIL are kept (only listed in qc/assay_qc_flags.csv)", "flag only",
    "proteins", "missing values in the main models", "", "> 20% missing", "R/models.R, R/correlate.R", "a protein with more than 20% missing values in a model's samples is not tested; the remaining missing values are replaced by the protein median", "removes data",
    "proteins", "assays without any value", "", "0 values in a matrix", "scripts/02_import_qc.R", "assay/matrix combinations without any NPX value (Olink EXCLUDED) are removed", "removes data",
    "proteins", "detection filter (main pipeline)", "qc$min_detect_frac", "", "R/qc.R assay_detection(), step 02 keep", "assay kept in a matrix if above LOD in >= this share of the samples of at least one detection group; assays without LOD kept", "removes data",
    "proteins", "detection filter (per comparison)", "qc$min_detect_frac", "", "R/visit3.R de_limma()", "dedicated analyses: protein tested if above LOD in >= this share of the samples of one of the two compared groups; assays without LOD are not tested", "removes data",
    "proteins", "detection tiers (dISF profile)", "", "> 10%, >= 90% above LOD", "scripts/09_isf_profile.R", "proteins described as sporadic / robust in the dISF profile", "labels only",
    "proteins", "QC overview rule", "qc_overview$rule", "", "scripts/02b_qc_overview.R", "which share of samples counts as measurable in the QC overview (serum cohorts: qc_overview$serum_cohorts)", "decides counts",
    # significance
    "significance", "false discovery rate", "stats$fdr", "", "steps 03-18, dedicated analyses", "significant = Benjamini-Hochberg FDR below this", "decides significance",
    "significance", "nominal p for candidate lists", "", "p < 0.05", "relapse_v1/, relapse_followup/, visits_v1_v4/, steps 12-17", "candidate and screening lists (not corrected for many proteins)", "decides significance",
    "significance", "exact relabelling p", "", "exact p < 0.05", "relapse_v1/, relapse_followup/", "relapse candidates need an exact p (all ways to choose the relapse group) below 0.05", "decides significance",
    "significance", "relapse target: LOD in the group", "qc$min_detect_frac", "", "relapse_v1/ (candidate_up / _down)", "a relapse target must be above LOD in >= this share of the group with the higher values", "decides significance",
    "significance", "relapse signature over visits", "", "p < 0.05 at >= 2 visits", "visits_v1_v4/", "proteins in the relapse signature: same direction with p < 0.05 at two or more tested visits", "decides significance",
    "significance", "marker groups (relapse follow-up)", "", ">= half of detected markers (>= 2), exact p < 0.05", "relapse_followup/", "a marker group (T cells, activation) counts as lower / higher only if at least half of its detected markers differ", "decides verdict",
    "significance", "below-LOD caution", "", "< 50% above LOD", "relapse_followup/, visits_v1_v4/", "results of proteins mostly below LOD are marked for caution", "flag only",
    "significance", "firm verdict in the key questions", "", ">= 10 patients per group, p < 0.05, CI excluding 0.5", "scripts/15_key_questions.R", "a relapse prediction is called firm only with these; otherwise at most 'possible'", "decides verdict",
    "significance", "dISF vs serum enrichment", "stats$min_rel_log2", "", "scripts/10_matrix_comparison.R", "minimum relative dISF/serum difference (log2) to call a protein enriched in dISF", "decides significance",
    # group sizes
    "group sizes", "minimum group size", "stats$min_group_n", "", "R/models.R fit_contrasts(), R/correlate.R, visits_v1_v4/", "a contrast is skipped if a group has fewer samples (main models count samples, visits_v1_v4 counts patients; some modules use a fixed 3, step 17 a fixed 2)", "removes tests",
    "group sizes", "visits in the visit course", "stats$visit_min_subjects", "", "R/design.R visit_specs(), step 13", "a visit is analysed only with at least this many AD patients with dISF", "removes tests",
    "group sizes", "minimum pairs for a correlation", "tnfrsf9_correlation$min_n", "", "scripts/16_tnfrsf9_correlation.R", "correlations need at least this many sample pairs", "removes tests",
    "group sizes", "visit 3 groups", "", ">= 3 samples per group", "visit3/visit3_analysis.R", "the visit 3 analysis stops if a group has fewer than 3 samples", "removes tests",
    "group sizes", "minimum samples per model", "", "4 (models, xL-NL pairs), 6 (step 17), 8 (step 16)", "R/models.R, R/correlate.R, R/design.R, steps 04, 10, 16, 17", "a model or test with fewer samples is skipped", "removes tests",
    "group sizes", "minimum pairs for other correlations", "", "3, 5, 8 or 10 by step", "steps 06, 08, 10, 12, 15, 16, 17", "correlations need 5 pairs (steps 06, 15, 17), 8 (steps 10, 16), 10 (steps 08, 12) or 3 (step 08)", "removes tests",
    "group sizes", "visits in steps 14 and 17", "", "no minimum", "scripts/14_serum_vs_disf.R, scripts/17_disf_serum_signatures.R", "every visit with AD samples is analysed (stats$visit_min_subjects applies to step 13 only)", "removes tests",
    # gene sets
    "gene sets", "gene set size (GSEA)", "enrichment$min_size", "", "steps 09, 02b; dedicated GSEA", "gene sets with fewer measured proteins are not tested", "removes tests",
    "gene sets", "gene set size (step 07)", "", "5 (min of enrichment$min_size and 5)", "scripts/07_enrichment.R", "step 07 tests gene sets from 5 measured proteins; sets of 5-9 are marked small_set but stay in the FDR", "removes tests",
    "gene sets", "gene set size (relapse follow-up re-test)", "", "3", "relapse_followup/relapse_followup_analysis.R", "the immune gene sets re-tested after removing the sample level need 3 measured proteins", "removes tests",
    "gene sets", "largest gene set", "enrichment$max_size", "", "steps 07, 09, 17; dedicated GSEA", "gene sets with more measured proteins are not tested", "removes tests",
    "gene sets", "gene set size (signatures)", "signatures$min_size", "", "scripts/17_disf_serum_signatures.R", "smallest gene set in the over-representation analysis of step 17", "removes tests",
    "gene sets", "gene set size within a protein class", "qc_overview$pathway_min_size_class", "", "scripts/02b_qc_overview.R", "smallest gene set for the dISF-only pathways within a protein class", "removes tests")
  inv |> mutate(value = if_else(config_key == "", fixed_value, map_chr(config_key, \(k) if (k == "") "" else val(k))), config_key = if_else(config_key == "", "hard-coded", config_key)) |>
    select(area, threshold, value, config_key, where, effect, type)
}

#' Places where the same kind of threshold is defined differently in different steps (from a review of the code).
threshold_inconsistencies <- function() tribble(
  ~topic, ~finding,
  "assays without LOD", "kept by the step 02 filter, but never tested in the dedicated analyses (visit 3, relapse screens, per visit), dropped from the dISF-only lists (02b) and the step 09 / 10 sets, and counted as above LOD in step 16 and the relapse follow-up",
  "what 'detected' means", "step 02: >= qc$min_detect_frac in any detection group (serum pools all cohorts); 02b: of all MicroAD samples; step 09: of all dISF samples; step 10: of all MicroAD samples per matrix; dedicated analyses: in one of the two compared groups; relapse targets: in the relapse group - so steps 14 and 17 can call different proteins 'not measurable in serum'",
  "sample QC", "one failing block makes the whole sample FAIL, and the whole sample (all assays) is removed",
  "assay QC", "values with AssayQC FAIL are never removed; with qc$drop_assay_qc_warn = true, WARN values become missing but still count in the detection fractions",
  "missing values", "the main models drop proteins with > 20% missing values and impute the rest; the dedicated analyses do neither",
  "minimum group size", "the main models count samples per group (a relapse group of 1-2 patients with >= 3 samples passes), visits_v1_v4 counts patients; visit 3, the relapse screens and the relapse follow-up use a fixed 3 and step 17 a fixed 2 instead of stats$min_group_n; the relapse screens test V1 non-lesional, V2 and serum V1 without a minimum",
  "firm verdicts", "'>= 10 patients per group for a firm verdict' exists only in key question Q4; the relapse screens state answers with 5 vs 5 patients (backed by the exact relabelling check)",
  "minimum pairs for correlations", "3, 5, 8, 10 or none depending on the step; tnfrsf9_correlation$min_n is used only in step 16 (not for its partial correlation, which needs 8)",
  "minimum patients per visit", "stats$visit_min_subjects applies only to step 13; steps 14 and 17 analyse every visit",
  "before-relapse samples", "visits_v1_v4 leaves relapsers out from their relapse visit on (configurable); steps 14 and 17 use visits before the relapse (fixed); the relapse models of steps 04, 05 and 11 have no before-relapse filter; the V2 'cleared' sets of the relapse screens do not check the skin state",
  "dropouts", "removed from the non-lesional time model but not from the ex-lesional one; step 11 uses the relapse label instead of the dropout column; the AD vs healthy comparisons keep dropouts",
  "gene set size", "5 in step 07 (sets of 5-9 only marked), 10 in steps 09 / 02b and the dedicated GSEA, 5 in step 17, 3 in the relapse follow-up re-test",
  "nominal p", "p < 0.05 is fixed in the code everywhere, while the FDR comes from stats$fdr",
  "fixed texts", "step 18 always says outliers are 'flagged, not removed' (ignores qc$drop_outliers); '50%' and 'FDR < 0.05' are written out in a few texts and ignore the settings",
  "below-LOD cautions", "the caution cut-offs (50%, 80%, 90% / 10%) are not linked to qc$min_detect_frac",
  "skin sites", "step 06 counts healthy skin as non-lesional and includes CPUO patients; steps 10 and 16 exclude CPUO",
  "step 14", "a dISF result at a visit is dropped when the serum model at that visit cannot be estimated")
