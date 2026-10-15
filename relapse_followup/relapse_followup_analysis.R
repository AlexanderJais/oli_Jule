# Relapse follow-up (run with run_relapse_followup.R): eight questions on relapse vs non-relapse, MicroAD AD patients with a
# relapse / non-relapse label (dropouts excluded; the patients of the V1 relapse screen).
#   1 Is relapse vs non-relapse calculated correctly? Labels, direction of log2FC / NES, and the immune gene sets lower in
#     relapsers: the proteins behind them (each tested on its own), overall sample level, positive control (AD lesional vs healthy)
#   2 How the V1 lesion site comparison is calculated, in plain words, with a worked example
#   3 Leukocyte content: PTPRC / CD45          4 T cells: LCK, CD3 and the other selected markers, each on its own; share of CD45
#   5 IL4, CCL27, CCL5, CCL8 over the visits: heatmaps (group means; per patient, clustered)
#   6 / 7 Volcano dISF lesion site V1 / serum V1, every protein with p < 0.05 labelled, with the plotted data
#   8 Druggable targets: higher in relapsers, and higher in lesional and non-lesional AD skin than in healthy skin
# In:  output of steps 01-02 (data/npx_clean.rds); MSigDB / KEGG gene sets; Human Protein Atlas (optional);
#      relapse_followup/ad_drugs.csv (known drugs)
# Out: output/relapse_followup/ - Relapse_followup_executive_summary.pdf, 1_*.xlsx ... 8_*.xlsx with the figures (.pdf / .png)

source("R/utils.R")
source("R/qc_overview.R")
source("R/visit3.R")
source("R/relapse_v1.R")
source("R/relapse_followup.R")
cfg <- load_config()
rf <- cfg$relapse_followup %||% list()
fdr <- cfg$stats$fdr %||% 0.05; min_f <- cfg$qc$min_detect_frac %||% 0.5
colours <- cfg$qc_overview$colours %||% list(above = "#69005F", below = "#FF506E")
p_lab <- rf$volcano_label_p %||% 0.05
g_leuk <- unlist(rf$leukocyte_markers %||% "PTPRC")
g_leuk2 <- unlist(rf$leukocyte_support %||% c("CD48", "ITGB2", "LCP1", "CORO1A", "CD53", "SPN"))
g_t <- unlist(rf$t_cell_markers %||% c("LCK", "CD3E", "CD3D", "CD3G", "CD247", "ZAP70", "LAT", "ITK", "CD2", "CD5", "CD6", "CD28"))
g_act <- unlist(rf$t_cell_activation %||% c("CD69", "IL2RA", "IL2", "IFNG", "TNFRSF9", "TNFRSF4", "ICOS", "CD40LG", "CTLA4",
                                             "PDCD1", "LAG3", "HAVCR2", "GZMB", "CRTAM", "CD70"))
cytokines <- unlist(rf$cytokines %||% c("IL4", "CCL27", "CCL5", "CCL8"))
check_rx <- regex(paste(unlist(rf$check_gene_sets %||% c("T cell receptor", "\\bTCR\\b", "death receptor", "NF.?kappa.?B|NF.?KB",
                                                         "\\bTNF", "MAPK|MAP kinase|kinase cascade|I.?kappa.?B kinase")),
                        collapse = "|"), ignore_case = TRUE)
n_check <- rf$n_check_gene_sets %||% 12
out <- "relapse_followup"
clear_outputs(cfg, out)
set.seed(1)

clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
info <- clean |> distinct(OlinkID, Assay, UniProt)
ad <- clean |> filter(cohort == "MicroAD", group == "AD", relapse %in% c("relapse", "non-relapse"))
samples_where <- function(...) ad |> filter(...) |> distinct(SampleID, SubjectID, relapse) |>
  mutate(grp = factor(relapse, c("non-relapse", "relapse")))
s <- list(V1_lesion = samples_where(matrix == "ISF", visit_num == 1, site == "L"),
          V1_nonlesional = samples_where(matrix == "ISF", visit_num == 1, site == "NL"),
          V2_cleared_lesion = samples_where(matrix == "ISF", visit_num == 2, site == "L", is.na(relapse_visit) | visit_num < relapse_visit),
          serum_V1 = samples_where(matrix == "Serum", visit_num == 1))
nice <- c(V1_lesion = "dISF lesion site V1", V1_nonlesional = "dISF non-lesional V1", V2_cleared_lesion = "dISF cleared lesion site V2",
          serum_V1 = "serum V1")
if (min(table(s$V1_lesion$grp)) < 3) stop("Fewer than 3 patients in a group at the V1 lesion site.", call. = FALSE)
n_rel <- sum(s$V1_lesion$grp == "relapse"); n_non <- sum(s$V1_lesion$grp == "non-relapse")
rel_visits <- ad |> filter(relapse == "relapse") |> distinct(SubjectID, relapse_visit) |> count(relapse_visit) |>
  mutate(t = if_else(is.na(relapse_visit), sprintf("unknown: %d", n), sprintf("V%d: %d", relapse_visit, n))) |> pull(t) |> paste(collapse = ", ")

# ---- differential proteins: relapse vs non-relapse (as in the V1 relapse screen), and AD skin vs healthy skin --------------
res <- imap(s, \(x, nm) de_limma(clean, x, min_f = min_f, fdr = fdr))
walk2(names(res), res, \(nm, r) msg("%s: %d proteins tested, %d with p < 0.05 (higher in relapse %d), %d with FDR < %.2f",
                                    nm, nrow(r), sum(r$p < 0.05), sum(r$p < 0.05 & r$logFC > 0), sum(r$significant), fdr))
hc <- clean |> filter(cohort == "MicroAD", group == "HC", matrix == "ISF") |> distinct(SampleID, SubjectID)
ad_v1 <- \(st) clean |> filter(cohort == "MicroAD", group == "AD", matrix == "ISF", visit_num == 1, site == st) |> distinct(SampleID, SubjectID)
vs_healthy <- \(x, lb) bind_rows(hc |> mutate(grp = "healthy"), x |> mutate(grp = lb)) |> mutate(grp = factor(grp, c("healthy", lb)))
res_h <- if (nrow(hc) >= 3 && nrow(ad_v1("L")) >= 3 && nrow(ad_v1("NL")) >= 3)
  list(lesional_vs_healthy = de_limma(clean, vs_healthy(ad_v1("L"), "AD lesional"), min_f = min_f, fdr = fdr),
       nonlesional_vs_healthy = de_limma(clean, vs_healthy(ad_v1("NL"), "AD non-lesional"), min_f = min_f, fdr = fdr)) else NULL
if (is.null(res_h)) msg("Fewer than 3 healthy, lesional or non-lesional dISF samples: AD skin vs healthy skin is skipped")
n_h <- c(healthy = nrow(hc), lesional = nrow(ad_v1("L")), nonlesional = nrow(ad_v1("NL")))

# per-patient matrices (one sample per patient; columns = patients in the same order in every matrix)
mk_all <- assay_ids(info, unique(c(g_leuk, g_leuk2, g_t, g_act, cytokines)))
pm <- map(s, \(x) { M <- patient_matrix(clean, x, unique(c(unlist(map(res, "OlinkID")), mk_all$OlinkID))); M[, sort(colnames(M)), drop = FALSE] })
rel_of <- \(nm) s[[nm]]$grp[match(colnames(pm[[nm]]), s[[nm]]$SubjectID)] == "relapse"
pct_lod <- map(s, \(x) clean |> filter(SampleID %in% x$SampleID) |> group_by(OlinkID) |>
                 summarise(pct = 100 * mean(!below_lod, na.rm = TRUE), .groups = "drop") |> deframe())
r1 <- res$V1_lesion; M1 <- pm$V1_lesion[r1$OlinkID, , drop = FALSE]; y1 <- rel_of("V1_lesion")
grp_lab <- \(y) factor(if_else(y, "relapse", "non-relapse"), c("non-relapse", "relapse"))
fmt_p <- \(x) ifelse(is.na(x), "-", formatC(x, digits = 2, format = "g"))
sig3 <- \(x) x |> mutate(across(where(is.double), \(v) signif(v, 3)))
rc <- relabel_check(M1, y1)          # exact relabelling of the V1 lesion site (sections 1, 2, 8)

# ---- figures ----------------------------------------------------------------------------------------------------------------
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
save_fig <- function(p, file, width = 11.69, height = 8.27) {
  for (type in c("pdf", "png")) {
    open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
    print(set_family(p, fi$family))
    close_device(fi)
  }
  invisible(p)
}
save_grob <- function(g, file, width = 11.69, height = 8.27) for (type in c("pdf", "png")) {
  open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
  grid::grid.newpage(); grid::grid.draw(g)
  close_device(fi)
}
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
grp_cols <- c(relapse = colours$above, `non-relapse` = "grey55")
# per-patient values (columns SubjectID, group, value, label = panel title), one panel per protein or measure
dots <- function(d, title, subtitle = NULL, ncol = 4, ylab = "NPX") {
  ggplot(d, aes(group, value)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey70") +
    geom_point(aes(colour = group), size = 2.2, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
    scale_colour_manual(values = grp_cols, guide = "none") +
    facet_wrap(~label, scales = "free_y", ncol = ncol) + labs(x = NULL, y = ylab, title = title, subtitle = subtitle) +
    th + theme(strip.text = element_text(size = 7.5))
}
# one panel per protein (rows of pm[[cmp]]), panel title with log2FC, exact p and % above LOD
protein_dots <- function(ids, cmp, title, subtitle, file, ncol = 6) {
  M <- pm[[cmp]]; y <- rel_of(cmp); ids <- intersect(ids, rownames(M))
  if (!length(ids)) return(NULL)
  st <- map(ids, \(o) compare_groups(M[o, ], y)) |> bind_rows() |> mutate(OlinkID = ids, Assay = info$Assay[match(ids, info$OlinkID)])
  d <- map(seq_along(ids), \(i) tibble(SubjectID = colnames(M), value = M[ids[i], ], group = grp_lab(y),
                                       label = sprintf("%s\nlog2FC %+.2f, exact p %s\nabove LOD %.0f%%", st$Assay[i], st$log2FC[i], fmt_p(st$p_exact[i]),
                                                       coalesce(pct_lod[[cmp]][ids[i]], NA_real_)))) |> bind_rows() |>
    mutate(label = factor(label, unique(label)))
  save_fig(dots(d, title, subtitle, ncol), file, height = min(8.27, 1.6 + 2.2 * ceiling(length(ids) / ncol)))
}

# =============================================================================================================================
# 1  Is relapse vs non-relapse calculated correctly?
# =============================================================================================================================
pat <- clean |> filter(cohort == "MicroAD", group == "AD") |>
  distinct(SampleID, SubjectID, matrix, visit_num, site, state, relapse, relapse_raw, time_to_relapse, relapse_visit, volume_ul)
course_txt <- pat |> filter(matrix == "ISF", site == "L") |> arrange(SubjectID, visit_num) |> group_by(SubjectID) |>
  summarise(lesion_site_course = paste0("V", visit_num, " ", coalesce(state, "?"), collapse = ", "),
            V1_lesion_state = state[visit_num %in% 1][1], .groups = "drop")
labels <- pat |> group_by(SubjectID) |>
  summarise(label = paste(sort(unique(na.omit(relapse))), collapse = " / "),
            label_in_manifest = paste(sort(unique(na.omit(relapse_raw))), collapse = " / "),
            time_to_relapse = paste(sort(unique(na.omit(time_to_relapse))), collapse = " / "),
            relapse_visit_from_skin_course = first(relapse_visit), last_visit = suppressWarnings(max(visit_num, na.rm = TRUE)), .groups = "drop") |>
  left_join(course_txt, by = "SubjectID") |>
  mutate(in_comparison = SubjectID %in% s$V1_lesion$SubjectID,
         check = pmap_chr(list(label, relapse_visit_from_skin_course, last_visit, V1_lesion_state), \(l, rv, lv, st) {
           p <- c(if (str_detect(l, " / ")) "different labels on different samples",
                  if (l == "") "no relapse label",
                  if (l == "relapse" && is.na(rv)) "relapse label, but the lesion site does not become lesional again after clearing",
                  if (l == "non-relapse" && !is.na(rv)) sprintf("non-relapse label, but the lesion site becomes lesional again at V%d", rv),
                  if (l == "relapse" && !is.na(rv) && is.finite(lv) && rv != lv) sprintf("relapse at V%d, but samples until V%d", rv, lv),
                  if (!is.na(st) && st != "lesional") sprintf("lesion site at V1 is %s, not lesional", st))
           if (length(p)) paste(p, collapse = "; ") else "ok" })) |>
  arrange(desc(in_comparison), label, SubjectID)
lab_cmp <- labels |> filter(in_comparison)
labels_ok <- all(lab_cmp$check == "ok")

# direction: limma log2FC must equal mean(relapse) - mean(non-relapse) of the raw NPX
raw_d <- rowMeans(M1[, y1, drop = FALSE], na.rm = TRUE) - rowMeans(M1[, !y1, drop = FALSE], na.rm = TRUE)
dir_same <- sum(sign(raw_d) == sign(r1$logFC), na.rm = TRUE); dir_n <- sum(!is.na(raw_d) & !is.na(r1$logFC))
dir_maxdiff <- max(abs(raw_d - r1$logFC), na.rm = TRUE)
direction_ok <- dir_same == dir_n && dir_maxdiff < 0.01

# GSEA as in the V1 relapse screen; does the NES sign agree with the fold changes of the proteins in the set?
sets <- v3_gene_sets(cfg); gmin <- cfg$enrichment$min_size %||% 10; gmax <- cfg$enrichment$max_size %||% 500
g1 <- gsea_keyed(r1, sets, gmin, gmax)
if (nrow(g1)) g1 <- g1 |> mutate(mean_log2FC_measured_members = round(set_mean_logfc(r1, sets, pathway), 3))
g1_sig <- if (nrow(g1)) g1 |> filter(p < 0.05) else tibble()
nes_agree <- if (nrow(g1_sig)) sum(sign(g1_sig$NES) == sign(g1_sig$mean_log2FC_measured_members), na.rm = TRUE) else 0
n_lower_fdr <- if (nrow(g1)) sum(g1$NES < 0 & g1$FDR < fdr) else 0; n_higher_fdr <- if (nrow(g1)) sum(g1$NES > 0 & g1$FDR < fdr) else 0
chk <- if (nrow(g1)) g1 |> filter(NES < 0, p < 0.05, str_detect(gene_set, check_rx)) |> arrange(p) |> slice_head(n = n_check) else tibble()
chk_found <- nrow(chk) > 0
if (!chk_found && nrow(g1)) chk <- g1 |> filter(str_detect(gene_set, check_rx)) |> arrange(p) |> slice_head(n = n_check)
n_imm_lower <- if (nrow(g1)) sum(g1$NES < 0 & g1$p < 0.05 & str_detect(g1$gene_set, check_rx)) else 0
n_imm_higher <- if (nrow(g1)) sum(g1$NES > 0 & g1$p < 0.05 & str_detect(g1$gene_set, check_rx)) else 0

# the checked gene sets re-tested: each sample's overall level removed; AD lesional vs healthy skin (positive control)
lim <- function(M, design, coef) {
  fit <- limma::eBayes(limma::lmFit(M, design))
  tibble(OlinkID = rownames(M), logFC = fit$coefficients[, coef], t = fit$t[, coef], p = fit$p.value[, coef]) |> left_join(info, by = "OlinkID")
}
sens <- list(sample_level_removed = lim(sweep(M1, 2, apply(M1, 2, median, na.rm = TRUE)), cbind(1, y1), 2))
if (!is.null(res_h)) sens$AD_lesional_vs_healthy <- res_h$lesional_vs_healthy
sens_g <- if (nrow(chk)) imap(sens, \(r, nm) {
  g <- gsea_keyed(r, sets[na.omit(chk$pathway)], 3, gmax)
  (if (nrow(g)) g |> transmute(pathway, NES, p) else tibble(pathway = character(), NES = numeric(), p = numeric())) |>
    rename_with(\(v) paste0(v, "_", nm), c(NES, p))
}) else list()
robust <- if (nrow(chk)) reduce(sens_g, \(a, b) left_join(a, b, by = "pathway"),
                                .init = chk |> select(pathway, collection, gene_set, measured_proteins = measured_proteins_in_set, NES, p, FDR,
                                                      mean_log2FC_measured_members, leading_edge)) |>
  mutate(n_leading_edge = map_int(str_split(leading_edge, ", "), length)) |> relocate(leading_edge, .after = last_col()) else tibble()
k_chk <- nrow(robust)
sets_title <- if (chk_found) "Immune signalling gene sets lower in relapsers" else "Immune signalling gene sets (none lower in relapsers with p < 0.05; lowest p shown, either direction)"
k_level <- if (k_chk) sum(robust$NES_sample_level_removed < 0 & robust$p_sample_level_removed < 0.05, na.rm = TRUE) else NA_integer_
k_pos <- if (k_chk && "NES_AD_lesional_vs_healthy" %in% names(robust)) sum(robust$NES_AD_lesional_vs_healthy > 0, na.rm = TRUE) else NA_integer_

# which proteins drive these gene sets: the leading-edge proteins, each tested on its own
le_genes <- if (k_chk) unique(unlist(str_split(robust$leading_edge, ", "))) else character()
le <- r1 |> filter(map_lgl(str_split(Assay, "_"), \(g) any(g %in% le_genes))) |>
  mutate(gene_sets = map_chr(str_split(Assay, "_"), \(g) paste(str_to_sentence(robust$gene_set[map_lgl(str_split(robust$leading_edge, ", "), \(x) any(g %in% x))]), collapse = "; ")),
         cg = map(OlinkID, \(o) compare_groups(M1[o, ], y1)), p_exact = map_dbl(cg, "p_exact"),
         same_sign_without_any_one_patient = map_lgl(cg, "same_sign_without_any_one_patient")) |>
  transmute(Assay, OlinkID, log2FC = logFC, p_limma = p, p_exact, FDR, same_sign_without_any_one_patient, mean_NPX_relapse, `mean_NPX_non-relapse`,
            pct_above_LOD_relapse, `pct_above_LOD_non-relapse`, gene_sets) |> arrange(p_limma)
le_low <- le |> filter(log2FC < 0, p_exact < 0.05)
le_one <- le |> filter(log2FC < 0, p_limma < 0.05, !same_sign_without_any_one_patient)

# overall level of each V1 lesion sample: a lower level of all proteins would push many gene sets down at once
glob <- tibble(SubjectID = colnames(M1), group = grp_lab(y1), median_NPX_all_tested_proteins = apply(M1, 2, median, na.rm = TRUE)) |>
  left_join(clean |> filter(SampleID %in% s$V1_lesion$SampleID) |> group_by(SubjectID) |>
              summarise(pct_all_proteins_above_LOD = 100 * mean(!below_lod, na.rm = TRUE), volume_ul = first(volume_ul),
                        sample_QC_warning = any(SampleQC %in% "WARN"), .groups = "drop"), by = "SubjectID")
glob_tests <- map(c("median_NPX_all_tested_proteins", "pct_all_proteins_above_LOD", "volume_ul"),
                  \(v) compare_groups(as.numeric(glob[[v]]), y1) |> mutate(measure = v, .before = 1)) |> bind_rows()
sev <- read_severity(cfg$paths$severity)
sev_tests <- if (!is.null(sev) && length(severity_scores(sev))) {
  v1s <- sev |> filter(visit == "V1")
  map(severity_scores(sev), \(sc) compare_groups(v1s[[sc]][match(colnames(M1), v1s$SubjectID)], y1) |> mutate(measure = paste("V1", sc), .before = 1)) |> bind_rows()
} else tibble()

q1_verdict <- c(
  sprintf("Labels: %s", if (labels_ok) sprintf("all %d patients have one label that agrees with their skin course (relapsers: lesion site lesional again after clearing; non-relapsers: not).", nrow(lab_cmp))
          else sprintf("%d of %d patients need a look (sheet labels, column check): %s.", sum(lab_cmp$check != "ok"), nrow(lab_cmp),
                       paste(sprintf("%s: %s", lab_cmp$SubjectID[lab_cmp$check != "ok"], lab_cmp$check[lab_cmp$check != "ok"]), collapse = "; "))),
  sprintf("Direction: log2FC and NES are relapse minus non-relapse (positive = higher in relapsers). %s: the limma log2FC equals the difference of the raw group means for %d of %d proteins (largest difference %.4f); the NES sign agrees with the mean fold change of the set's proteins in %d of %d gene sets with p < 0.05.",
          if (direction_ok) "Correct" else "PROBLEM", dir_same, dir_n, dir_maxdiff, nes_agree, nrow(g1_sig)),
  sprintf("Gene sets at the V1 lesion site: %d lower and %d higher in relapsers with FDR < %.2f. Immune signalling sets (TCR, death receptor, NF-kB, TNF, kinase cascades) with p < 0.05: %d lower, %d higher.",
          n_lower_fdr, n_higher_fdr, fdr, n_imm_lower, n_imm_higher),
  if (k_chk) sprintf("Proteins behind %s: the leading edge of the %d checked sets holds %d proteins; %d are lower in relapsers on their own (exact p < 0.05)%s. %s",
                     if (chk_found) "them" else "the closest immune sets (none is lower with p < 0.05)", k_chk, nrow(le), nrow(le_low), if (nrow(le_low)) paste0(": ", paste(head(sprintf("%s %+.2f", le_low$Assay, le_low$log2FC), 15), collapse = ", ")) else "",
                     if (nrow(le_one)) sprintf("For %s the difference changes sign when one patient is left out.", paste(le_one$Assay, collapse = ", ")) else "No protein depends on one single patient.")
  else "No immune signalling gene set was found to check.",
  if (k_chk && !is.na(k_pos)) sprintf("Positive control (AD lesional vs healthy skin, same pipeline): %d of the %d sets are higher in lesional skin%s.", k_pos, k_chk,
                                      if (k_pos == k_chk) ", as expected for immune signalling" else
                                        " - mixed; the direction check above is the direct evidence that the signs are right"),
  if (k_chk) sprintf("After removing each sample's overall level (median NPX), %d of the %d sets are lower in relapsers (p < 0.05).", k_level, k_chk),
  sprintf("Overall level of the V1 lesion samples: median NPX of all tested proteins %+.2f in relapsers (exact p = %s); proteins above LOD %+.1f percentage points (exact p = %s).",
          glob_tests$log2FC[1], fmt_p(glob_tests$p_exact[1]), glob_tests$log2FC[2], fmt_p(glob_tests$p_exact[2])),
  if (nrow(sev_tests)) paste0("Severity at V1, relapse - non-relapse: ", paste(sprintf("%s %+.1f (exact p = %s)", sev_tests$measure, sev_tests$log2FC, fmt_p(sev_tests$p_exact)), collapse = "; "), "."))
q1_short <- paste0(
  if (labels_ok && direction_ok) "Yes, the calculation is correct: labels match the skin course and log2FC / NES are relapse minus non-relapse (checked on every protein). "
  else "A problem was found - see section 1. ",
  if (!k_chk) "No immune signalling gene set could be checked."
  else if (!chk_found) "No immune signalling gene set is lower in relapsers with p < 0.05 in this run."
  else sprintf("The lower immune signalling in relapsers is carried by %d individually lower proteins (exact p < 0.05: %s); after removing each sample's overall level it %s.",
               nrow(le_low), if (nrow(le_low)) paste(head(le_low$Assay, 8), collapse = ", ") else "none",
               if (k_level == k_chk) "holds in all checked gene sets" else if (k_level == 0) "disappears - interpret it with caution"
               else sprintf("holds in %d of %d checked gene sets - interpret with some caution", k_level, k_chk)))

le_wide <- if (nrow(le)) tibble(Assay = le$Assay, as_tibble(round(M1[le$OlinkID, , drop = FALSE], 3))) else tibble()
writexl::write_xlsx(list(
  README = tibble(item = c("question", "answer", paste("finding", seq_along(q1_verdict)), "labels", "direction", "gene_sets_checked", "leading_edge_proteins",
                           "leading_edge_NPX", "sample_level", "GSEA_V1_lesion"),
                  note = c("Was relapse vs non-relapse calculated correctly? Why are whole immune signalling blocks lower in relapsers?", q1_short, q1_verdict,
                           "one row per MicroAD AD patient: relapse label (manifest), the relapse visit derived from the skin course (first lesional visit of the lesion site after it had cleared); check = ok or the problem",
                           "per protein: limma log2FC vs the difference of the raw group means at the V1 lesion site",
                           "immune signalling gene sets lower in relapsers (if none: the closest, either direction): NES and p as in the V1 relapse screen, after removing each sample's overall level (median NPX), and AD lesional vs healthy skin (positive control: expected NES > 0); leading_edge = the proteins that drive the set",
                           "every leading-edge protein of these sets tested on its own: log2FC, limma p, exact relabelling p (difference of the group means over all ways to choose the relapse group), whether one patient drives it, and the gene sets it belongs to",
                           "NPX per patient (columns, V1 lesion site) of the leading-edge proteins",
                           "per patient: median NPX of all tested proteins, % of all proteins above LOD, sample volume - relapse vs non-relapse",
                           "full GSEA of the V1 lesion site (as in the V1 relapse screen) with the mean log2FC of the measured proteins of each set")),
  labels = labels, direction = tibble(OlinkID = r1$OlinkID, Assay = r1$Assay, log2FC_limma = r1$logFC, mean_relapse_minus_mean_nonrelapse = round(raw_d, 4),
                                      same_sign = sign(raw_d) == sign(r1$logFC)),
  gene_sets_checked = if (nrow(robust)) robust |> select(-pathway) |> sig3() else tibble(note = "no immune signalling gene set matched relapse_followup$check_gene_sets"), leading_edge_proteins = le |> sig3(), leading_edge_NPX = le_wide,
  sample_level = glob |> sig3(), sample_level_tests = bind_rows(glob_tests, sev_tests) |> sig3(), GSEA_V1_lesion = if (nrow(g1)) g1 |> select(-pathway) else tibble()),
  out_path(cfg, out, "1_check_relapse_vs_nonrelapse.xlsx"))

p1_le <- if (nrow(le)) protein_dots(head(le$OlinkID, 18), "V1_lesion",
                                    sprintf("1  The proteins behind the %s (leading edge), dISF lesion site V1 (%d relapse vs %d non-relapse)",
                                            if (chk_found) "immune signalling gene sets lower in relapsers" else "closest immune signalling gene sets", n_rel, n_non),
                                    "Each protein on its own, the 18 with the lowest p. NPX per patient", "1_immune_leading_edge_proteins") else NULL
p1_level <- save_fig(glob |> select(SubjectID, group, median_NPX_all_tested_proteins, pct_all_proteins_above_LOD) |>
                       pivot_longer(c(median_NPX_all_tested_proteins, pct_all_proteins_above_LOD), values_to = "value") |>
                       mutate(label = recode(name, median_NPX_all_tested_proteins = "median NPX of all tested proteins", pct_all_proteins_above_LOD = "% of all proteins above LOD")) |>
                       dots("1  Overall level of each V1 lesion sample (a lower level would lower many gene sets at once)",
                            sprintf("median NPX: exact p = %s; %% above LOD: exact p = %s", fmt_p(glob_tests$p_exact[1]), fmt_p(glob_tests$p_exact[2])), ncol = 2, ylab = NULL),
                     "1_sample_level_check", width = 9, height = 5)

# =============================================================================================================================
# 2  How the V1 lesion site comparison is calculated (plain words, worked example)
# =============================================================================================================================
ex <- r1 |> slice_min(p, n = 1, with_ties = FALSE)
ex_v <- tibble(SubjectID = colnames(M1), group = grp_lab(y1), NPX = round(M1[ex$OlinkID, ], 2)) |> arrange(desc(group), desc(NPX))
ex_p_exact <- exact_p_diff(M1[ex$OlinkID, ], y1)
q2_text <- c(
  "## In one sentence",
  sprintf("For every protein we compare the V1 dISF values of the lesion site of the %d patients who relapsed later with those of the %d who did not: how much higher or lower is the relapse group, and could that be chance?", n_rel, n_non),
  "## Step by step",
  sprintf("1. Samples: one dISF sample per patient, from the lesion site at V1 (before any treatment). Relapse = the cleared lesion became lesional again later (%s); non-relapse = it did not until the last visit. Dropouts are left out.", rel_visits),
  sprintf("2. Proteins: %d proteins were tested - those above the detection limit (LOD) in at least %.0f%% of the samples of at least one group.", nrow(r1), 100 * min_f),
  "3. Difference: average NPX of the relapsers minus average NPX of the non-relapsers. NPX is already on a log2 scale, so this difference is the log2 fold change: +1 = twice as high in relapsers, -1 = half as high, 0 = no difference.",
  "4. Is it real or chance? A t-test (limma) compares the difference with the spread inside the groups. limma borrows information on the spread from all proteins, which makes the test steadier with so few patients. Result: the p-value. Small p = a difference this big rarely happens by chance.",
  sprintf("5. Many proteins, many chances: with %d proteins, about %d would get p < 0.05 by luck alone. The FDR corrects for that: here %d proteins have p < 0.05 and %d have FDR < %.2f.",
          nrow(r1), round(0.05 * nrow(r1)), sum(r1$p < 0.05), sum(r1$FDR < fdr), fdr),
  sprintf("6. Shuffle check: there are exactly %d ways to choose %d of the %d patients as 'relapsers'. We repeated the whole test for every one of them. Real labels gave %d proteins higher in relapsers with p < 0.05; random labels gave a median of %.0f. Per protein, the exact p = share of the %d labellings with a difference at least as large (smallest possible: 1/%d).",
          rc$n_labellings, n_rel, n_rel + n_non, rc$observed_hits, median(rc$hits), rc$n_labellings, rc$n_labellings),
  "7. Gene sets (GSEA): all tested proteins are sorted by their t-value, from 'most higher in relapsers' to 'most lower'. A gene set (e.g. T cell receptor signalling) counts as lower in relapsers when its proteins pile up at the 'lower' end. NES < 0 = lower in relapsers, NES > 0 = higher.",
  "## Worked example (the protein with the lowest p)",
  sprintf("%s: mean NPX relapsers %.2f, non-relapsers %.2f; difference (log2FC) %+.2f = %.2f-fold; limma p = %s; exact p (shuffle) = %s. The values of the patients are in the table on the next page.",
          ex$Assay, ex$mean_NPX_relapse, ex$`mean_NPX_non-relapse`, ex$logFC, 2^ex$logFC, fmt_p(ex$p), fmt_p(ex_p_exact)))
q2_short <- sprintf("Per protein: mean V1 lesion-site dISF NPX of the %d relapsers minus that of the %d non-relapsers (= log2 fold change), tested with limma (p, FDR), and checked by repeating the test for all %d ways to pick the relapse group (exact p). GSEA sorts all proteins by the limma t-value.",
                    n_rel, n_non, rc$n_labellings)

# =============================================================================================================================
# 3 / 4  Leukocytes (CD45 / PTPRC) and T cells: each selected marker on its own (config relapse_followup), then summary scores
# =============================================================================================================================
mk <- bind_rows(assay_ids(info, g_leuk) |> mutate(set = "leukocytes (CD45)"), assay_ids(info, g_leuk2) |> mutate(set = "leukocytes (supporting)"),
                assay_ids(info, g_t) |> mutate(set = "T cells"), assay_ids(info, g_act) |> mutate(set = "T-cell activation"))
not_measured <- list(leukocytes = setdiff(c(g_leuk, g_leuk2), mk$gene), t_cells = setdiff(g_t, mk$gene), activation = setdiff(g_act, mk$gene))
marker_res <- imap(s, \(x, nm) {
  M <- pm[[nm]]; y <- rel_of(nm)
  mk |> filter(OlinkID %in% rownames(M)) |>
    mutate(comparison = nice[[nm]], pct_above_LOD = round(pct_lod[[nm]][OlinkID], 1), detected = OlinkID %in% res[[nm]]$OlinkID,
           st = map(OlinkID, \(o) compare_groups(M[o, ], y))) |> unnest(st) |>
    left_join(res[[nm]] |> select(OlinkID, p_limma = p), by = "OlinkID") |> relocate(comparison, set)
}) |> bind_rows()
cd45_id <- head(mk$OlinkID[mk$set == "leukocytes (CD45)"], 1)
# T cells relative to all leukocytes, marker by marker: NPX of the T-cell marker minus NPX of CD45 (log2 ratio)
ratio_res <- imap(s, \(x, nm) {
  M <- pm[[nm]]; y <- rel_of(nm); tm <- mk |> filter(set == "T cells", OlinkID %in% rownames(M))
  if (!length(cd45_id) || !cd45_id %in% rownames(M) || !nrow(tm)) return(tibble())
  map(seq_len(nrow(tm)), \(i) compare_groups(M[tm$OlinkID[i], ] - M[cd45_id, ], y) |>
        mutate(comparison = nice[[nm]], ratio = sprintf("%s / CD45", tm$Assay[i]), OlinkID = tm$OlinkID[i], .before = 1)) |> bind_rows()
}) |> bind_rows()
if (nrow(ratio_res)) ratio_res <- ratio_res |> left_join(marker_res |> distinct(comparison, OlinkID, marker_detected = detected), by = c("comparison", "OlinkID"))
mk_v1 <- \(st) marker_res |> filter(set %in% st, comparison == nice[["V1_lesion"]])
sig_list <- \(d, dir) { x <- d |> filter(coalesce(if ("detected" %in% names(d)) detected else marker_detected, FALSE), p_exact < 0.05, sign(log2FC) == dir)
  if (nrow(x)) paste(sprintf("%s %+.2f", x[[if ("Assay" %in% names(x)) "Assay" else "ratio"]], x$log2FC), collapse = ", ") else "none" }
one_marker <- \(d) sprintf("%s %+.2f (exact p %s, limma p %s, above LOD %.0f%%)", d$Assay, d$log2FC, fmt_p(d$p_exact), fmt_p(d$p_limma), d$pct_above_LOD)
cd45_v1 <- mk_v1("leukocytes (CD45)") |> filter(OlinkID %in% cd45_id) |> slice_head(n = 1); leuk2_v1 <- mk_v1("leukocytes (supporting)")
cd45_ok <- nrow(cd45_v1) == 1 && isTRUE(cd45_v1$detected) && !is.na(cd45_v1$log2FC)
cd45_other <- marker_res |> filter(OlinkID %in% cd45_id, comparison != nice[["V1_lesion"]])
q3_short <- if (!nrow(cd45_v1) || is.na(cd45_v1$log2FC)) sprintf("%s (CD45) is not measured in the dISF data, so leukocyte content cannot be compared with it.", g_leuk[1]) else
  paste0(sprintf("CD45 (%s) at the V1 lesion site: %s in relapsers, log2FC %+.2f (exact p = %s)%s.", cd45_v1$Assay, if (cd45_v1$log2FC < 0) "lower" else "higher",
                 cd45_v1$log2FC, fmt_p(cd45_v1$p_exact), if (!isTRUE(cd45_v1$detected)) " - mostly below LOD, not reliable" else if (isTRUE(cd45_v1$p_exact >= 0.05)) " - no clear difference" else ""),
         sprintf(" Other leukocyte proteins (%s): lower in relapsers: %s; higher: %s.", paste(leuk2_v1$Assay, collapse = ", "), sig_list(leuk2_v1, -1), sig_list(leuk2_v1, 1)))
q3_lines <- c(if (nrow(cd45_v1)) paste("V1 lesion site:", one_marker(cd45_v1)),
              if (nrow(cd45_other)) paste0("CD45 in the other samples: ", paste(sprintf("%s: %+.2f (exact p %s)", cd45_other$comparison, cd45_other$log2FC, fmt_p(cd45_other$p_exact)), collapse = "; "), "."),
              if (nrow(leuk2_v1)) paste0("Other leukocyte proteins, V1 lesion site: ", paste(one_marker(leuk2_v1), collapse = "; "), "."),
              if (length(not_measured$leukocytes)) sprintf("Not in the data: %s.", paste(not_measured$leukocytes, collapse = ", ")))
tm1 <- mk_v1("T cells"); am1 <- mk_v1("T-cell activation"); rt1 <- if (nrow(ratio_res)) ratio_res |> filter(comparison == nice[["V1_lesion"]]) else tibble()
lck <- tm1 |> filter(gene == "LCK") |> slice_head(n = 1)
if (nrow(rt1)) rt1 <- rt1 |> filter(coalesce(marker_detected, FALSE), cd45_ok)        # ratios only of markers and CD45 above LOD
# per marker group: how many detected markers are lower / higher (exact p < 0.05); "consistent" = at least half of them, and at least 2
grp_call <- function(d) {
  d <- d |> filter(coalesce(if ("detected" %in% names(d)) detected else marker_detected, FALSE), !is.na(p_exact))
  lo <- sum(d$p_exact < 0.05 & d$log2FC < 0); hi <- sum(d$p_exact < 0.05 & d$log2FC > 0); n <- nrow(d)
  list(n = n, lo = lo, hi = hi, chance = 0.05 * n,
       call = if (n && lo >= max(2, ceiling(n / 2)) && hi == 0) "lower" else if (n && hi >= max(2, ceiling(n / 2)) && lo == 0) "higher" else "no consistent difference")
}
gt <- grp_call(tm1); ga <- grp_call(am1); gr <- if (nrow(rt1)) grp_call(rt1) else list(n = 0, call = "not available")
q4_short <- paste0(
  sprintf("T-cell markers at the V1 lesion site (%d of %d measured above LOD): lower in relapsers (exact p < 0.05): %s; higher: %s. %s",
          sum(tm1$detected), nrow(tm1), sig_list(tm1, -1), sig_list(tm1, 1),
          if (nrow(lck)) sprintf("LCK: %+.2f (exact p = %s%s).", lck$log2FC, fmt_p(lck$p_exact), if (!lck$detected) ", mostly below LOD" else "") else "LCK is not in the data."),
  if (gr$n) sprintf(" Relative to CD45 (T-cell share of leukocytes): %d of %d markers lower, %d higher (exact p < 0.05).", gr$lo, gr$n, gr$hi)
  else if (!cd45_ok) " The T-cell share (marker minus CD45) is not interpreted: CD45 is not reliably measured." else "",
  sprintf(" Activation markers (%d of %d above LOD): lower: %s; higher: %s. ", sum(am1$detected), nrow(am1), sig_list(am1, -1), sig_list(am1, 1)),
  sprintf("(By chance about %.1f T-cell and %.1f activation markers would reach p < 0.05.) ", gt$chance, ga$chance),
  if (gt$call == "lower" && ga$call == "lower") "Consistent: fewer T-cell proteins and lower activation markers - this fits a less engaged, less activated T-cell response in relapse-prone lesions."
  else if (gt$call == "lower") "T-cell proteins are consistently lower (a less engaged T-cell response), the activation markers not consistently."
  else if (ga$call == "lower") "Activation markers are consistently lower (less activated), the T-cell proteins not consistently."
  else if ("higher" %in% c(gt$call, ga$call)) "This does not fit a weaker T-cell response: T-cell markers are consistently higher in relapsers."
  else "No consistent difference (at least half of a marker group in the same direction): the data neither confirm nor rule out a weaker T-cell response.")
q4_lines <- c(
  if (nrow(tm1)) paste0("T-cell markers: ", paste(one_marker(tm1), collapse = "; "), ".") else "No T-cell marker is in the data.",
  if (nrow(rt1)) paste0("Marker minus CD45 (log2 ratio): ", paste(sprintf("%s %+.2f (exact p %s)", rt1$ratio, rt1$log2FC, fmt_p(rt1$p_exact)), collapse = "; "), "."),
  if (nrow(am1)) paste0("Activation markers: ", paste(one_marker(am1), collapse = "; "), "."),
  if (length(c(not_measured$t_cells, not_measured$activation))) sprintf("Not in the data: %s.", paste(c(not_measured$t_cells, not_measured$activation), collapse = ", ")))
cell_readme <- tibble(item = c("question 3", "answer 3", "question 4", "answer 4", "markers", "T_cells_per_CD45", "NPX_per_patient", "selection", "caution"),
  note = c("Leukocyte content relapse vs non-relapse, PTPRC / CD45 as pan-immune marker", q3_short,
           "T cells relapse vs non-relapse: LCK, CD3 and the other selected T-cell markers, their share of CD45, activation markers", q4_short,
           "each selected protein on its own, per sample set: n, means, log2FC (relapse - non-relapse), Welch t-test p, exact relabelling p (all ways to choose the relapse group), AUC, limma p (as in the volcano), % above LOD; detected = above LOD in >= 50% of at least one group; same_sign_without_any_one_patient = FALSE: one patient drives it",
           "each T-cell marker minus CD45 (log2 ratio = the T-cell share of leukocytes), same statistics; interpreted only when marker and CD45 are above LOD (marker_detected)",
           "NPX per patient of every selected protein, per sample set",
           "the markers are set in config.yml -> relapse_followup (leukocyte_markers, leukocyte_support, t_cell_markers, t_cell_activation)",
           "Olink measures proteins in the fluid (secreted, shed or released), not cells: these are proxies of cell content, not cell counts. LCK and CD3 chains are intracellular / membrane proteins."))
npx_pp <- imap(s, \(x, nm) { M <- pm[[nm]]; ids <- intersect(mk$OlinkID, rownames(M))
  tibble(comparison = nice[[nm]], Assay = rep(info$Assay[match(ids, info$OlinkID)], each = ncol(M)), SubjectID = rep(colnames(M), length(ids)),
         group = rep(as.character(grp_lab(rel_of(nm))), length(ids)), NPX = as.vector(t(M[ids, , drop = FALSE]))) }) |> bind_rows()
writexl::write_xlsx(list(README = cell_readme, markers = marker_res |> sig3(), T_cells_per_CD45 = ratio_res |> sig3(),
                         NPX_per_patient = npx_pp |> sig3()),
                    out_path(cfg, out, "3_4_leukocytes_T_cells.xlsx"))
p3 <- protein_dots(mk$OlinkID[mk$set %in% c("leukocytes (CD45)", "leukocytes (supporting)")], "V1_lesion",
                   "3  Leukocyte content: CD45 (PTPRC) and other leukocyte proteins, dISF lesion site V1", "Each protein on its own; NPX per patient",
                   "3_leukocyte_markers")
p3b <- if (length(cd45_id)) {
  d <- imap(s, \(x, nm) if (cd45_id %in% rownames(pm[[nm]])) {
    r <- marker_res |> filter(OlinkID == cd45_id, comparison == nice[[nm]]) |> slice_head(n = 1)
    tibble(SubjectID = colnames(pm[[nm]]), value = pm[[nm]][cd45_id, ], group = grp_lab(rel_of(nm)),
           label = sprintf("%s\nlog2FC %+.2f, exact p %s\nabove LOD %.0f%%", nice[[nm]], r$log2FC, fmt_p(r$p_exact), r$pct_above_LOD)) }) |> bind_rows() |>
    mutate(label = factor(label, unique(label)))
  save_fig(dots(d, "3  CD45 (PTPRC) in every sample set", "NPX per patient", ncol = 4), "3_CD45_all_samples", height = 4.5)
} else NULL
p4t <- protein_dots(mk$OlinkID[mk$set == "T cells"], "V1_lesion", "4  T-cell markers (LCK, CD3 ...), dISF lesion site V1",
                    "Each marker on its own; NPX per patient. Values below LOD are kept as measured (Olink recommendation): a low % above LOD = less reliable", "4_T_cell_markers")
p4a <- protein_dots(mk$OlinkID[mk$set == "T-cell activation"], "V1_lesion", "4  T-cell activation markers, dISF lesion site V1",
                    "Each marker on its own; NPX per patient. Values below LOD are kept as measured: a low % above LOD = less reliable", "4_T_cell_activation_markers")
p4r <- if (nrow(rt1)) {
  d <- map(seq_len(nrow(rt1)), \(i) tibble(SubjectID = colnames(pm$V1_lesion), value = pm$V1_lesion[rt1$OlinkID[i], ] - pm$V1_lesion[cd45_id, ], group = grp_lab(y1),
                                           label = sprintf("%s\nlog2FC %+.2f, exact p %s", rt1$ratio[i], rt1$log2FC[i], fmt_p(rt1$p_exact[i])))) |> bind_rows() |>
    mutate(label = factor(label, unique(label)))
  save_fig(dots(d, "4  T-cell share of leukocytes: each T-cell marker minus CD45, dISF lesion site V1", "log2 ratio per patient (NPX of the marker - NPX of CD45)",
                ncol = 6, ylab = "NPX marker - NPX CD45"), "4_T_cell_markers_per_CD45", height = min(8.27, 1.6 + 2.2 * ceiling(nrow(rt1) / 6)))
} else NULL
mk_table <- \(d) d |> transmute(marker = Assay, n = sprintf("%d/%d", n_relapse, n_nonrelapse), mean_relapse = round(mean_relapse, 2), mean_non = round(mean_nonrelapse, 2),
                                log2FC = round(log2FC, 2), p_limma = fmt_p(p_limma), p_exact = fmt_p(p_exact), AUC = round(AUC, 2), above_LOD = sprintf("%.0f%%", pct_above_LOD),
                                one_patient = if_else(same_sign_without_any_one_patient %in% FALSE, "yes", ""))
mk_note <- "log2FC = relapse - non-relapse (NPX); p_exact = exact relabelling p; p_limma = as in the volcano (- = not tested, below LOD); one_patient = the sign changes when one patient is left out"

# =============================================================================================================================
# 5  IL4, CCL27, CCL5, CCL8 over the visits
# =============================================================================================================================
cy <- assay_ids(info, cytokines)
q5_missing <- setdiff(cytokines, cy$gene)
cyd <- ad |> filter(OlinkID %in% cy$OlinkID, (matrix == "ISF" & site == "L") | matrix == "Serum", !is.na(visit_num)) |>
  transmute(SubjectID, Assay = factor(Assay, rev(cy$Assay)), matrix = if_else(matrix == "ISF", "dISF lesion site", "serum"), visit_num,
            value, below_lod, group = factor(relapse, c("relapse", "non-relapse")), at_relapse = !is.na(relapse_visit) & visit_num == relapse_visit)
vlev <- paste0("V", sort(unique(cyd$visit_num)))
cyd <- cyd |> mutate(visit = factor(paste0("V", visit_num), vlev))
n_grp <- ad |> distinct(SubjectID, relapse) |> count(relapse) |> deframe()
cm <- cyd |> group_by(matrix, Assay, group, visit) |>
  summarise(mean_NPX = mean(value, na.rm = TRUE), n = sum(!is.na(value)), pct_above_LOD = 100 * mean(!below_lod, na.rm = TRUE),
            at_relapse = sum(at_relapse), .groups = "drop") |>
  group_by(matrix, Assay) |> mutate(centred = mean_NPX - mean(mean_NPX, na.rm = TRUE)) |> ungroup() |>
  mutate(pct_above_LOD = coalesce(pct_above_LOD, 100),
         group = factor(sprintf("%s (%d patients)", group, n_grp[as.character(group)]), sprintf("%s (%d patients)", c("relapse", "non-relapse"), n_grp[c("relapse", "non-relapse")])))
pv <- cyd |> group_by(matrix, Assay, visit) |> group_modify(\(d, k) compare_groups(d$value, d$group == "relapse")) |> ungroup() |>
  mutate(Assay = as.character(Assay), visit = as.character(visit))
p5a <- if (nrow(cm)) save_fig(heatmap_group_means(cm, colours) +
  labs(title = sprintf("5  %s over the visits: relapse vs non-relapse (group means)", paste(cy$Assay, collapse = ", ")),
       subtitle = paste0("Colour = group mean NPX minus the protein's mean (log2); text = mean NPX and number of patients; * = most values below LOD.\n",
                         "Relapsers leave the study at relapse: their last visit is the relapse (lesion site lesional again).")),
  "5_cytokines_heatmap_group_means", width = 11.69, height = 4 + 0.5 * nrow(cy) * n_distinct(cm$matrix)) else NULL
# per patient (dISF lesion site), clustered: rows = patients, columns = each cytokine x visit
ph <- NULL; clus_txt <- "not available"
pp <- cyd |> filter(matrix == "dISF lesion site") |> group_by(Assay) |> mutate(centred = value - mean(value, na.rm = TRUE)) |> ungroup()
if (nrow(pp) && n_distinct(pp$SubjectID) >= 3) {
  cols <- as.vector(t(outer(rev(levels(pp$Assay)), vlev, paste)))
  pts <- sort(unique(pp$SubjectID))
  mat <- matrix(NA_real_, length(pts), length(cols), dimnames = list(pts, cols)); rmk <- matrix("", length(pts), length(cols), dimnames = list(pts, cols))
  key <- cbind(match(pp$SubjectID, pts), match(paste(pp$Assay, pp$visit), cols))
  mat[key] <- pp$centred; rmk[key[pp$at_relapse, , drop = FALSE]] <- "R"
  keep_c <- colSums(!is.na(mat)) > 0; mat <- mat[, keep_c, drop = FALSE]; rmk <- rmk[, keep_c, drop = FALSE]
  grp_pt <- ad |> distinct(SubjectID, relapse) |> filter(SubjectID %in% pts)
  ann_r <- data.frame(group = grp_pt$relapse[match(pts, grp_pt$SubjectID)], row.names = pts)
  cy_col <- str_remove(colnames(mat), " V\\d+$")
  ann_c <- data.frame(protein = factor(cy_col, unique(cy_col)), row.names = colnames(mat))
  dd <- suppressWarnings(dist(mat))
  hc_tree <- if (!anyNA(dd) && nrow(mat) > 2) hclust(dd, "average") else FALSE
  if (inherits(hc_tree, "hclust")) {
    cl <- cutree(hc_tree, 2)
    tab <- table(cluster = cl, group = ann_r$group)
    clus_txt <- paste(sprintf("cluster %s: %s", rownames(tab), apply(tab, 1, \(r) paste(sprintf("%d %s", r, colnames(tab)), collapse = ", "))), collapse = "; ")
  }
  lim5 <- max(abs(mat), na.rm = TRUE)
  ph <- pheatmap::pheatmap(mat, cluster_rows = hc_tree, cluster_cols = FALSE, gaps_col = which(diff(as.integer(ann_c$protein)) != 0),
                           annotation_row = ann_r, annotation_col = ann_c,
                           annotation_colors = list(group = c(relapse = colours$above, `non-relapse` = "grey60")),
                           color = grDevices::colorRampPalette(c(colours$below, "white", colours$above))(101), breaks = seq(-lim5, lim5, length.out = 102),
                           labels_col = str_extract(colnames(mat), "V\\d+$"), display_numbers = rmk, number_color = "black", na_col = "grey92",
                           border_color = "white", fontsize = 9, angle_col = "0",
                           main = sprintf("5  dISF lesion site per patient: %s (NPX minus the protein's mean; R = relapse visit; grey = no sample)", paste(cy$Assay, collapse = ", ")),
                           silent = TRUE)
  save_grob(ph$gtable, "5_cytokines_heatmap_per_patient", width = 11.69, height = 2.5 + 0.3 * nrow(mat))
}
sig_visits <- pv |> filter(p_exact < 0.05)
q5_lines <- c(
  sprintf("Proteins: %s%s. dISF lesion site and serum, every visit; relapsers %d, non-relapsers %d patients.", paste(cy$Assay, collapse = ", "),
          if (length(q5_missing)) sprintf(" (not in the data: %s)", paste(q5_missing, collapse = ", ")) else "", n_grp["relapse"], n_grp["non-relapse"]),
  if (nrow(sig_visits)) paste0("Visits with a difference (exact p < 0.05): ", paste(sprintf("%s %s %s: %s in relapsers (%+.2f, exact p = %s)", sig_visits$Assay, sig_visits$matrix, sig_visits$visit,
                                                                                 if_else(sig_visits$log2FC > 0, "higher", "lower"), sig_visits$log2FC, fmt_p(sig_visits$p_exact)), collapse = "; "), ".")
  else "No visit shows a difference between relapsers and non-relapsers with exact p < 0.05.",
  sprintf("Clustering of the patients (dISF lesion site, all four proteins over the visits, 2 clusters): %s.", clus_txt),
  { lowlod <- cm |> filter(pct_above_LOD < 50) |> distinct(matrix, Assay)
    if (nrow(lowlod)) sprintf("Mostly below LOD (*) in: %s - read these with caution.", paste(unique(sprintf("%s (%s)", lowlod$Assay, lowlod$matrix)), collapse = ", ")) })
q5_short <- paste(if (nrow(sig_visits)) sprintf("%d protein-visit combinations differ (exact p < 0.05): %s.", nrow(sig_visits),
                                               paste(head(sprintf("%s %s %s %s", sig_visits$Assay, sig_visits$visit, if_else(sig_visits$matrix == "serum", "serum", "dISF"),
                                                                  if_else(sig_visits$log2FC > 0, "up", "down")), 6), collapse = ", "))
                  else "No visit differs with exact p < 0.05.", sprintf("Patient clustering: %s.", clus_txt))
writexl::write_xlsx(list(README = tibble(item = c("question", "answer", "group_means", "per_visit_tests", "per_patient"),
                                         note = c("IL4, CCL27, CCL5, CCL8 in relapsers vs non-relapsers over the visits (dISF lesion site and serum)", q5_short,
                                                  "mean NPX per group and visit; centred = mean minus the protein's mean over all cells (the heatmap colour); pct_above_LOD; at_relapse = patients at their relapse visit",
                                                  "relapse vs non-relapse at each visit: n, means, log2FC, Welch t-test p, exact relabelling p, AUC",
                                                  "NPX per patient and visit (the per-patient heatmap); at_relapse = relapse visit")),
                         group_means = cm |> sig3(), per_visit_tests = pv |> sig3(),
                         per_patient = cyd |> arrange(matrix, Assay, group, SubjectID, visit_num) |> select(-visit_num) |> sig3()),
                    out_path(cfg, out, "5_cytokines_over_visits.xlsx"))

# =============================================================================================================================
# 6 / 7  Volcano plots with every protein p < 0.05 labelled, and the plotted data
# =============================================================================================================================
n_ser <- table(s$serum_V1$grp)
vol <- list(dISF = volcano_all_p(r1, sprintf("6  dISF lesion site at V1: %d later relapse vs %d non-relapse - every protein with p < %.2g labelled", n_rel, n_non, p_lab),
                                 colours, p_lab, fdr),
            serum = volcano_all_p(res$serum_V1, sprintf("7  Serum at V1: %d later relapse vs %d non-relapse - every protein with p < %.2g labelled", n_ser[["relapse"]], n_ser[["non-relapse"]], p_lab),
                                  colours, p_lab, fdr))
save_fig(vol$dISF$plot, "6_volcano_dISF_lesion_V1", width = 16.54, height = 11.69)
save_fig(vol$serum$plot, "7_volcano_serum_V1", width = 16.54, height = 11.69)
gsea_answer <- sprintf("Yes. The GSEA of the V1 relapse screen (V1_relapse_GSEA.xlsx, sheet V1_lesion) and of section 1 uses this same limma result: all %d tested proteins - not only the labelled ones - are ranked by column t_used_for_GSEA (the moderated t: its sign is the direction, its size grows with the fold change and the certainty). The volcano shows log2FC and p of the same fit.", nrow(r1))
vol_readme <- \(what, n1, n0, v, extra = NULL) tibble(item = c("plot", "samples", "labels", "columns", extra[[1]]),
  note = c(sprintf("Volcano %s at V1, relapse vs non-relapse (limma); log2FC = relapse - non-relapse", what), sprintf("%d later relapse, %d non-relapse (dropouts excluded)", n1, n0),
           sprintf("every protein with p < %.2g is labelled (%d: %d higher, %d lower in relapse); black ring = FDR < %.2f", p_lab, sum(v$labelled), sum(v$labelled & v$log2FC > 0), sum(v$labelled & v$log2FC < 0), fdr),
           "log2FC, p, minus_log10_p (the y axis), FDR, labelled, class (colour), t_used_for_GSEA, mean NPX and % above LOD per group", extra[[2]]))
r3 <- \(x) x |> mutate(across(c(log2FC, minus_log10_p, t_used_for_GSEA, starts_with("mean_NPX_"), starts_with("pct_above_LOD_")), \(v) signif(v, 3)))   # p, FDR as computed
writexl::write_xlsx(list(README = vol_readme("dISF lesion site", n_rel, n_non, vol$dISF$data, list("used for GSEA?", gsea_answer)), volcano_data = vol$dISF$data |> r3()),
                    out_path(cfg, out, "6_volcano_dISF_lesion_V1_data.xlsx"))
writexl::write_xlsx(list(README = vol_readme("serum", n_ser[["relapse"]], n_ser[["non-relapse"]], vol$serum$data,
                                             list("used for GSEA?", "Yes: the serum GSEA of the V1 relapse screen (sheet serum_V1) ranks all tested proteins by column t_used_for_GSEA of this same fit.")),
                         volcano_data = vol$serum$data |> r3()),
                    out_path(cfg, out, "7_volcano_serum_V1_data.xlsx"))
vol_short <- \(v) sprintf("%d proteins with p < %.2g labelled (%d higher, %d lower in relapsers; FDR < %.2f: %d)", sum(v$labelled, na.rm = TRUE), p_lab,
                         sum(v$labelled & v$log2FC > 0), sum(v$labelled & v$log2FC < 0), fdr, sum(v$FDR_significant, na.rm = TRUE))

# =============================================================================================================================
# 8  Druggable targets: higher in relapsers, and higher in lesional and non-lesional AD skin than in healthy skin
# =============================================================================================================================
hpa <- read_hpa(cfg$paths$hpa)
anno <- if (is.null(hpa)) tibble(OlinkID = character()) else
  assay_protein_classes(info |> filter(OlinkID %in% r1$OlinkID), hpa,
                        classes = list(approved_drug_target = "FDA approved drug targets", potential_drug_target = "Potential drug targets",
                                       secreted = "Predicted secreted proteins", membrane = "Predicted membrane proteins"))$classes |>
  mutate(v = TRUE) |> pivot_wider(names_from = class, values_from = v, values_fill = FALSE)
for (cl in c("approved_drug_target", "potential_drug_target", "secreted", "membrane")) if (!cl %in% names(anno)) anno[[cl]] <- logical(nrow(anno))
drugs <- read_csv("relapse_followup/ad_drugs.csv", show_col_types = FALSE)
drug_of <- info |> mutate(gene = str_split(Assay, "_")) |> unnest(gene) |> inner_join(drugs, by = "gene", relationship = "many-to-many") |>
  group_by(OlinkID) |> summarise(known_drugs = paste(unique(drug), collapse = "; "), drug_action = paste(unique(action), collapse = "; "),
                                 drug_status = paste(unique(status), collapse = "; "), .groups = "drop")
auc <- auc_rows(M1, y1)
# the same AD skin vs healthy comparison with the non-relapsers only: is the protein raised in AD skin without the relapsers?
nr_v1 <- \(st) ad_v1(st) |> filter(SubjectID %in% s$V1_lesion$SubjectID[s$V1_lesion$grp == "non-relapse"])
res_hn <- if (!is.null(res_h) && nrow(nr_v1("L")) >= 3 && nrow(nr_v1("NL")) >= 3)
  list(lesional = de_limma(clean, vs_healthy(nr_v1("L"), "non-relapse lesional"), min_f = min_f, fdr = fdr),
       nonlesional = de_limma(clean, vs_healthy(nr_v1("NL"), "non-relapse non-lesional"), min_f = min_f, fdr = fdr)) else NULL
hn_side <- \(nm, pre) if (is.null(res_hn)) tibble(OlinkID = character(), log2FC = numeric(), p = numeric()) |> rename_with(\(x) paste0(x, "_", pre), -OlinkID) else
  res_hn[[nm]] |> transmute(OlinkID, log2FC = logFC, p) |> rename_with(\(x) paste0(x, "_", pre), -OlinkID)
hside <- \(nm, pre) (if (is.null(res_h)) tibble(OlinkID = character(), logFC = numeric(), p = numeric(), FDR = numeric())
                     else res_h[[nm]] |> select(OlinkID, logFC, p, FDR, starts_with("mean_NPX_"))) |> rename_with(\(x) paste0(x, "_", pre), c(logFC, p, FDR))
tg <- r1 |> transmute(OlinkID, Assay, UniProt, log2FC_relapse_V1_lesion = logFC, p_relapse_V1_lesion = p, FDR_relapse_V1_lesion = FDR,
                      mean_NPX_relapse_V1_lesion = mean_NPX_relapse, mean_NPX_nonrelapse_V1_lesion = `mean_NPX_non-relapse`, pct_above_LOD_relapse) |>
  mutate(p_exact_up = map2_dbl(OlinkID, log2FC_relapse_V1_lesion > 0 & p_relapse_V1_lesion < 0.05,     # one-sided, only for candidates
                               \(o, cand) if (isTRUE(cand)) exact_p_diff(M1[o, ], y1, "greater") else NA_real_),
         AUC = round(auc[OlinkID], 2)) |>
  left_join(res$V1_nonlesional |> select(OlinkID, log2FC_relapse_V1_nonlesional = logFC, p_relapse_V1_nonlesional = p), by = "OlinkID") |>
  left_join(res$V2_cleared_lesion |> select(OlinkID, log2FC_relapse_V2_cleared = logFC, p_relapse_V2_cleared = p), by = "OlinkID") |>
  left_join(hside("lesional_vs_healthy", "lesional_vs_healthy") |> select(-any_of("mean_NPX_AD non-lesional")), by = "OlinkID") |>
  left_join(hside("nonlesional_vs_healthy", "nonlesional_vs_healthy") |> select(-any_of("mean_NPX_healthy")), by = "OlinkID") |>
  left_join(hn_side("lesional", "nonrelapsers_lesional_vs_healthy"), by = "OlinkID") |>
  left_join(hn_side("nonlesional", "nonrelapsers_nonlesional_vs_healthy"), by = "OlinkID") |>
  left_join(anno, by = "OlinkID") |> left_join(drug_of, by = "OlinkID") |>
  mutate(across(c(approved_drug_target, potential_drug_target, secreted, membrane), \(v) coalesce(v, FALSE)),
         up_in_relapse = log2FC_relapse_V1_lesion > 0 & p_relapse_V1_lesion < 0.05,
         up_lesional_vs_healthy = coalesce(logFC_lesional_vs_healthy > 0 & p_lesional_vs_healthy < 0.05, FALSE),
         up_nonlesional_vs_healthy = coalesce(logFC_nonlesional_vs_healthy > 0 & p_nonlesional_vs_healthy < 0.05, FALSE),
         also_up_in_relapse_nonlesional = coalesce(log2FC_relapse_V1_nonlesional > 0 & p_relapse_V1_nonlesional < 0.05, FALSE),
         up_vs_healthy_in_nonrelapsers_too = coalesce(log2FC_nonrelapsers_lesional_vs_healthy > 0 & p_nonrelapsers_lesional_vs_healthy < 0.05 &
                                                        log2FC_nonrelapsers_nonlesional_vs_healthy > 0 & p_nonrelapsers_nonlesional_vs_healthy < 0.05, FALSE),
         druggable = approved_drug_target | potential_drug_target | !is.na(known_drugs),
         tier = case_when(up_in_relapse & up_lesional_vs_healthy & up_nonlesional_vs_healthy & druggable ~ "1 all criteria, druggable",
                          up_in_relapse & up_lesional_vs_healthy & up_nonlesional_vs_healthy ~ "2 all criteria, no drug / drug-target flag",
                          up_in_relapse & (up_lesional_vs_healthy | up_nonlesional_vs_healthy) & druggable ~ "3 higher in relapse and in one AD skin type, druggable",
                          TRUE ~ NA_character_)) |>
  filter(!is.na(tier)) |>
  arrange(tier, desc(p_exact_up < 0.05), desc(also_up_in_relapse_nonlesional), p_exact_up, p_relapse_V1_lesion) |>
  mutate(rank = row_number(), .before = 1)
t1 <- tg |> filter(str_starts(tier, "1")); t2 <- tg |> filter(str_starts(tier, "2")); t3 <- tg |> filter(str_starts(tier, "3"))
q8_short <- if (is.null(res_h)) "AD skin vs healthy skin could not be compared (too few healthy samples)." else if (nrow(t1)) {
  sprintf("%d druggable candidates are higher in relapsers and higher in both lesional and non-lesional skin than in healthy skin: %s.", nrow(t1),
                        paste(head(sprintf("%s (%s)", t1$Assay, coalesce(t1$known_drugs, if_else(t1$approved_drug_target, "approved drug target (HPA)", "potential drug target (HPA)"))), 8), collapse = "; "))
} else {
  sprintf("No druggable protein meets all criteria. Proteins meeting all criteria without a drug flag: %d%s. Druggable proteins higher in relapse and in one of the two AD skin types: %d.",
          nrow(t2), if (nrow(t2)) paste0(" (", paste(head(t2$Assay, 8), collapse = ", "), ")") else "", nrow(t3))
}
q8_lines <- c(
  sprintf("Criteria: (a) higher in relapsers at the V1 lesion site (p < 0.05); (b) higher in AD lesional dISF than in healthy dISF (V1, all %d AD patients vs %d healthy volunteers, p < 0.05); (c) higher in AD non-lesional dISF than in healthy dISF (p < 0.05). Druggable: FDA-approved or potential drug target (Human Protein Atlas%s) or a known drug (relapse_followup/ad_drugs.csv).",
          n_h[["lesional"]], n_h[["healthy"]], if (is.null(hpa)) ": file missing - only the known-drug list was used" else ""),
  "Rationale: a protein that is higher in relapsers before the relapse, and higher in lesional and non-lesional AD skin than in healthy skin, is a candidate whose blocking might lengthen the relapse-free time while sparing healthy skin - to be tested, not shown here.",
  "Note: the AD vs healthy comparison uses all AD patients, relapsers included; column up_vs_healthy_in_nonrelapsers_too says whether the protein is also raised in the non-relapsers' skin alone.",
  if (nrow(t1)) paste0("Tier 1 in detail: ", paste(sprintf("%s (exact p %s%s%s)", head(t1$Assay, 10), fmt_p(head(t1$p_exact_up, 10)),
                                                          if_else(head(t1$also_up_in_relapse_nonlesional, 10), "; also higher in relapsers' non-lesional skin", ""),
                                                          if_else(head(t1$up_vs_healthy_in_nonrelapsers_too, 10), "; raised vs healthy in non-relapsers too", "")), collapse = "; "), "."),
  sprintf("Tier 2 (all criteria, no drug flag): %d. Tier 3 (relapse and one AD skin type, druggable): %d. Full list: 8_druggable_targets.xlsx.", nrow(t2), nrow(t3)))
writexl::write_xlsx(list(README = tibble(item = c("question", "answer", "criteria", "tiers", "columns", "caution"),
                                         note = c("Relapse vs non-relapse, dISF: druggable targets that could lengthen the relapse-free time; only upregulated in lesional and non-lesional AD skin, not in healthy skin",
                                                  q8_short, q8_lines[1],
                                                  "1 = all criteria and druggable; 2 = all criteria, no drug flag; 3 = higher in relapse and in one of the two AD skin types, druggable",
                                                  "p_exact_up = exact relabelling p, one-sided (difference of the group means, higher in relapse); AUC; relapse at V1 non-lesional and V2 cleared lesion site; AD lesional / non-lesional vs healthy (all AD patients: log2FC, p, FDR, means) and the same with non-relapsers only; HPA drug-target classes; known drugs",
                                                  "4 vs 6 patients: candidates to test, not proven targets. Known drugs are a short curated list (check status before use).")),
                         targets = tg |> sig3()),
                    out_path(cfg, out, "8_druggable_targets.xlsx"))
p8 <- {
  top <- head(if (nrow(t1)) t1 else tg, 8)
  if (nrow(top)) {
    d8 <- clean |> filter(cohort == "MicroAD", matrix == "ISF", OlinkID %in% top$OlinkID,
                          group == "HC" | (group == "AD" & visit_num == 1 & site %in% c("L", "NL") & relapse %in% c("relapse", "non-relapse"))) |>
      transmute(Assay = factor(Assay, unique(top$Assay)), value,
                where = factor(case_when(group == "HC" ~ "healthy", site == "NL" ~ paste("non-lesional", relapse), TRUE ~ paste("lesional", relapse)),
                               c("healthy", "non-lesional non-relapse", "non-lesional relapse", "lesional non-relapse", "lesional relapse")))
    save_fig(ggplot(d8, aes(where, value)) + geom_boxplot(outlier.shape = NA, colour = "grey70", width = 0.5) +
               geom_point(aes(colour = str_detect(where, "relapse$") & !str_detect(where, "non-relapse")), size = 1.6, position = position_jitter(width = 0.1, height = 0, seed = 1)) +
               scale_colour_manual(values = c(`TRUE` = colours$above, `FALSE` = "grey50"), guide = "none") +
               facet_wrap(~Assay, scales = "free_y", ncol = 4) + labs(x = NULL, y = "NPX, dISF at V1 (healthy: single visit)",
                    title = if (nrow(t1)) "8  Druggable targets: higher in relapsers and in lesional and non-lesional AD skin than in healthy skin" else "8  No druggable target meets all criteria - the best candidates",
                    subtitle = "violet: patients who relapse later") + th + theme(axis.text.x = element_text(angle = 30, hjust = 1)), "8_druggable_targets")
  } else NULL }

# =============================================================================================================================
# Executive summary PDF
# =============================================================================================================================
source("R/report.R")
pdf_open(out_path(cfg, out, "Relapse_followup_executive_summary.pdf"))
section("Summary", page_text("Relapse vs non-relapse: eight follow-up questions - executive summary", c(
  sprintf("MicroAD, AD patients with a relapse label (dropouts excluded). dISF lesion site at V1: %d later relapse (%s) vs %d non-relapse. log2FC = relapse - non-relapse. With %d vs %d patients, findings are candidates to test in new patients; Olink measures proteins in the fluid, not cells (CD45, LCK, CD3 are proxies of cell content).",
          n_rel, rel_visits, n_non, n_rel, n_non),
  "## 1  Was relapse vs non-relapse calculated correctly?", q1_short,
  "## 2  How was the V1 lesion site comparison calculated?", q2_short,
  "## 3  Leukocyte content (CD45 / PTPRC)", q3_short,
  "## 4  T cells (LCK, CD3 and the other selected markers, each on its own)", q4_short,
  "## 5  IL4, CCL27, CCL5, CCL8 over the visits", q5_short,
  "## 6  Volcano dISF lesion site V1", paste0(vol_short(vol$dISF$data), ". Used for the GSEA? Yes - the GSEA ranks all tested proteins of this same limma result by t."),
  "## 7  Volcano serum V1", paste0(vol_short(vol$serum$data), "."),
  "## 8  Druggable targets (higher in relapse; up in lesional and non-lesional AD skin, not healthy)", q8_short),
  subtitle = sprintf("O-MicroAD, Olink Explore HT, dISF and serum. Generated %s. Details per question on the following pages; tables in the 1_ ... 8_ Excel files.", format(Sys.Date())),
  size = 9, width = 150))
section("1 Check", {
  page_text("1  Was relapse vs non-relapse calculated correctly?", c("## Answer", q1_short, "## Findings", q1_verdict,
    "## What was checked", "Labels: manifest label vs the skin course of the lesion site. Direction: limma log2FC vs raw group means. Immune gene sets lower in relapsers: the proteins that drive them (leading edge), each tested on its own; the sets re-tested after removing each sample's overall level; AD lesional vs healthy skin as positive control."))
  page_table("1  Patients and labels", lab_cmp |> transmute(patient = SubjectID, label, relapse_visit = if_else(is.na(relapse_visit_from_skin_course), "-", paste0("V", relapse_visit_from_skin_course)),
                                                         last_visit = paste0("V", last_visit), check = str_trunc(check, 38)),
             subtitle = "relapse_visit = first visit at which the cleared lesion site is lesional again (from the skin states in the manifest)", rows_per_page = 25)
  if (nrow(robust)) page_table(paste("1 ", sets_title), robust |>
    transmute(gene_set = str_to_sentence(gene_set), source = collection, proteins = measured_proteins, leading_edge = n_leading_edge,
              NES = round(NES, 2), p = fmt_p(p), FDR = fmt_p(FDR), level_NES = round(NES_sample_level_removed, 2), level_p = fmt_p(p_sample_level_removed),
              L_vs_H = if ("NES_AD_lesional_vs_healthy" %in% names(robust)) round(NES_AD_lesional_vs_healthy, 2) else NA),
    subtitle = "NES < 0 = lower in relapsers. leading_edge = number of proteins driving the set (listed in the Excel file); level_: each sample's overall level removed",
    note = "L_vs_H: NES AD lesional vs healthy skin (positive control, expected > 0)", rows_per_page = 20)
  if (nrow(le)) page_table("1  The proteins behind these gene sets (leading edge), each on its own", le |> head(30) |>
    transmute(protein = Assay, log2FC = round(log2FC, 2), p_limma = fmt_p(p_limma), p_exact = fmt_p(p_exact), FDR = fmt_p(FDR),
              one_patient = if_else(same_sign_without_any_one_patient %in% FALSE, "yes", ""), gene_sets = str_trunc(gene_sets, 38)),
    subtitle = "dISF lesion site V1; log2FC = relapse - non-relapse; one_patient = the sign changes when one patient is left out", rows_per_page = 30)
  page_plot(p1_le); page_plot(p1_level) })
section("2 Method", {
  page_text("2  How the dISF lesion site comparison at V1 was calculated - in plain words", q2_text)
  page_table(sprintf("2  Worked example: %s, V1 lesion site, NPX per patient", ex$Assay), ex_v,
             subtitle = sprintf("mean relapse %.2f - mean non-relapse %.2f = log2FC %+.2f; limma p = %s; exact p = %s", ex$mean_NPX_relapse, ex$`mean_NPX_non-relapse`, ex$logFC, fmt_p(ex$p), fmt_p(ex_p_exact)), rows_per_page = 25) })
section("3 Leukocytes", {
  page_text("3  Leukocyte content: CD45 (PTPRC)", c("## Answer", q3_short, "## Details", q3_lines,
    "## Keep in mind", "Olink measures CD45 protein in the fluid (shed or released), not the number of cells: it is a proxy of leukocyte content. Full statistics: 3_4_leukocytes_T_cells.xlsx."))
  page_table("3  Leukocyte proteins, dISF lesion site V1 (each on its own)", mk_table(mk_v1(c("leukocytes (CD45)", "leukocytes (supporting)"))), subtitle = mk_note, rows_per_page = 30)
  page_plot(p3); page_plot(p3b) })
section("4 T cells", { page_text("4  T cells: LCK, CD3 and the other selected markers", c("## Answer", q4_short, "## Details", q4_lines,
    "## How to read it", "Each marker is shown and tested on its own (log2FC, exact p). Marker / CD45 = NPX of the T-cell marker minus NPX of CD45: the T-cell share of all leukocytes. A marker group counts as lower only if at least half of its detected markers (and at least 2) are lower with exact p < 0.05 and none higher."))
  page_table("4  T-cell and activation markers, dISF lesion site V1 (each on its own)",
             bind_rows(mk_v1("T cells") |> mutate(Assay = paste0(Assay, " (T cell)")), mk_v1("T-cell activation") |> mutate(Assay = paste0(Assay, " (activation)"))) |> mk_table(),
             subtitle = mk_note, rows_per_page = 30)
  page_plot(p4t); page_plot(p4r); page_plot(p4a) })
section("5 Cytokines", { page_text("5  IL4, CCL27, CCL5, CCL8 over the visits", c("## Answer", q5_short, "## Details", q5_lines,
    "## Figures", "Next pages: group means per visit (y = protein, x = visit, relapse and non-relapse side by side), then each patient clustered (rows = patients, columns = each protein over the visits). Values per visit: 5_cytokines_over_visits.xlsx."))
  page_plot(p5a); if (!is.null(ph)) { grid::grid.newpage(); grid::grid.draw(ph$gtable) }
  page_table("5  Relapse vs non-relapse per visit", pv |> transmute(matrix, protein = as.character(Assay), visit, n = sprintf("%d/%d", n_relapse, n_nonrelapse),
                                                                    log2FC = round(log2FC, 2), p_t = fmt_p(p_t), p_exact = fmt_p(p_exact)) |> arrange(matrix, protein, visit),
             subtitle = "n = relapse / non-relapse patients with a sample at that visit", rows_per_page = 30) })
section("6 Volcano dISF", { page_plot(vol$dISF$plot)
  page_text("6  Volcano dISF lesion site V1: the data and the GSEA", c(vol_short(vol$dISF$data), "## Was the volcano data used for the GSEA?", gsea_answer,
    "The plotted data (every protein, with the labelled ones marked) are in 6_volcano_dISF_lesion_V1_data.xlsx.")) })
section("7 Volcano serum", { page_plot(vol$serum$plot)
  page_text("7  Volcano serum V1", c(vol_short(vol$serum$data), "The plotted data are in 7_volcano_serum_V1_data.xlsx; the serum GSEA of the V1 relapse screen uses the t column of the same fit.")) })
section("8 Targets", { page_text("8  Druggable targets to lengthen the relapse-free time", c("## Answer", q8_short, "## Details", q8_lines))
  page_table("8  Candidate targets (dISF)", tg |> head(30) |>
               transmute(rank, protein = Assay, tier = str_sub(tier, 1, 1), log2FC_rel = round(log2FC_relapse_V1_lesion, 2), p_rel = fmt_p(p_relapse_V1_lesion), exact_p = fmt_p(p_exact_up),
                         NL_rel = if_else(also_up_in_relapse_nonlesional, "yes", ""), L_vs_H = round(logFC_lesional_vs_healthy, 2), NL_vs_H = round(logFC_nonlesional_vs_healthy, 2),
                         drug = str_trunc(coalesce(known_drugs, if_else(approved_drug_target, "approved target (HPA)", if_else(potential_drug_target, "potential target (HPA)", ""))), 32)),
             subtitle = "tier 1 = all criteria and druggable; 2 = all criteria, no drug flag; 3 = relapse and one AD skin type, druggable. NL_rel: also higher in relapsers' non-lesional skin",
             note = "L_vs_H / NL_vs_H: log2FC AD lesional / non-lesional vs healthy dISF. Known drugs: relapse_followup/ad_drugs.csv (curated, check before use)", rows_per_page = 30)
  page_plot(p8) })
section("Methods", page_text("Methods", c(
  "Data: Olink Explore HT, PCNormalizedNPX; MicroAD AD patients with a relapse / non-relapse label (dropouts excluded), the same patients and tests as the V1 relapse screen (run_relapse_v1.R).",
  sprintf("Differential tests: limma per protein (proteins above LOD in >= %.0f%% of the samples of at least one group), relapse vs non-relapse at the V1 lesion site, V1 non-lesional skin, V2 cleared lesion site and V1 serum; AD lesional / non-lesional (V1) vs healthy dISF for section 8.", 100 * min_f),
  "Exact relabelling: the test repeated for all ways to choose the relapse group; exact p = share of labellings with a difference at least as large.",
  "Sections 3, 4: each selected marker (config relapse_followup) tested on its own: Welch t-test, exact relabelling p (difference of the group means over all ways to choose the relapse group; the same definition in all sections), AUC, leave-one-patient-out sign stability, and the limma p of the volcano. Marker / CD45 = NPX difference (log2 ratio).",
  "Section 1: GSEA (fgsea) on the limma t as in the V1 relapse screen; gene sets matching the immune signalling terms (config relapse_followup$check_gene_sets) re-tested after subtracting each sample's median NPX and in AD lesional vs healthy skin; their leading-edge proteins tested one by one.",
  "Section 5: group means per visit; per-patient heatmap clustered by average linkage on the centred NPX of all four proteins over the visits.",
  sprintf("Sections 6, 7: volcano of the limma result; every protein with p < %.2g labelled.", p_lab),
  "Section 8: Human Protein Atlas drug-target classes (if downloaded) and a short curated list of drugs relevant to AD / type 2 inflammation.")))
invisible(grDevices::dev.off())
msg("Relapse follow-up: %s", file.path(cfg$paths$output, out))
