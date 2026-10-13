# V1 relapse screen (run with run_relapse_v1.R): which dISF proteins are already higher at baseline (V1, before the
# topical steroids at V2) in the AD patients who later relapse? MicroAD, AD patients with a relapse / non-relapse label
# (dropouts excluded). Groups are tiny (4 vs 6), so the list is a ranked screen with an exact relabelling check.
#   V1 lesion site (primary), V1 non-lesional skin, V1 lesion minus non-lesional (same patient),
#   V2 cleared lesion site (still higher after treatment?), V1 serum (also visible in blood?)
# In:  output of steps 01-02 (data/npx_clean.rds); Human Protein Atlas (drug targets, location); MSigDB / KEGG gene sets
# Out: output/relapse_V1/ - V1_relapse_targets.xlsx, V1_relapse_GSEA.xlsx, V1_relapse_volcano.pdf/.png,
#      V1_relapse_targets_course.pdf/.png, V1_relapse_chance_check.pdf/.png, V1_relapse_GSEA.pdf/.png

source("R/utils.R")
source("R/qc_overview.R")
source("R/visit3.R")
source("R/relapse_v1.R")
cfg <- load_config()
rv <- cfg$relapse_v1 %||% list()
n_lab <- rv$top_labels %||% 25; n_course <- rv$n_course %||% 12
fdr <- cfg$stats$fdr %||% 0.05; min_f <- cfg$qc$min_detect_frac %||% 0.5
colours <- cfg$qc_overview$colours %||% list(above = "#69005F", below = "#FF506E")
out <- "relapse_V1"
clear_outputs(cfg, out)
set.seed(1)

clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
ad <- clean |> filter(cohort == "MicroAD", group == "AD", relapse %in% c("relapse", "non-relapse"))
samples_where <- function(...) ad |> filter(...) |> distinct(SampleID, SubjectID, relapse, plate) |>
  mutate(grp = factor(relapse, c("non-relapse", "relapse")))
s <- list(V1_lesion = samples_where(matrix == "ISF", visit_num == 1, site == "L"),
          V1_nonlesional = samples_where(matrix == "ISF", visit_num == 1, site == "NL"),
          V2_cleared_lesion = samples_where(matrix == "ISF", visit_num == 2, site == "L", is.na(relapse_visit) | visit_num < relapse_visit),
          serum_V1 = samples_where(matrix == "Serum", visit_num == 1))
n_txt <- imap_chr(s, \(x, nm) sprintf("%s: %d relapse, %d non-relapse", nm, sum(x$grp == "relapse"), sum(x$grp == "non-relapse")))
walk(n_txt, \(x) msg(x))
if (min(table(s$V1_lesion$grp)) < 3) stop("Fewer than 3 patients in a group at the V1 lesion site.", call. = FALSE)
rel_visits <- ad |> filter(relapse == "relapse") |> distinct(SubjectID, relapse_visit) |> count(relapse_visit) |>
  mutate(t = paste0("V", relapse_visit, ": ", n)) |> pull(t) |> paste(collapse = ", ")

# lesion minus non-lesional skin of the same patient at V1 (detection judged at the lesion site)
delta <- clean |> filter(SampleID %in% s$V1_lesion$SampleID) |> select(SubjectID, OlinkID, Assay, UniProt, L = value, below_lod) |>
  inner_join(clean |> filter(SampleID %in% s$V1_nonlesional$SampleID) |> select(SubjectID, OlinkID, NL = value), by = c("SubjectID", "OlinkID")) |>
  transmute(SampleID = paste0(SubjectID, "_LminusNL"), SubjectID, OlinkID, Assay, UniProt, value = L - NL, below_lod)
s$V1_lesion_minus_nonlesional <- s$V1_lesion |> mutate(SampleID = paste0(SubjectID, "_LminusNL"))

# ---- differential proteins per comparison ----------------------------------------------------------------------------
res <- imap(s, \(x, nm) de_limma(if (nm == "V1_lesion_minus_nonlesional") delta else clean, x, min_f = min_f, fdr = fdr))
walk2(names(res), res, \(nm, r) msg("%s: %d proteins tested, %d with p < 0.05 (higher in relapse %d), %d with FDR < %.2f",
                                    nm, nrow(r), sum(r$p < 0.05), sum(r$p < 0.05 & r$logFC > 0), sum(r$significant), fdr))

# ---- the screen on the V1 lesion site: exact relabelling check, leave-one-out, AUC ------------------------------------
v1 <- res$V1_lesion
M <- patient_matrix(clean, s$V1_lesion, v1$OlinkID)
is_rel <- s$V1_lesion$grp[match(colnames(M), s$V1_lesion$SubjectID)] == "relapse"
rc <- relabel_check(M, is_rel)
loo <- leave_one_out(M, is_rel)
auc <- auc_rows(M, is_rel)
chance <- tibble(item = c("patients", "relabellings (all ways to choose the relapse group)", "proteins higher in relapse with p < 0.05 (observed)",
                          "the same under random labels: median (90th percentile)", "share of random labellings with at least as many",
                          "estimated chance hits among the observed"),
                 value = c(sprintf("%d relapse, %d non-relapse", sum(is_rel), sum(!is_rel)), as.character(rc$n_labellings), as.character(rc$observed_hits),
                           sprintf("%.0f (%.0f)", median(rc$hits), quantile(rc$hits, 0.9)), sprintf("%.2f", mean(rc$hits >= rc$observed_hits)),
                           sprintf("about %.0f of %d", min(median(rc$hits), rc$observed_hits), rc$observed_hits)))
msg("Chance check: %d proteins higher in relapse (p < 0.05); random labels give a median of %.0f; P = %.2f",
    rc$observed_hits, median(rc$hits), mean(rc$hits >= rc$observed_hits))

# ---- can V1 dISF predict relapse? (honest test of the whole marker search, and per-marker prediction) -------------------
combs <- utils::combn(length(is_rel), sum(is_rel))
lab_of <- \(j) seq_along(is_rel) %in% combs[, j]
ps_obs <- perfect_separators(M, matrix(which(is_rel)))
ps_null <- perfect_separators(M, combs)
pred_k <- unlist(rv$predict_k %||% c(1, 5))
nested <- map(pred_k, \(k) {
  obs <- nested_loo_correct(M, is_rel, k)
  null <- map_int(seq_len(ncol(combs)), \(j) nested_loo_correct(M, lab_of(j), k))
  list(k = k, obs = obs, null = null, p = mean(null >= obs))
})
walk(nested, \(x) msg("Prediction with the top %d protein(s), chosen without the predicted patient: %d of %d patients correct; random labels: median %.0f; P = %.3f",
                      x$k, x$obs, length(is_rel), median(x$null), x$p))
predict_tbl <- tibble(item = c(sprintf("proteins separating relapse and non-relapse perfectly (AUC 1 or 0)"),
                               "the same under random labels: median (90th percentile)",
                               map_chr(nested, \(x) sprintf("patients predicted correctly with the top %d protein(s), picked on the other patients (leave one out)", x$k)),
                               map_chr(nested, \(x) sprintf("the same under random labels (top %d): median; P (at least as many correct)", x$k))),
                      value = c(as.character(ps_obs), sprintf("%.0f (%.0f)", median(ps_null), quantile(ps_null, 0.9)),
                                map_chr(nested, \(x) sprintf("%d of %d", x$obs, length(is_rel))),
                                map_chr(nested, \(x) sprintf("%.0f; P = %.3f", median(x$null), x$p))))
time_score <- s$V1_lesion |> distinct(SubjectID, relapse) |>
  left_join(ad |> distinct(SubjectID, relapse_visit), by = "SubjectID") |>
  mutate(t = if_else(relapse == "relapse", relapse_visit, max(ad$visit_num, na.rm = TRUE) + 1))
tvec <- time_score$t[match(colnames(M), time_score$SubjectID)]
rho_time <- apply(M, 1, \(v) suppressWarnings(cor(v, tvec, method = "spearman", use = "pairwise.complete.obs")))
p_time <- apply(M, 1, \(v) tryCatch(suppressWarnings(cor.test(v, tvec, method = "spearman", exact = FALSE)$p.value), error = \(e) NA_real_))
per_marker <- map(rownames(M), \(oid) loo_predict_marker(M[oid, ], is_rel) |> mutate(OlinkID = oid)) |> bind_rows() |>
  mutate(rho_time_to_relapse = round(rho_time[OlinkID], 2), p_time_to_relapse = signif(p_time[OlinkID], 3))

# ---- annotation: drug targets and location (Human Protein Atlas) --------------------------------------------------------
hpa <- read_hpa(cfg$paths$hpa)
anno <- if (is.null(hpa)) tibble(OlinkID = character()) else {
  a <- assay_protein_classes(v1 |> distinct(OlinkID, Assay, UniProt), hpa,
                             classes = list(approved_drug_target = "FDA approved drug targets", potential_drug_target = "Potential drug targets",
                                            secreted = "Predicted secreted proteins", membrane = "Predicted membrane proteins"))$classes
  a |> mutate(v = TRUE) |> pivot_wider(names_from = class, values_from = v, values_fill = FALSE)
}
for (cl in c("approved_drug_target", "potential_drug_target", "secreted", "membrane")) if (!cl %in% names(anno)) anno[[cl]] <- logical(nrow(anno))

# same-plate check: the relapsers' V1 lesion samples may share one Olink plate; a target should not simply be higher on it
rel_plate <- names(which.max(table(s$V1_lesion$plate[s$V1_lesion$grp == "relapse"])))
nr <- s$V1_lesion |> filter(grp == "non-relapse") |> mutate(same = plate == rel_plate)
plate_check <- function(oid) {
  if (sum(nr$same) < 2 || sum(!nr$same) < 2 || n_distinct(s$V1_lesion$plate[s$V1_lesion$grp == "relapse"]) > 1) return("not needed / not testable")
  v <- M[oid, nr$SubjectID]; p <- tryCatch(t.test(v[nr$same], v[!nr$same])$p.value, error = \(e) NA)
  if (!is.na(p) && p < 0.05 && mean(v[nr$same], na.rm = TRUE) > mean(v[!nr$same], na.rm = TRUE)) "also higher in non-relapsers on that plate" else "ok"
}

# ---- ranked lists -------------------------------------------------------------------------------------------------------
side <- \(nm, pre) res[[nm]] |> select(OlinkID, logFC, p) |> rename_with(\(x) paste0(x, "_", pre), c(logFC, p))
scr <- v1 |> left_join(rc$per_protein |> select(OlinkID, p_exact_up, p_exact_down, p_exact_two_sided), by = "OlinkID") |>
  left_join(loo, by = "OlinkID") |> mutate(AUC = round(auc[OlinkID], 2)) |> left_join(per_marker, by = "OlinkID") |>
  left_join(side("V1_nonlesional", "V1_nonlesional"), by = "OlinkID") |>
  left_join(side("V1_lesion_minus_nonlesional", "V1_LminusNL"), by = "OlinkID") |>
  left_join(side("V2_cleared_lesion", "V2_cleared"), by = "OlinkID") |>
  left_join(side("serum_V1", "serum_V1"), by = "OlinkID") |>
  left_join(anno, by = "OlinkID") |>
  mutate(across(c(approved_drug_target, potential_drug_target, secreted, membrane), \(v) coalesce(v, FALSE)),
         persists_V2_up = coalesce(logFC_V2_cleared > 0 & p_V2_cleared < 0.05, FALSE),
         persists_V2_down = coalesce(logFC_V2_cleared < 0 & p_V2_cleared < 0.05, FALSE),
         also_serum_up = coalesce(logFC_serum_V1 > 0 & p_serum_V1 < 0.05, FALSE),
         also_serum_down = coalesce(logFC_serum_V1 < 0 & p_serum_V1 < 0.05, FALSE),
         detected_in_relapse = `pct_above_LOD_relapse` >= 100 * min_f,
         detected_in_nonrelapse = `pct_above_LOD_non-relapse` >= 100 * min_f,
         candidate_up = logFC > 0 & p < 0.05 & p_exact_up < 0.05 & detected_in_relapse,
         candidate_down = logFC < 0 & p < 0.05 & p_exact_down < 0.05 & detected_in_nonrelapse)
cand_up <- scr |> filter(candidate_up) |> arrange(desc(loo_correct), desc(AUC), desc(persists_V2_up), desc(loo_always_higher), desc(logFC))
cand_up <- cand_up |> mutate(rank = row_number(), same_plate_check = map_chr(OlinkID, plate_check))
cand_down <- scr |> filter(candidate_down) |> arrange(desc(loo_correct), AUC, desc(persists_V2_down), logFC) |> mutate(rank = row_number(), same_plate_check = map_chr(OlinkID, plate_check))
tcols <- c("rank", "Assay", "OlinkID", "UniProt", "loo_correct", "loo_n", "loo_sensitivity", "loo_specificity", "cutoff_NPX", "AUC",
           "rho_time_to_relapse", "p_time_to_relapse", "logFC", "p", "p_exact_up", "p_exact_down", "loo_always_higher", "loo_max_p",
           "persists_V2_up", "persists_V2_down", "logFC_V2_cleared", "p_V2_cleared", "logFC_V1_nonlesional", "p_V1_nonlesional",
           "logFC_V1_LminusNL", "p_V1_LminusNL", "also_serum_up", "also_serum_down", "logFC_serum_V1", "p_serum_V1",
           "approved_drug_target", "potential_drug_target", "secreted", "membrane",
           "mean_NPX_relapse", "mean_NPX_non-relapse", "pct_above_LOD_relapse", "pct_above_LOD_non-relapse", "FDR", "same_plate_check")
msg("Targets: %d higher, %d lower in relapsers at V1 (lesion site)", nrow(cand_up), nrow(cand_down))

answer <- if (any(map_lgl(nested, \(x) x$p < 0.05))) "V1 dISF predicts relapse better than chance in these patients (sheet can_V1_predict); the ranked candidates are in targets_up / targets_down" else
  "With these patients, predicting relapse from V1 dISF is not better than chance once the search over all proteins is accounted for (sheet can_V1_predict); the targets are candidates to test in new patients, not proven predictors"
readme <- tibble(item = c("question", "answer", "patients", "relapse visits", "comparisons", "target (higher)", "ranking", "prediction columns", "flags", "chance check", "caution"),
  note = c("Which dISF markers at V1 (baseline, before the topical steroids) predict a later relapse?", answer,
           paste(n_txt, collapse = "; "), rel_visits,
           "limma per protein, relapse vs non-relapse: V1 lesion site (primary), V1 non-lesional skin, V1 lesion minus non-lesional (same patient), V2 cleared lesion site (before any relapse), V1 serum. log2FC = relapse - non-relapse",
           sprintf("V1 lesion site: log2FC > 0, p < 0.05, exact relabelling p < 0.05 (all %d ways to choose the relapse group), above LOD in >= %.0f%% of the relapsers", rc$n_labellings, 100 * min_f),
           "patients predicted correctly when each is left out (loo_correct), then AUC, then still higher at V2 after treatment, then higher in every leave-one-patient-out refit, then fold change",
           "loo_correct / loo_sensitivity / loo_specificity: each patient classified with a cut-off set on the other patients (midpoint of the group means); cutoff_NPX: that cut-off on all patients (V1 lesion site, NPX); rho_time_to_relapse: Spearman with the relapse visit (non-relapse = after the last visit), < 0 = higher value, earlier relapse",
           "persists_V2: still higher in the cleared lesion site at V2; also_serum: also higher in V1 serum; approved / potential drug target, secreted, membrane: Human Protein Atlas",
           "sheet chance_check: how many proteins random labels would give; if the observed number is not above that, the list is mostly chance",
           "4 vs 6 patients: a screen to pick candidates, not a confirmation; the smallest possible exact p is 1/210"))
writexl::write_xlsx(list(README = readme, can_V1_predict = predict_tbl, chance_check = chance,
                         targets_up = cand_up |> select(any_of(tcols)) |> mutate(across(where(is.double), \(v) signif(v, 3))),
                         targets_down = cand_down |> select(any_of(tcols)) |> mutate(across(where(is.double), \(v) signif(v, 3))),
                         V1_lesion = res$V1_lesion, V1_nonlesional = res$V1_nonlesional, V1_lesion_minus_nonlesional = res$V1_lesion_minus_nonlesional,
                         V2_cleared_lesion = res$V2_cleared_lesion, serum_V1 = res$serum_V1),
                    out_path(cfg, out, "V1_relapse_targets.xlsx"))

# ---- GSEA ---------------------------------------------------------------------------------------------------------------
sets <- v3_gene_sets(cfg)
gsea <- map(res[c("V1_lesion", "V1_lesion_minus_nonlesional", "V2_cleared_lesion")],
            \(r) run_gsea(r, sets, cfg$enrichment$min_size %||% 10, cfg$enrichment$max_size %||% 500))
writexl::write_xlsx(c(list(README = tibble(note = "GSEA (fgsea) on the limma t statistics, relapse vs non-relapse; NES > 0 = higher in relapsers. MSigDB collections of config enrichment$collections and KEGG")),
                      gsea), out_path(cfg, out, "V1_relapse_GSEA.xlsx"))

# ---- figures --------------------------------------------------------------------------------------------------------------
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
save_fig <- function(p, file, width = 11.69, height = 8.27) for (type in c("pdf", "png")) {
  open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
  print(set_family(p, fi$family))
  close_device(fi)
}
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
vol <- scr |> mutate(cl = case_when(candidate_up ~ "higher in relapse (target)", candidate_down ~ "lower in relapse", TRUE ~ "other"))
lab <- vol |> arrange(desc(candidate_up | candidate_down), p) |> slice_head(n = n_lab)
save_fig(ggplot(vol, aes(logFC, -log10(p))) +
           geom_point(aes(colour = cl), size = 1.1) +
           ggrepel::geom_text_repel(data = lab, aes(label = Assay), size = 2.6, max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60", segment.size = 0.3) +
           scale_colour_manual(values = c(`higher in relapse (target)` = colours$above, `lower in relapse` = colours$below, other = "grey75"), name = NULL) +
           labs(x = "log2 fold change (relapse - non-relapse)", y = "-log10 p",
                title = sprintf("dISF lesion site at V1: %d later relapse vs %d non-relapse - %d proteins higher, %d lower (targets)",
                                sum(is_rel), sum(!is_rel), nrow(cand_up), nrow(cand_down)),
                subtitle = sprintf("Target: p < 0.05 and exact relabelling p < 0.05. %d proteins higher with p < 0.05; random labels give a median of %.0f. Labelled: %d proteins (targets first, then lowest p)",
                                   rc$observed_hits, median(rc$hits), n_lab)) + th,
         "V1_relapse_volcano")

show <- (if (nrow(cand_up)) cand_up else scr |> filter(logFC > 0) |> arrange(p)) |> slice_head(n = n_course)
course <- clean |> filter(cohort == "MicroAD", matrix == "ISF", site == "L", SubjectID %in% s$V1_lesion$SubjectID, OlinkID %in% show$OlinkID) |>
  select(SubjectID, OlinkID, Assay, visit_num, value, relapse, relapse_visit) |>
  mutate(Assay = factor(Assay, unique(show$Assay)), at_relapse = !is.na(relapse_visit) & visit_num == relapse_visit)
save_fig(ggplot(course, aes(visit_num, value, group = SubjectID, colour = relapse)) +
           geom_line(alpha = 0.8) + geom_point(size = 1.2) +
           geom_point(data = \(d) filter(d, at_relapse), shape = 4, size = 3, stroke = 1.1, colour = "black") +
           scale_colour_manual(values = c(relapse = colours$above, `non-relapse` = "grey60"), name = NULL) +
           scale_x_continuous(breaks = 1:6, labels = paste0("V", 1:6)) + facet_wrap(~Assay, scales = "free_y", ncol = 4) +
           labs(x = NULL, y = "NPX, dISF lesion site",
                title = if (nrow(cand_up)) sprintf("Top %d targets: dISF lesion site per patient over the visits", nrow(show)) else
                  sprintf("No target passed - the %d proteins with the lowest p (higher in relapse)", nrow(show)),
                subtitle = "violet: patients who relapse (x = relapse visit); grey: non-relapse. V1 baseline, V2 after topical steroids") + th,
         "V1_relapse_targets_course")

save_fig(ggplot(tibble(hits = rc$hits), aes(hits)) + geom_histogram(binwidth = 1, fill = "grey75") +
           geom_vline(xintercept = rc$observed_hits, colour = colours$above, linewidth = 1) +
           labs(x = "proteins higher in the 'relapse' group with p < 0.05", y = "relabellings",
                title = sprintf("Chance check: observed %d (violet) vs all %d ways to choose %d of the %d patients as relapse",
                                rc$observed_hits, rc$n_labellings, sum(is_rel), length(is_rel)),
                subtitle = sprintf("%.0f%% of the random labellings give at least as many proteins", 100 * mean(rc$hits >= rc$observed_hits))) + th,
         "V1_relapse_chance_check", width = 9, height = 5)
pn <- map(nested, \(x) tibble(k = sprintf("top %d protein(s)", x$k), correct = x$null, obs = x$obs)) |> bind_rows()
save_fig(ggplot(pn, aes(correct)) + geom_histogram(binwidth = 1, fill = "grey75") +
           geom_vline(data = distinct(pn, k, obs), aes(xintercept = obs), colour = colours$above, linewidth = 1) + facet_wrap(~k) +
           scale_x_continuous(breaks = 0:length(is_rel)) +
           labs(x = sprintf("patients predicted correctly (of %d), markers picked without the predicted patient", length(is_rel)), y = "relabellings",
                title = "Can V1 dISF predict relapse? Observed (violet) vs all random labellings",
                subtitle = paste(map_chr(nested, \(x) sprintf("top %d: %d correct, P = %.3f", x$k, x$obs, x$p)), collapse = "; ")) + th,
         "V1_relapse_prediction_test", width = 9, height = 5)
if (nrow(gsea$V1_lesion)) save_fig(gsea_plot(gsea$V1_lesion, "GSEA, dISF lesion site at V1, relapse vs non-relapse", colours, fdr), "V1_relapse_GSEA")
msg("V1 relapse screen: %s", file.path(cfg$paths$output, out))
