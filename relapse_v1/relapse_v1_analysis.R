# V1 relapse screen (run with run_relapse_v1.R): which dISF (and serum) proteins are already higher at baseline (V1, before the
# topical steroids at V2) in the AD patients who later relapse? MicroAD, AD patients with a relapse / non-relapse label
# (dropouts excluded). Groups are tiny (4 vs 6), so the list is a ranked screen with an exact relabelling check.
#   V1 lesion site (primary), V1 non-lesional skin, V1 lesion minus non-lesional (same patient),
#   V2 cleared lesion site (still higher after treatment?), V1 serum (also visible in blood?)
# In:  output of steps 01-02 (data/npx_clean.rds); Human Protein Atlas (drug targets, location); MSigDB / KEGG gene sets
# Out: output/relapse_V1/ - V1_relapse_targets.xlsx, V1_relapse_GSEA.xlsx, V1_relapse_volcano.pdf/.png,
#      V1_relapse_targets_course.pdf/.png, V1_relapse_chance_check.pdf/.png, V1_relapse_GSEA.pdf/.png,
#      V1_relapse_volcano_serum.pdf/.png, V1_relapse_dISF_vs_serum.pdf/.png, V1_relapse_markers_dISF_and_serum_course.pdf/.png,
#      V1_relapse_executive_summary.pdf (answer, predictors in dISF / serum / both, gene sets, methods)

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
samples_where <- function(...) ad |> filter(...) |> distinct(SampleID, SubjectID, relapse) |>
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

# ---- the screen (dISF V1 lesion site, and the same for V1 serum): exact relabelling, leave-one-out, AUC, prediction ------
pred_k <- unlist(rv$predict_k %||% c(1, 5))
time_of <- ad |> distinct(SubjectID, relapse, relapse_visit) |>
  mutate(t = if_else(relapse == "relapse", relapse_visit, max(ad$visit_num, na.rm = TRUE) + 1))
screen_matrix <- function(r, smp, label) {
  M <- patient_matrix(clean, smp, r$OlinkID)
  M <- M[, sort(colnames(M)), drop = FALSE]                    # same patient order in every matrix (joint relabelling)
  is_rel <- smp$grp[match(colnames(M), smp$SubjectID)] == "relapse"
  rc <- relabel_check(M, is_rel)
  lab_of <- \(j) seq_along(is_rel) %in% rc$combs[, j]
  nested <- map(pred_k, \(k) {
    obs <- nested_loo_correct(M, is_rel, k)
    null <- map_int(seq_len(ncol(rc$combs)), \(j) nested_loo_correct(M, lab_of(j), k))
    list(k = k, obs = obs, null = null, p = mean(null >= obs))
  })
  walk(nested, \(x) msg("%s - prediction with the top %d protein(s), chosen without the predicted patient: %d of %d correct; random labels: median %.0f; P = %.3f",
                        label, x$k, x$obs, length(is_rel), median(x$null), x$p))
  ps_obs <- perfect_separators(M, matrix(which(is_rel))); ps_null <- perfect_separators(M, rc$combs)
  tvec <- time_of$t[match(colnames(M), time_of$SubjectID)]
  rho <- apply(M, 1, \(v) suppressWarnings(cor(v, tvec, method = "spearman", use = "pairwise.complete.obs")))
  ptime <- apply(M, 1, \(v) tryCatch(suppressWarnings(cor.test(v, tvec, method = "spearman", exact = FALSE)$p.value), error = \(e) NA_real_))
  per_marker <- map(rownames(M), \(oid) loo_predict_marker(M[oid, ], is_rel) |> mutate(OlinkID = oid)) |> bind_rows() |>
    mutate(rho_time_to_relapse = round(rho[OlinkID], 2), p_time_to_relapse = signif(ptime[OlinkID], 3))
  auc <- auc_rows(M, is_rel)
  per <- r |> left_join(rc$per_protein |> select(OlinkID, p_exact_up, p_exact_down, p_exact_two_sided), by = "OlinkID") |>
    left_join(leave_one_out(M, is_rel), by = "OlinkID") |> mutate(AUC = round(auc[OlinkID], 2)) |> left_join(per_marker, by = "OlinkID") |>
    mutate(detected_in_relapse = `pct_above_LOD_relapse` >= 100 * min_f,
           detected_in_nonrelapse = `pct_above_LOD_non-relapse` >= 100 * min_f,
           candidate_up = logFC > 0 & p < 0.05 & p_exact_up < 0.05 & detected_in_relapse,
           candidate_down = logFC < 0 & p < 0.05 & p_exact_down < 0.05 & detected_in_nonrelapse)
  chance <- tibble(item = c("patients", "relabellings (all ways to choose the relapse group)", "proteins higher in relapse with p < 0.05 (observed)",
                            "the same under random labels: median (90th percentile)", "share of random labellings with at least as many",
                            "estimated chance hits among the observed"),
                   value = c(sprintf("%d relapse, %d non-relapse", sum(is_rel), sum(!is_rel)), as.character(rc$n_labellings), as.character(rc$observed_hits),
                             sprintf("%.0f (%.0f)", median(rc$hits), quantile(rc$hits, 0.9)), sprintf("%.2f", mean(rc$hits >= rc$observed_hits)),
                             sprintf("about %.0f of %d", min(median(rc$hits), rc$observed_hits), rc$observed_hits)))
  predict_tbl <- tibble(item = c("proteins separating relapse and non-relapse perfectly (AUC 1 or 0)", "the same under random labels: median (90th percentile)",
                                 map_chr(nested, \(x) sprintf("patients predicted correctly with the top %d protein(s), picked on the other patients (leave one out)", x$k)),
                                 map_chr(nested, \(x) sprintf("the same under random labels (top %d): median; P (at least as many correct)", x$k))),
                        value = c(as.character(ps_obs), sprintf("%.0f (%.0f)", median(ps_null), quantile(ps_null, 0.9)),
                                  map_chr(nested, \(x) sprintf("%d of %d", x$obs, length(is_rel))),
                                  map_chr(nested, \(x) sprintf("%.0f; P = %.3f", median(x$null), x$p))))
  msg("%s - chance check: %d proteins higher in relapse (p < 0.05); random labels give a median of %.0f; P = %.2f",
      label, rc$observed_hits, median(rc$hits), mean(rc$hits >= rc$observed_hits))
  list(M = M, is_rel = is_rel, rc = rc, nested = nested, per = per, chance = chance, predict_tbl = predict_tbl,
       predicts = any(map_lgl(nested, \(x) x$p < 0.05)))
}
sd <- screen_matrix(res$V1_lesion, s$V1_lesion, "dISF")
ss <- screen_matrix(res$serum_V1, s$serum_V1, "serum")
v1 <- res$V1_lesion; M <- sd$M; is_rel <- sd$is_rel; rc <- sd$rc; nested <- sd$nested
chance <- sd$chance; predict_tbl <- sd$predict_tbl

# ---- annotation: drug targets and location (Human Protein Atlas) --------------------------------------------------------
hpa <- read_hpa(cfg$paths$hpa)
anno <- if (is.null(hpa)) tibble(OlinkID = character()) else {
  a <- assay_protein_classes(bind_rows(v1, res$serum_V1) |> distinct(OlinkID, Assay, UniProt), hpa,
                             classes = list(approved_drug_target = "FDA approved drug targets", potential_drug_target = "Potential drug targets",
                                            secreted = "Predicted secreted proteins", membrane = "Predicted membrane proteins"))$classes
  a |> mutate(v = TRUE) |> pivot_wider(names_from = class, values_from = v, values_fill = FALSE)
}
for (cl in c("approved_drug_target", "potential_drug_target", "secreted", "membrane")) if (!cl %in% names(anno)) anno[[cl]] <- logical(nrow(anno))

# ---- ranked lists -------------------------------------------------------------------------------------------------------
side <- \(nm, pre) res[[nm]] |> select(OlinkID, logFC, p) |> rename_with(\(x) paste0(x, "_", pre), c(logFC, p))
flags <- \(x) x |> left_join(anno, by = "OlinkID") |> mutate(across(c(approved_drug_target, potential_drug_target, secreted, membrane), \(v) coalesce(v, FALSE)))
scr <- sd$per |>
  left_join(side("V1_nonlesional", "V1_nonlesional"), by = "OlinkID") |>
  left_join(side("V1_lesion_minus_nonlesional", "V1_LminusNL"), by = "OlinkID") |>
  left_join(side("V2_cleared_lesion", "V2_cleared"), by = "OlinkID") |>
  left_join(side("serum_V1", "serum_V1"), by = "OlinkID") |> flags() |>
  mutate(persists_V2_up = coalesce(logFC_V2_cleared > 0 & p_V2_cleared < 0.05, FALSE),
         persists_V2_down = coalesce(logFC_V2_cleared < 0 & p_V2_cleared < 0.05, FALSE),
         also_serum_up = coalesce(logFC_serum_V1 > 0 & p_serum_V1 < 0.05, FALSE),
         also_serum_down = coalesce(logFC_serum_V1 < 0 & p_serum_V1 < 0.05, FALSE))
rank_up <- \(x) x |> filter(candidate_up) |> arrange(desc(loo_correct), desc(AUC), desc(loo_always_higher), desc(logFC)) |> mutate(rank = row_number())
rank_down <- \(x) x |> filter(candidate_down) |> arrange(desc(loo_correct), AUC, logFC) |> mutate(rank = row_number())
cand_up <- scr |> filter(candidate_up) |> arrange(desc(loo_correct), desc(AUC), desc(persists_V2_up), desc(loo_always_higher), desc(logFC))
cand_up <- cand_up |> mutate(rank = row_number())
cand_down <- scr |> filter(candidate_down) |> arrange(desc(loo_correct), AUC, desc(persists_V2_down), logFC) |>
  mutate(rank = row_number())
tcols <- c("rank", "Assay", "OlinkID", "UniProt", "loo_correct", "loo_n", "loo_sensitivity", "loo_specificity", "cutoff_NPX", "AUC",
           "rho_time_to_relapse", "p_time_to_relapse", "logFC", "p", "p_exact_up", "p_exact_down", "loo_always_higher", "loo_max_p",
           "persists_V2_up", "persists_V2_down", "logFC_V2_cleared", "p_V2_cleared", "logFC_V1_nonlesional", "p_V1_nonlesional",
           "logFC_V1_LminusNL", "p_V1_LminusNL", "also_serum_up", "also_serum_down", "logFC_serum_V1", "p_serum_V1",
           "also_dISF_up", "also_dISF_down", "logFC_dISF_V1", "p_dISF_V1",
           "approved_drug_target", "potential_drug_target", "secreted", "membrane",
           "mean_NPX_relapse", "mean_NPX_non-relapse", "pct_above_LOD_relapse", "pct_above_LOD_non-relapse", "FDR")
msg("dISF targets: %d higher, %d lower in relapsers at V1 (lesion site)", nrow(cand_up), nrow(cand_down))

# serum at V1: the same screen, with the dISF lesion site result next to it
scr_s <- ss$per |> left_join(side("V1_lesion", "dISF_V1"), by = "OlinkID") |> flags() |>
  mutate(also_dISF_up = coalesce(logFC_dISF_V1 > 0 & p_dISF_V1 < 0.05, FALSE), also_dISF_down = coalesce(logFC_dISF_V1 < 0 & p_dISF_V1 < 0.05, FALSE))
ser_up <- rank_up(scr_s)
ser_down <- rank_down(scr_s)
msg("Serum targets: %d higher, %d lower in relapsers at V1", nrow(ser_up), nrow(ser_down))

# ---- dISF and serum together: markers in both matrices ------------------------------------------------------------------
pick_cols <- \(x, pre) x |> transmute(OlinkID, Assay, UniProt, up = candidate_up, down = candidate_down, logFC, p, p_exact_up, p_exact_down,
                                      correct = loo_correct, n = loo_n, AUC, cutoff_NPX) |>
  rename_with(\(v) paste0(pre, "_", v), -c(OlinkID, Assay, UniProt))
ovl <- inner_join(pick_cols(scr, "dISF"), pick_cols(scr_s, "serum") |> select(-Assay, -UniProt), by = "OlinkID") |>
  left_join(scr |> select(OlinkID, persists_V2_up, persists_V2_down), by = "OlinkID") |> flags() |>
  mutate(category = case_when(
    dISF_up & serum_up ~ "higher in relapse: marker in dISF and serum",
    dISF_up & serum_logFC > 0 & serum_p < 0.05 ~ "higher in relapse: marker in dISF, supported in serum",
    serum_up & dISF_logFC > 0 & dISF_p < 0.05 ~ "higher in relapse: marker in serum, supported in dISF",
    dISF_down & serum_down ~ "lower in relapse: marker in dISF and serum",
    dISF_down & serum_logFC < 0 & serum_p < 0.05 ~ "lower in relapse: marker in dISF, supported in serum",
    serum_down & dISF_logFC < 0 & dISF_p < 0.05 ~ "lower in relapse: marker in serum, supported in dISF",
    TRUE ~ NA_character_)) |>
  filter(!is.na(category)) |>
  arrange(factor(category, c("higher in relapse: marker in dISF and serum", "higher in relapse: marker in dISF, supported in serum",
                             "higher in relapse: marker in serum, supported in dISF", "lower in relapse: marker in dISF and serum",
                             "lower in relapse: marker in dISF, supported in serum", "lower in relapse: marker in serum, supported in dISF")),
          desc(dISF_correct + serum_correct), desc(dISF_AUC + serum_AUC))
both_up <- ovl |> filter(category == "higher in relapse: marker in dISF and serum")
# how many proteins would be "higher with p < 0.05 in both matrices" by chance? the same relabelling for both (same patients)
common <- intersect(rownames(sd$rc$up_hits), rownames(ss$rc$up_hits))
joint_ok <- identical(colnames(sd$M), colnames(ss$M))
ovl_obs <- sum(scr$OlinkID %in% common & scr$p < 0.05 & scr$logFC > 0 &
               scr$OlinkID %in% scr_s$OlinkID[scr_s$p < 0.05 & scr_s$logFC > 0])
ovl_null <- if (joint_ok) colSums(sd$rc$up_hits[common, , drop = FALSE] & ss$rc$up_hits[common, , drop = FALSE]) else NA_integer_
overlap_chance <- tibble(item = c("proteins measured in both matrices", "higher in relapse with p < 0.05 in both dISF and serum (observed)",
                                  "the same under random labels (same labels for both matrices): median (90th percentile)",
                                  "share of random labellings with at least as many", "markers in both matrices (all criteria in both)"),
                         value = c(as.character(length(common)), as.character(ovl_obs),
                                   if (joint_ok) sprintf("%.0f (%.0f)", median(ovl_null), quantile(ovl_null, 0.9)) else "not computed (different patients)",
                                   if (joint_ok) sprintf("%.2f", mean(ovl_null >= ovl_obs)) else "", as.character(nrow(both_up))))
msg("dISF and serum: %d proteins higher in relapse (p < 0.05) in both; random labels give a median of %s; %d markers in both", ovl_obs,
    if (joint_ok) sprintf("%.0f", median(ovl_null)) else "-", nrow(both_up))

answer <- if (sd$predicts) "V1 dISF predicts relapse better than chance in these patients (sheet can_V1_predict); the ranked candidates are in targets_up / targets_down" else
  "With these patients, predicting relapse from V1 dISF is not better than chance once the search over all proteins is accounted for (sheet can_V1_predict); the targets are candidates to test in new patients, not proven predictors"
answer_serum <- if (ss$predicts) "V1 serum predicts relapse better than chance in these patients (sheet serum_can_V1_predict)" else
  "V1 serum does not predict relapse better than chance in these patients once the search over all proteins is accounted for (sheet serum_can_V1_predict)"
answer_both <- if (nrow(both_up)) sprintf("%d proteins are markers (higher in relapsers at V1) in both dISF and serum: %s (sheet dISF_and_serum)",
                                          nrow(both_up), paste(head(both_up$Assay, 10), collapse = ", ")) else
  "No protein is a marker in both dISF and serum; proteins supported in the other matrix are in sheet dISF_and_serum"
readme <- tibble(item = c("question", "answer", "answer serum", "answer dISF and serum", "patients", "relapse visits", "comparisons", "target (higher)", "ranking", "prediction columns", "flags", "chance check", "caution"),
  note = c("Which dISF markers at V1 (baseline, before the topical steroids) predict a later relapse - and which also in serum?", answer, answer_serum, answer_both,
           paste(n_txt, collapse = "; "), rel_visits,
           "limma per protein, relapse vs non-relapse: V1 lesion site (primary), V1 non-lesional skin, V1 lesion minus non-lesional (same patient), V2 cleared lesion site (before any relapse), V1 serum. log2FC = relapse - non-relapse",
           sprintf("V1 lesion site: log2FC > 0, p < 0.05, exact relabelling p < 0.05 (all %d ways to choose the relapse group), above LOD in >= %.0f%% of the relapsers", rc$n_labellings, 100 * min_f),
           "patients predicted correctly when each is left out (loo_correct), then AUC, then still higher at V2 after treatment, then higher in every leave-one-patient-out refit, then fold change",
           "loo_correct / loo_sensitivity / loo_specificity: each patient classified with a cut-off set on the other patients (midpoint of the group means); cutoff_NPX: that cut-off on all patients (V1 lesion site, NPX); rho_time_to_relapse: Spearman with the relapse visit (non-relapse = after the last visit), < 0 = higher value, earlier relapse",
           "persists_V2: still higher in the cleared lesion site at V2; also_serum: also higher in V1 serum; approved / potential drug target, secreted, membrane: Human Protein Atlas",
           "sheet chance_check: how many proteins random labels would give; if the observed number is not above that, the list is mostly chance",
           "4 vs 6 patients: a screen to pick candidates, not a confirmation; the smallest possible exact p is 1/210"))
sig3 <- \(x) x |> mutate(across(where(is.double), \(v) signif(v, 3)))
writexl::write_xlsx(list(README = readme, can_V1_predict = predict_tbl, chance_check = chance,
                         targets_up = cand_up |> select(any_of(tcols)) |> sig3(),
                         targets_down = cand_down |> select(any_of(tcols)) |> sig3(),
                         dISF_and_serum = ovl |> relocate(category, .after = UniProt) |> sig3(), overlap_chance = overlap_chance,
                         serum_can_V1_predict = ss$predict_tbl, serum_chance_check = ss$chance,
                         serum_targets_up = ser_up |> select(any_of(tcols)) |> sig3(), serum_targets_down = ser_down |> select(any_of(tcols)) |> sig3(),
                         V1_lesion = res$V1_lesion, V1_nonlesional = res$V1_nonlesional, V1_lesion_minus_nonlesional = res$V1_lesion_minus_nonlesional,
                         V2_cleared_lesion = res$V2_cleared_lesion, serum_V1 = res$serum_V1),
                    out_path(cfg, out, "V1_relapse_targets.xlsx"))

# ---- GSEA ---------------------------------------------------------------------------------------------------------------
sets <- v3_gene_sets(cfg)
gsea <- map(res[c("V1_lesion", "V1_lesion_minus_nonlesional", "V2_cleared_lesion", "serum_V1")],
            \(r) run_gsea(r, sets, cfg$enrichment$min_size %||% 10, cfg$enrichment$max_size %||% 500))
gs_both <- if (nrow(gsea$V1_lesion) && nrow(gsea$serum_V1))
  inner_join(gsea$V1_lesion |> select(collection, gene_set, NES_dISF = NES, p_dISF = p, FDR_dISF = FDR, leading_edge_dISF = leading_edge),
             gsea$serum_V1 |> select(collection, gene_set, NES_serum = NES, p_serum = p, FDR_serum = FDR, leading_edge_serum = leading_edge),
             by = c("collection", "gene_set")) |>
  mutate(category = case_when(sign(NES_dISF) != sign(NES_serum) ~ NA_character_,
                              FDR_dISF < fdr & FDR_serum < fdr ~ "FDR < 0.05 in both",
                              FDR_dISF < fdr & p_serum < 0.05 ~ "FDR < 0.05 in dISF, p < 0.05 in serum",
                              FDR_serum < fdr & p_dISF < 0.05 ~ "FDR < 0.05 in serum, p < 0.05 in dISF",
                              p_dISF < 0.05 & p_serum < 0.05 ~ "p < 0.05 in both"),
         direction = if_else(NES_dISF > 0, "higher in relapsers", "lower in relapsers")) |>
  filter(!is.na(category)) |> arrange(direction, factor(category, c("FDR < 0.05 in both", "FDR < 0.05 in dISF, p < 0.05 in serum",
                                                                    "FDR < 0.05 in serum, p < 0.05 in dISF", "p < 0.05 in both")), p_dISF) |>
  relocate(direction, category, .after = gene_set) else tibble()
msg("Gene sets changed in the same direction in dISF and serum (V1): %d", nrow(gs_both))
writexl::write_xlsx(c(list(README = tibble(note = c("GSEA (fgsea) on the limma t statistics, relapse vs non-relapse; NES > 0 = higher in relapsers. MSigDB collections of config enrichment$collections and KEGG",
                                                    "gene_sets_dISF_and_serum: gene sets changed in the same direction at V1 in the dISF lesion site and in serum"))),
                      gsea, list(gene_sets_dISF_and_serum = gs_both)), out_path(cfg, out, "V1_relapse_GSEA.xlsx"))

# ---- figures --------------------------------------------------------------------------------------------------------------
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
save_fig <- function(p, file, width = 11.69, height = 8.27) {
  for (type in c("pdf", "png")) {
    open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
    print(set_family(p, fi$family))
    close_device(fi)
  }
  invisible(p)
}
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot")
volcano_screen <- function(x, what, scx, file) {
  rcx <- scx$rc
  v <- x |> mutate(cl = case_when(candidate_up ~ "higher in relapse (marker)", candidate_down ~ "lower in relapse (marker)", TRUE ~ "other"))
  lab <- v |> arrange(desc(candidate_up | candidate_down), p) |> slice_head(n = n_lab)
  save_fig(ggplot(v, aes(logFC, -log10(p))) +
             geom_point(aes(colour = cl), size = 1.1) +
             ggrepel::geom_text_repel(data = lab, aes(label = Assay), size = 2.6, max.overlaps = Inf, min.segment.length = 0, segment.colour = "grey60", segment.size = 0.3) +
             scale_colour_manual(values = c(`higher in relapse (marker)` = colours$above, `lower in relapse (marker)` = colours$below, other = "grey75"), name = NULL) +
             labs(x = "log2 fold change (relapse - non-relapse)", y = "-log10 p",
                  title = sprintf("%s at V1: %d later relapse vs %d non-relapse - %d proteins higher, %d lower (markers)",
                                  what, sum(scx$is_rel), sum(!scx$is_rel), sum(v$candidate_up), sum(v$candidate_down)),
                  subtitle = sprintf("Marker: p < 0.05 and exact relabelling p < 0.05. %d proteins higher with p < 0.05; random labels give a median of %.0f. Labelled: %d proteins (markers first, then lowest p)",
                                     rcx$observed_hits, median(rcx$hits), n_lab)) + th, file)
}
p_volcano <- volcano_screen(scr, "dISF lesion site", sd, "V1_relapse_volcano")
p_volcano_s <- volcano_screen(scr_s, "Serum", ss, "V1_relapse_volcano_serum")

show <- (if (nrow(cand_up)) cand_up else scr |> filter(logFC > 0) |> arrange(p)) |> slice_head(n = n_course)
course <- clean |> filter(cohort == "MicroAD", matrix == "ISF", site == "L", SubjectID %in% s$V1_lesion$SubjectID, OlinkID %in% show$OlinkID) |>
  select(SubjectID, OlinkID, Assay, visit_num, value, relapse, relapse_visit) |>
  mutate(Assay = factor(Assay, unique(show$Assay)), at_relapse = !is.na(relapse_visit) & visit_num == relapse_visit)
p_course <- save_fig(ggplot(course, aes(visit_num, value, group = SubjectID, colour = relapse)) +
           geom_line(alpha = 0.8) + geom_point(size = 1.2) +
           geom_point(data = \(d) filter(d, at_relapse), shape = 4, size = 3, stroke = 1.1, colour = "black") +
           scale_colour_manual(values = c(relapse = colours$above, `non-relapse` = "grey60"), name = NULL) +
           scale_x_continuous(breaks = 1:6, labels = paste0("V", 1:6)) + facet_wrap(~Assay, scales = "free_y", ncol = 4) +
           labs(x = NULL, y = "NPX, dISF lesion site",
                title = if (nrow(cand_up)) sprintf("Top %d targets: dISF lesion site per patient over the visits", nrow(show)) else
                  sprintf("No target passed - the %d proteins with the lowest p (higher in relapse)", nrow(show)),
                subtitle = "violet: patients who relapse (x = relapse visit); grey: non-relapse. V1 baseline, V2 after topical steroids") + th,
         "V1_relapse_targets_course")

p_chance <- save_fig(ggplot(tibble(hits = rc$hits), aes(hits)) + geom_histogram(binwidth = 1, fill = "grey75") +
           geom_vline(xintercept = rc$observed_hits, colour = colours$above, linewidth = 1) +
           labs(x = "proteins higher in the 'relapse' group with p < 0.05", y = "relabellings",
                title = sprintf("Chance check: observed %d (violet) vs all %d ways to choose %d of the %d patients as relapse",
                                rc$observed_hits, rc$n_labellings, sum(is_rel), length(is_rel)),
                subtitle = sprintf("%.0f%% of the random labellings give at least as many proteins", 100 * mean(rc$hits >= rc$observed_hits))) + th,
         "V1_relapse_chance_check", width = 9, height = 5)
pn <- bind_rows(map(sd$nested, \(x) tibble(m = "dISF lesion site", k = sprintf("top %d protein(s)", x$k), correct = x$null, obs = x$obs)),
                map(ss$nested, \(x) tibble(m = "serum", k = sprintf("top %d protein(s)", x$k), correct = x$null, obs = x$obs)))
p_pred <- save_fig(ggplot(pn, aes(correct)) + geom_histogram(binwidth = 1, fill = "grey75") +
           geom_vline(data = distinct(pn, m, k, obs), aes(xintercept = obs), colour = colours$above, linewidth = 1) + facet_grid(m ~ k) +
           scale_x_continuous(breaks = 0:length(is_rel)) +
           labs(x = sprintf("patients predicted correctly (of %d), markers picked without the predicted patient", length(is_rel)), y = "relabellings",
                title = "Can V1 dISF / serum predict relapse? Observed (violet) vs all random labellings",
                subtitle = paste(c(map_chr(sd$nested, \(x) sprintf("dISF top %d: %d correct, P = %.3f", x$k, x$obs, x$p)),
                                   map_chr(ss$nested, \(x) sprintf("serum top %d: %d correct, P = %.3f", x$k, x$obs, x$p))), collapse = "; ")) + th,
         "V1_relapse_prediction_test", width = 10, height = 6.5)

# dISF vs serum: fold changes of all proteins measured in both, markers coloured, the overlap labelled
sc <- inner_join(scr |> select(OlinkID, Assay, logFC_d = logFC, up_d = candidate_up, down_d = candidate_down),
                 scr_s |> select(OlinkID, logFC_s = logFC, up_s = candidate_up, down_s = candidate_down), by = "OlinkID") |>
  mutate(cl = case_when((up_d & up_s) | (down_d & down_s) ~ "marker in both", up_d | down_d ~ "marker in dISF", up_s | down_s ~ "marker in serum",
                        TRUE ~ "other"))
p_scatter <- save_fig(ggplot(sc, aes(logFC_s, logFC_d)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
           geom_point(data = \(d) filter(d, cl == "other"), colour = "grey80", size = 0.9) +
           geom_point(data = \(d) filter(d, cl != "other"), aes(colour = cl), size = 1.8) +
           ggrepel::geom_text_repel(data = \(d) filter(d, cl == "marker in both"), aes(label = Assay), size = 2.8, max.overlaps = Inf, min.segment.length = 0) +
           scale_colour_manual(values = c(`marker in both` = "black", `marker in dISF` = colours$above, `marker in serum` = colours$below), name = NULL) +
           labs(x = "log2FC serum V1 (relapse - non-relapse)", y = "log2FC dISF lesion site V1 (relapse - non-relapse)",
                title = sprintf("dISF vs serum at V1: %d markers in both matrices (labelled)", sum(sc$cl == "marker in both")),
                subtitle = sprintf("%d proteins measured in both. Higher with p < 0.05 in both: %d; random labels give a median of %s", nrow(sc), ovl_obs,
                                   if (joint_ok) sprintf("%.0f", median(ovl_null)) else "-")) + th,
         "V1_relapse_dISF_vs_serum")
show_b <- (if (nrow(both_up)) both_up else ovl |> filter(str_detect(category, "^higher"))) |> slice_head(n = 8)
p_course_b <- if (nrow(show_b)) {
  cb <- clean |> filter(cohort == "MicroAD", SubjectID %in% s$V1_lesion$SubjectID, OlinkID %in% show_b$OlinkID,
                        (matrix == "ISF" & site == "L") | matrix == "Serum") |>
    transmute(SubjectID, Assay = factor(Assay, unique(show_b$Assay)), m = if_else(matrix == "ISF", "dISF lesion site", "serum"),
              visit_num, value, relapse, at_relapse = !is.na(relapse_visit) & visit_num == relapse_visit)
  save_fig(ggplot(cb, aes(visit_num, value, group = SubjectID, colour = relapse)) + geom_line(alpha = 0.8) + geom_point(size = 1) +
             geom_point(data = \(d) filter(d, at_relapse), shape = 4, size = 2.6, stroke = 1, colour = "black") +
             scale_colour_manual(values = c(relapse = colours$above, `non-relapse` = "grey60"), name = NULL) +
             scale_x_continuous(breaks = 1:6, labels = paste0("V", 1:6)) + facet_grid(m ~ Assay, scales = "free_y") +
             labs(x = NULL, y = "NPX", title = if (nrow(both_up)) "Markers in dISF and serum: per patient over the visits" else
                    "No marker in both matrices - markers supported in the other matrix",
                  subtitle = "violet: patients who relapse (x = relapse visit); grey: non-relapse") + th,
           "V1_relapse_markers_dISF_and_serum_course")
} else NULL
p_gsea <- if (nrow(gsea$V1_lesion)) save_fig(gsea_plot(gsea$V1_lesion, "GSEA, dISF lesion site at V1, relapse vs non-relapse", colours, fdr), "V1_relapse_GSEA") else NULL

# ---- executive summary PDF ------------------------------------------------------------------------------------------------
source("R/report.R")
fmt_p <- \(x) ifelse(is.na(x), "", formatC(x, digits = 2, format = "g"))
yn <- \(x) ifelse(x %in% TRUE, "yes", "")
pretty_gs <- \(x) if_else(x == toupper(x), str_to_sentence(x), x)        # MSigDB names are upper case
gs_up <- \(g) if (is.null(g) || !nrow(g)) tibble() else g |> filter(NES > 0) |> arrange(p) |> mutate(gene_set = pretty_gs(gene_set))
g1 <- gs_up(gsea$V1_lesion); g2 <- gs_up(gsea$V2_cleared_lesion)
top_t <- head(cand_up, 5)
pred_line <- map_chr(nested, \(x) sprintf("top %d protein(s): %d of %d patients correct (random labels: median %.0f; P = %.3f)",
                                          x$k, x$obs, length(is_rel), median(x$null), x$p))
drug_n <- sum(cand_up$approved_drug_target | cand_up$potential_drug_target)
pdf_open(out_path(cfg, out, "V1_relapse_executive_summary.pdf"))
section("Summary", page_text("Markers at V1 that predict relapse (dISF and serum) - executive summary", c(
  "## Question",
  "Which proteins at V1 (baseline, before the topical steroids at V2) predict that an AD patient relapses later - in dISF, in serum, and in both?",
  "## Patients",
  sprintf("MicroAD, dISF lesion site at V1: %d patients who relapse later (relapse at %s) vs %d who do not (dropouts excluded).",
          sum(is_rel), rel_visits, sum(!is_rel)),
  "## Answer",
  answer,
  paste("Can V1 dISF predict relapse? Each patient was predicted with markers picked and cut-offs set on the other patients:",
        paste(pred_line, collapse = "; "), "."),
  sprintf("Chance check: %d proteins are higher in relapsers with p < 0.05; random labels give a median of %.0f, so about %.0f of them are expected by chance.",
          rc$observed_hits, median(rc$hits), min(median(rc$hits), rc$observed_hits)),
  "## Markers (higher in relapsers at V1)",
  if (nrow(cand_up)) sprintf("%d candidate markers (%d lower in relapsers). Top: %s.", nrow(cand_up), nrow(cand_down),
                             paste(sprintf("%s (%d/%d correct, AUC %.2f)", top_t$Assay, top_t$loo_correct, top_t$loo_n, top_t$AUC), collapse = ", "))
  else "No protein passed the marker criteria.",
  sprintf("Still higher after treatment (V2, cleared skin): %d; also higher in V1 serum: %d; approved or potential drug target: %d; secreted or membrane: %d.",
          sum(cand_up$persists_V2_up), sum(cand_up$also_serum_up), drug_n, sum(cand_up$secreted | cand_up$membrane)),
  "## Serum at V1",
  answer_serum,
  paste("Serum prediction:", paste(map_chr(ss$nested, \(x) sprintf("top %d protein(s): %d of %d correct (random labels: median %.0f; P = %.3f)",
                                                                   x$k, x$obs, length(ss$is_rel), median(x$null), x$p)), collapse = "; "), "."),
  if (nrow(ser_up)) sprintf("%d serum markers higher in relapsers. Top: %s.", nrow(ser_up),
                            paste(sprintf("%s (%d/%d correct, AUC %.2f)", head(ser_up$Assay, 5), head(ser_up$loo_correct, 5), head(ser_up$loo_n, 5), head(ser_up$AUC, 5)), collapse = ", "))
  else "No serum protein passed the marker criteria.",
  "## dISF and serum together",
  answer_both,
  sprintf("Higher in relapsers with p < 0.05 in both matrices: %d proteins; random labels (the same for both matrices) give a median of %s.", ovl_obs,
          if (joint_ok) sprintf("%.0f", median(ovl_null)) else "-"),
  "## Gene sets (higher in relapsers)",
  sprintf("V1 lesion site: %d gene sets with FDR < %.2f%s.", sum(g1$FDR < fdr), fdr,
          if (nrow(g1)) paste0("; top: ", paste(head(g1$gene_set, 5), collapse = "; ")) else ""),
  sprintf("V2 cleared lesion site: %d gene sets with FDR < %.2f%s.", sum(g2$FDR < fdr), fdr,
          if (nrow(g2)) paste0("; top: ", paste(head(g2$gene_set, 5), collapse = "; ")) else ""),
  { gb <- if (nrow(gs_both)) gs_both |> filter(direction == "higher in relapsers") else tibble()
    sprintf("Higher in relapsers in both dISF and serum (V1): %d gene sets%s.", nrow(gb),
            if (nrow(gb)) paste0("; top: ", paste(head(pretty_gs(gb$gene_set), 5), collapse = "; ")) else "") },
  "## Keep in mind",
  sprintf("%d vs %d patients: these are candidates to test in new patients, not confirmed predictors. The smallest possible exact p is 1/%d.",
          sum(is_rel), sum(!is_rel), rc$n_labellings),
  "Full tables: V1_relapse_targets.xlsx (dISF: targets_up / targets_down; serum: serum_targets_up / serum_targets_down; both: dISF_and_serum; prediction and chance sheets; all comparisons) and V1_relapse_GSEA.xlsx."),
  subtitle = sprintf("O-MicroAD, Olink Explore HT, dISF and serum. Generated %s.", format(Sys.Date()))))
section("Prediction test", { page_plot(p_pred); page_plot(p_chance) })
section("Markers", page_table("Predictors: markers higher in relapsers at V1 (ranked)",
  cand_up |> head(30) |> transmute(rank, marker = Assay, correct = sprintf("%d/%d", loo_correct, loo_n), sens = round(loo_sensitivity, 2),
                                   spec = round(loo_specificity, 2), cutoff_NPX = round(cutoff_NPX, 2), AUC, log2FC = round(logFC, 2), p = fmt_p(p),
                                   exact_p = fmt_p(p_exact_up), V2 = yn(persists_V2_up), serum = yn(also_serum_up),
                                   drug = if_else(approved_drug_target, "approved", if_else(potential_drug_target, "potential", "")),
                                   location = if_else(secreted, "secreted", if_else(membrane, "membrane", ""))),
  subtitle = "Ranked by patients predicted correctly, then AUC. correct: each patient predicted with a cut-off set on the other patients (leave one out)",
  note = "cutoff_NPX: V1 lesion site. V2: still higher at V2 (cleared skin). serum: also higher in V1 serum. drug / location: Human Protein Atlas", rows_per_page = 25))
section("Course", page_plot(p_course))
section("dISF and serum", {
  page_table("Markers in dISF and serum (relapse vs non-relapse at V1)",
    ovl |> head(30) |> transmute(marker = Assay, category = str_remove(category, "^(higher|lower) in relapse: "), dISF = sprintf("%d/%d, AUC %.2f", dISF_correct, dISF_n, dISF_AUC),
                                 serum = sprintf("%d/%d, AUC %.2f", serum_correct, serum_n, serum_AUC), log2FC_dISF = round(dISF_logFC, 2), log2FC_serum = round(serum_logFC, 2),
                                 cutoff_dISF = round(dISF_cutoff_NPX, 2), cutoff_serum = round(serum_cutoff_NPX, 2), V2 = yn(persists_V2_up),
                                 drug = if_else(approved_drug_target, "approved", if_else(potential_drug_target, "potential", ""))),
    subtitle = "dISF / serum: patients predicted correctly (leave one out) and AUC in each matrix; cut-offs in NPX at V1",
    note = "marker = p < 0.05 and exact p < 0.05 in that matrix; supported = p < 0.05 in the same direction in the other matrix", rows_per_page = 25)
  page_plot(p_scatter); if (!is.null(p_course_b)) page_plot(p_course_b) })
section("Serum markers", page_table("Serum: markers higher in relapsers at V1 (ranked)",
  ser_up |> head(25) |> transmute(rank, marker = Assay, correct = sprintf("%d/%d", loo_correct, loo_n), sens = round(loo_sensitivity, 2), spec = round(loo_specificity, 2),
                                  cutoff_NPX = round(cutoff_NPX, 2), AUC, log2FC = round(logFC, 2), p = fmt_p(p), exact_p = fmt_p(p_exact_up), dISF = yn(also_dISF_up)),
  subtitle = "dISF: also higher in the V1 dISF lesion site (p < 0.05)", rows_per_page = 25))
section("Serum volcano", page_plot(p_volcano_s))
section("Volcano", page_plot(p_volcano))
section("Gene sets V1", page_table("Gene sets higher in relapsers: dISF lesion site at V1 (GSEA)",
  g1 |> head(25) |> transmute(collection, gene_set = str_trunc(gene_set, 60), proteins = measured_proteins_in_set, NES, FDR = fmt_p(FDR),
                              leading_edge = str_trunc(leading_edge, 45)),
  subtitle = "NES > 0: higher in relapsers; sorted by p", rows_per_page = 25))
if (!is.null(p_gsea)) section("GSEA plot", page_plot(p_gsea))
if (nrow(gs_both)) section("Gene sets both", page_table("Gene sets changed in the same direction in dISF and serum at V1",
  gs_both |> head(30) |> transmute(direction = str_remove(direction, " in relapsers"), category, collection, gene_set = str_trunc(pretty_gs(gene_set), 55),
                                   NES_dISF, FDR_dISF = fmt_p(FDR_dISF), NES_serum, FDR_serum = fmt_p(FDR_serum)),
  subtitle = "GSEA relapse vs non-relapse in each matrix; NES > 0 = higher in relapsers", rows_per_page = 30))
section("Gene sets V2", page_table("Gene sets higher in relapsers after treatment: cleared lesion site at V2 (GSEA)",
  g2 |> head(20) |> transmute(collection, gene_set = str_trunc(gene_set, 60), proteins = measured_proteins_in_set, NES, FDR = fmt_p(FDR),
                              leading_edge = str_trunc(leading_edge, 45)),
  subtitle = "NES > 0: higher in relapsers; sorted by p", rows_per_page = 20))
if (nrow(cand_down)) section("Lower", page_table("Markers lower in relapsers at V1",
  cand_down |> head(20) |> transmute(rank, marker = Assay, correct = sprintf("%d/%d", loo_correct, loo_n), AUC, log2FC = round(logFC, 2), p = fmt_p(p),
                                     exact_p = fmt_p(p_exact_down), V2 = yn(persists_V2_down), serum = yn(also_serum_down)), rows_per_page = 20))
section("Methods", page_text("Methods", c(
  "Data: Olink Explore HT, PCNormalizedNPX, dISF of the MicroAD AD patients with a relapse / non-relapse label. Proteins above LOD in >= 50% of the samples of at least one group were tested.",
  "Differential test: limma, relapse vs non-relapse, at the V1 lesion site (primary), V1 non-lesional skin, V1 lesion minus non-lesional of the same patient, the V2 cleared lesion site (before any relapse) and V1 serum.",
  sprintf("Exact relabelling: the V1 lesion analysis was repeated for all %d ways to choose the relapse group; exact p per protein and the number of hits under random labels.", rc$n_labellings),
  "Marker criteria: higher in relapsers at the V1 lesion site with p < 0.05 and exact p < 0.05, above LOD in >= 50% of the relapsers. Ranked by patients predicted correctly (leave one out), AUC, persistence at V2, leave-one-out stability, fold change.",
  "Prediction: per marker, each patient classified with the cut-off (midpoint of the group means) set on the other patients. Whole-search test: the top 1 / top 5 proteins picked on the other patients, the left-out patient classified, repeated for every random labelling.",
  "Serum: the same screen on the V1 serum of the same patients. dISF and serum: markers in both matrices; the chance of a protein being higher with p < 0.05 in both is checked with the same random labels applied to both matrices.",
  "Gene sets: GSEA (fgsea) on the limma t statistics, MSigDB Hallmark, Reactome, GO biological processes and KEGG.",
  "Annotation: Human Protein Atlas protein classes (approved / potential drug targets, predicted secreted / membrane proteins).")))
invisible(grDevices::dev.off())
msg("V1 relapse screen: %s", file.path(cfg$paths$output, out))
