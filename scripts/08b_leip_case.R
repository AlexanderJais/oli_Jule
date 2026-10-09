# 08b - One Leipzig sample against the other Leipzig samples: which proteins differ, and what the proteome
#       says about the person (sex, age, BMI, body fat, lipids, CRP, kidney function, ...).
# Only Leipzig (LEIP) serum samples are used. The case is the Leipzig sample without clinical data (leip_case$case: auto)
# or the SampleID / SubjectID given in config.yml.
# In:  output of steps 01-02; reference/ (published protein models and marker lists, see reference/README.md)
# Out: output/leip_case/<case>_profile.pdf, <case>_profile.png, <case>_results.xlsx

source("R/utils.R")
source("R/report.R")
source("R/qc_overview.R")
source("R/leip_case.R")
suppressPackageStartupMessages(library(patchwork))
cfg <- load_config()
clear_outputs(cfg, "leip_case")
lc <- cfg$leip_case %||% list()
fdr <- lc$fdr %||% 0.05; min_diff <- lc$min_diff %||% 0.5; n_perm <- lc$n_perm %||% 999
min_det <- cfg$qc$min_detect_frac %||% 0.5
r2_estimate <- lc$r2_estimate %||% 0.3; r2_tertile <- lc$r2_tertile %||% 0.1
colours <- cfg$qc_overview$colours %||% list(above = "#69005F", below = "#FF506E")
set.seed(lc$seed %||% 1)

clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
leip <- clean |> filter(matrix == "Serum", cohort == "LEIP")
if (n_distinct(leip$SampleID) < 11) { msg("Fewer than 11 Leipzig serum samples - step 08b skipped."); quit(save = "no") }
samples <- leip |> distinct(SampleID, .keep_all = TRUE) |> select(SampleID, SubjectID, plate, any_of(clinical_cols)) |>
  mutate(across(any_of(setdiff(clinical_cols, "sex")), parse_clinical), across(any_of("sex"), parse_sex))
case_id <- resolve_case(samples, lc$case %||% "auto")
if (is.null(case_id)) quit(save = "no")
case_name <- coalesce(samples$SubjectID[samples$SampleID == case_id], case_id)
ref_ids <- setdiff(samples$SampleID, case_id)
clin <- samples |> filter(SampleID %in% ref_ids)
sex <- setNames(samples$sex %||% rep(NA_character_, nrow(samples)), samples$SampleID)
msg("Case %s (%s) against %d other Leipzig samples", case_name, case_id, length(ref_ids))

X <- case_matrix(leip)
B <- case_matrix(leip |> mutate(b = as.numeric(coalesce(below_lod, FALSE))), "b") >= 0.5
B[is.na(B)] <- TRUE
assay_info <- leip |> distinct(Assay, .keep_all = TRUE) |> select(Assay, OlinkID, UniProt)

# ---- 1. which proteins differ --------------------------------------------------------------------------------------
res <- case_test(X[case_id, ], X[ref_ids, ], B[ref_ids, ], min_det, fdr, min_diff) |>
  mutate(case_below_LOD = B[case_id, Assay])
msg("%d proteins tested (detected in >= %.0f%% of the reference); %d hits (FDR < %.2f, |difference| >= %.1f NPX, outside the reference range)",
    nrow(res), 100 * min_det, sum(res$hit), fdr, min_diff)
det_ref <- colMeans(!B[ref_ids, , drop = FALSE])
only_case <- tibble(Assay = colnames(X), case_value = X[case_id, ], case_detected = !B[case_id, ], ref_detected_pct = round(100 * det_ref)) |>
  filter((case_detected & ref_detected_pct < 20) | (!case_detected & ref_detected_pct >= 80)) |>
  mutate(pattern = if_else(case_detected, "detected only in the case", "not detected only in the case")) |> arrange(pattern, Assay)

# calibration: every reference sample tested against the others in the same way
calib <- map(ref_ids, \(j) {
  r <- case_test(X[j, ], X[setdiff(ref_ids, j), ], B[setdiff(ref_ids, j), ], min_det, fdr, min_diff)
  list(summary = tibble(SampleID = j, n_hits = sum(r$hit), max_abs_z = max(abs(r$z), na.rm = TRUE)), hits = r$Assay[r$hit])
})
calib_tbl <- bind_rows(map(calib, "summary")) |>
  bind_rows(tibble(SampleID = case_id, n_hits = sum(res$hit), max_abs_z = max(abs(res$z), na.rm = TRUE))) |>
  mutate(is_case = SampleID == case_id)
case_rank_p <- mean(calib_tbl$n_hits >= sum(res$hit))           # share of all samples with at least as many hits
hit_freq <- c(table(unlist(map(calib, "hits"))))
res <- res |> mutate(hit_in_other_samples = as.integer(coalesce(unname(hit_freq[Assay]), 0L)),
                     often_extreme = hit_in_other_samples >= 2)

# clinical context: correlation of each protein with the clinical values in the reference
assoc <- clinical_assoc(X[clin$SampleID, res$Assay, drop = FALSE], clin)
ctx <- assoc |> filter(fdr < 0.05 | (parameter != "sex (M - F)" & abs(rho) >= 0.5)) |>
  arrange(Assay, p) |> group_by(Assay) |>
  summarise(clinical_link = paste(sprintf("%s %s%.2f", parameter, if_else(parameter == "sex (M - F)", "diff ", "rho "), rho), collapse = "; "), .groups = "drop")
res <- res |> left_join(ctx, by = "Assay") |> left_join(assay_info, by = "Assay") |>
  relocate(OlinkID, UniProt, .after = Assay) |> arrange(p)

# Hallmark pathways on the case's t statistics
pw <- tryCatch({
  g <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  fgsea::fgsea(split(g$gene_symbol, g$gs_name), setNames(res$t, res$Assay), minSize = 10, maxSize = 500) |>
    as_tibble() |> transmute(pathway, size, NES = round(NES, 2), p = pval, fdr = padj,
                             leading_edge = map_chr(leadingEdge, \(x) paste(head(x, 10), collapse = ", "))) |> arrange(p)
}, error = \(e) { msg("Pathway analysis skipped: %s", conditionMessage(e)); tibble() })

# ---- 2. checks: sample quality, identity, handling -----------------------------------------------------------------
tested <- res$Assay
Xc <- X[, tested]; Xc <- apply(Xc, 2, \(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v })
Xc <- sweep(Xc, 2, apply(Xc[ref_ids, ], 2, median))
Xc <- Xc[, apply(Xc[ref_ids, ], 2, sd) > 0, drop = FALSE]
r_all <- cor(t(Xc))
r_case <- sort(r_all[case_id, ref_ids], decreasing = TRUE)
r_pairs <- r_all[ref_ids, ref_ids][upper.tri(r_all[ref_ids, ref_ids])]
# position in protein space: PCA of the reference (without each sample in turn), distance on the first 5 components
pc_dist <- function(train, x) {
  v <- apply(Xc[train, ], 2, sd) > 0
  pc <- prcomp(Xc[train, v], scale. = TRUE); k <- min(5, ncol(pc$x))
  s <- predict(pc, Xc[x, v, drop = FALSE])[, 1:k]; sqrt(sum(s^2 / pc$sdev[1:k]^2))
}
d_case <- pc_dist(ref_ids, case_id)
d_ref <- map_dbl(ref_ids, \(j) pc_dist(setdiff(ref_ids, j), j))
pca <- prcomp(Xc[ref_ids, ], scale. = TRUE)
pca_df <- as_tibble(rbind(pca$x[, 1:2], predict(pca, Xc[case_id, , drop = FALSE])[, 1:2])) |>
  mutate(SampleID = c(ref_ids, case_id), is_case = SampleID == case_id)

mk <- read_csv("reference/leip_case_markers.csv", show_col_types = FALSE)
wts <- function(t) { m <- mk |> filter(target == t); setNames(m$weight, m$assay) }
Z <- ref_z(X, ref_ids)
pct <- function(v, id = case_id) 100 * mean(v[ref_ids] < v[id], na.rm = TRUE)
score_text <- function(v) if (is.na(v[case_id])) c("markers not measured", "") else c(sprintf("%.2f", v[case_id]), sprintf("percentile %.0f", pct(v)))
check_score <- function(t) { s <- marker_score(Z, wts(t)); setNames(s, rownames(Z)) }
n_placenta <- { a <- intersect(names(wts("check_pregnancy")), colnames(X)); sum(X[case_id, a] > apply(X[ref_ids, a, drop = FALSE], 2, max, na.rm = TRUE) + 1, na.rm = TRUE) }
case_qc <- leip |> filter(SampleID == case_id) |> pull(SampleQC) |> unique()
checks <- tribble(
  ~check, ~case, ~reference, ~reading,
  "Olink sample QC", paste(case_qc, collapse = ", "), "", if (all(case_qc == "PASS")) "ok" else "warning in Olink QC",
  "plate", samples$plate[samples$SampleID == case_id], sprintf("%d of %d reference samples on the same plate", sum(samples$plate[samples$SampleID %in% ref_ids] == samples$plate[samples$SampleID == case_id]), length(ref_ids)), "",
  "proteins below LOD (%)", sprintf("%.1f", 100 * mean(B[case_id, ])), sprintf("reference %.1f-%.1f", 100 * min(rowMeans(B[ref_ids, ])), 100 * max(rowMeans(B[ref_ids, ]))),
    sprintf("percentile %.0f", pct(rowMeans(B))),
  "highest correlation with another sample", sprintf("%.2f (%s)", r_case[1], samples$SubjectID[match(names(r_case)[1], samples$SampleID)]),
    sprintf("highest between two reference samples %.2f", max(r_pairs)),
    if (r_case[1] > max(r_pairs)) "more similar to one sample than any two reference samples are: possible duplicate" else "no sign of a duplicate",
  "distance from the reference (PCA)", sprintf("%.1f", d_case), sprintf("reference %.1f-%.1f (each left out in turn)", min(d_ref), max(d_ref)),
    if (d_case > max(d_ref)) "further out than every reference sample" else sprintf("percentile %.0f", 100 * mean(d_ref < d_case)),
  "proteins differing (hits)", as.character(sum(res$hit)), sprintf("reference %d-%d (median %.0f), each tested against the others", min(calib_tbl$n_hits[!calib_tbl$is_case]), max(calib_tbl$n_hits[!calib_tbl$is_case]), median(calib_tbl$n_hits[!calib_tbl$is_case])),
    sprintf("%.0f%% of the samples have as many hits or more", 100 * case_rank_p),
  "platelet release score", score_text(check_score("check_platelet"))[1], "", score_text(check_score("check_platelet"))[2],
  "haemolysis score", score_text(check_score("check_haemolysis"))[1], "", score_text(check_score("check_haemolysis"))[2],
  "neutrophil release score", score_text(check_score("check_leukocyte"))[1], "", score_text(check_score("check_leukocyte"))[2],
  "placenta-only proteins > reference max + 1 NPX", as.character(n_placenta), "", if (n_placenta >= 3) "placental protein pattern" else "no placental pattern")

# ---- 3. what the proteome says about the person -----------------------------------------------------------------------
sx <- predict_sex(Z, B, ref_ids, case_id, sex, wts("sex"), min_det)
msg("Sex: %s (P(male) = %.2f; %d of %d reference samples correct leave-one-out)", sx$call, sx$p_male, sx$n_ref - sx$loo_errors, sx$n_ref)

usable <- function(w) w[names(w) %in% colnames(X) & names(w) %in% names(det_ref)[det_ref >= min_det]]
clock <- if (file.exists("reference/age_clock_goeminne2025.csv")) read_csv("reference/age_clock_goeminne2025.csv", show_col_types = FALSE) else NULL
bmi_w <- read_csv("reference/bmi_score_watanabe2023.csv", show_col_types = FALSE)
analyte <- apo_analyte(clin, lc$apo_analyte %||% "auto")
spec <- list(
  list(t = "age", label = "Age (years)", w = if (!is.null(clock)) c(tapply(clock$beta * clock$ukb_sd, clock$assay, sum)) else usable(wts("age_markers")),
       clock = !is.null(clock), model = if (!is.null(clock)) "published UK Biobank age clock (Goeminne 2025)" else "age marker score (clock file missing)",
       bands = list(c(40, 60), c("<40", "40-59", ">=60"))),
  list(t = "BMI", label = "BMI (kg/m2)", w = setNames(bmi_w$beta, bmi_w$assay), clock = TRUE, log = TRUE, sex = TRUE,
       model = "published protein BMI score (Watanabe 2023, 67-protein version of Wang 2024) + sex", bands = list(c(18.5, 25, 30), c("<18.5", "18.5-25", "25-30", ">=30"))),
  list(t = "c_fett", label = "Body fat (c_fett)", w = usable(wts("c_fett")), sex = TRUE, model = "leptin + sex"),
  list(t = "WHR", label = "Waist-hip ratio", type = "sex_only"),
  list(t = "HOMA_IR", label = "HOMA-IR", w = usable(wts("HOMA_IR")), log = TRUE, model = "insulin resistance markers"),
  list(t = "Gluc0_mg_dl", label = "Fasting glucose (mg/dl)", type = "none"),
  list(t = "c_CRP", label = "CRP (c_CRP)", w = usable(wts("c_CRP")), log = TRUE, model = "acute-phase proteins (CRP itself is not on the panel)",
       bands = list(c(1, 3), c("<1", "1-3", ">3"))),
  list(t = "C_CHOL", label = "Total cholesterol", w = usable(wts("C_CHOL")), model = "lipoprotein-carried proteins (UK Biobank weights)"),
  list(t = "C_HDL", label = "HDL cholesterol", w = usable(wts("C_HDL")), sex = TRUE, model = "HDL proteins (UK Biobank weights) + sex"),
  list(t = "C_LDL", label = "LDL cholesterol", w = usable(wts("C_LDL")), model = "LDL-carried proteins (UK Biobank weights)"),
  list(t = "C_TRIGLY", label = "Triglycerides", w = usable(wts("C_TRIGLY")), log = TRUE, model = "triglyceride-linked proteins (UK Biobank weights)"),
  list(t = "c_apo", label = sprintf("Apolipoprotein (c_apo%s)", if (is.na(analyte)) "" else paste0(", taken as ", analyte)),
       w = if (is.na(analyte)) NULL else usable(wts(paste0("apo_", analyte))), log = identical(analyte, "Lpa"),
       model = if (is.na(analyte)) "analyte unclear (set leip_case$apo_analyte)" else paste(analyte, "proteins"), type = if (is.na(analyte)) "none" else "model"),
  list(t = "MDRD_kurz", label = "eGFR (MDRD)", w = usable(wts("MDRD_kurz")), log = TRUE, sex = TRUE, model = "kidney filtration proteins + sex"))
spec <- keep(spec, \(s) s$t %in% names(clin) && sum(!is.na(clin[[s$t]])) >= 10)

fits <- list(); rows <- list(); loo_tbl <- list()
for (s in spec) {
  y <- setNames(clin[[s$t]], clin$SampleID)
  base <- tibble(parameter = s$t, label = s$label, n_ref = sum(!is.na(y)), ref_median = median(y, na.rm = TRUE),
                 ref_min = min(y, na.rm = TRUE), ref_max = max(y, na.rm = TRUE))
  type <- s$type %||% "model"
  if (type == "model" && length(s$w) == 0) type <- "none"
  if (type == "none") { rows[[s$t]] <- base |> mutate(verdict = "not predictable", how = s$model %||% "no protein model"); next }
  if (type == "sex_only") {
    by_sex <- clin |> filter(!is.na(.data[[s$t]]), !is.na(sex)) |> group_by(sex) |>
      summarise(m = median(.data[[s$t]]), lo = quantile(.data[[s$t]], 0.1), hi = quantile(.data[[s$t]], 0.9), .groups = "drop")
    pick <- if (sx$call == "male") "M" else if (sx$call == "female") "F" else c("M", "F")
    b <- by_sex |> filter(sex %in% pick)
    rows[[s$t]] <- base |> mutate(verdict = "sex-specific average only", how = "median and 10-90% range of the reference samples of the estimated sex",
                                  estimate = if (nrow(b) == 1) b$m else NA_real_, lo80 = min(b$lo), hi80 = max(b$hi))
    next
  }
  sc <- marker_score(Z, s$w, clock = isTRUE(s$clock))
  names(sc) <- rownames(Z)
  coverage <- if (isTRUE(s$clock)) sum(abs(s$w[names(s$w) %in% colnames(X)])) / sum(abs(s$w)) else 1
  fp <- fit_profile(y[clin$SampleID], sc[clin$SampleID], if (isTRUE(s$sex)) sex[clin$SampleID] else NULL,
                    sc[case_id], log = isTRUE(s$log), n_perm = n_perm)
  if (is.null(fp)) { rows[[s$t]] <- base |> mutate(verdict = "not predictable", how = paste(s$model, "(model could not be fitted)")); next }
  pr <- case_prediction(fp, sx$call, sx$p_male)
  fits[[s$t]] <- list(fp = fp, pr = pr, s = s)
  loo_tbl[[s$t]] <- fp$loo |> mutate(parameter = s$t)
  terc <- quantile(y, c(1, 2) / 3, na.rm = TRUE)
  rows[[s$t]] <- base |> mutate(
    how = s$model, n_markers = sum(names(s$w) %in% colnames(X)), marker_weight_measured_pct = round(100 * coverage),
    estimate = pr$est, lo80 = pr$lo, hi80 = pr$hi,
    third_of_reference = c("low", "middle", "high")[findInterval(pr$est, terc) + 1],
    bands = if (!is.null(s$bands)) band_probs(pr$dist, pr$w, s$bands[[1]], s$bands[[2]]) else NA_character_,
    r2_cv = fp$r2cv, mae = fp$mae, mae_guessing_mean = fp$mae_baseline, p_perm = fp$p_perm,
    case_score_percentile = round(fp$score_pct), case_score_outside_reference = fp$score_outside)
}
profile <- bind_rows(rows)
if (!"verdict" %in% names(profile)) profile$verdict <- NA_character_
if ("p_perm" %in% names(profile)) profile <- profile |>
  mutate(q = p.adjust(p_perm, "BH"),
         verdict = case_when(!is.na(verdict) ~ verdict,
                             q < 0.05 & r2_cv >= r2_estimate ~ "estimate",
                             q < 0.05 & r2_cv >= r2_tertile ~ "low / middle / high only",
                             TRUE ~ "not predictable"))
sex_row <- tibble(parameter = "sex", label = "Sex", verdict = sx$call,
                  how = sprintf("%d sex-specific proteins; P(male) = %.2f%s", nrow(sx$markers), sx$p_male, if (nzchar(sx$reason)) paste0("; ", sx$reason) else ""),
                  n_ref = sx$n_ref, r2_cv = NA_real_,
                  loo_correct = sprintf("%d of %d", sx$n_ref - coalesce(sx$loo_errors, sx$n_ref), sx$n_ref))
profile <- bind_rows(sex_row, profile)
# the shared table never shows a number for a value the proteins cannot predict
profile_out <- profile |>
  mutate(across(any_of(c("estimate", "lo80", "hi80", "third_of_reference", "bands")),
                \(v) if_else(verdict %in% c("estimate", "low / middle / high only", "sex-specific average only"), v, NA)),
         across(any_of(c("estimate", "lo80", "hi80")), \(v) if_else(verdict == "low / middle / high only", NA_real_, v)),
         across(where(is.double), \(v) signif(v, 3)))
for (i in seq_len(nrow(profile_out))) msg("  %-28s %s", profile_out$label[i], profile_out$verdict[i])

# single markers worth a look, as percentile among the reference
panel <- c("LEP", "FABP4", "IGFBP1", "SHBG", "FSHB", "LHB", "GDF15", "NPPB", "CST3", "LPA", "PCSK9", "IL6", "CGB3_CGB5_CGB8")
markers_tbl <- tibble(Assay = intersect(panel, colnames(X))) |>
  mutate(case_value = X[case_id, Assay], case_below_LOD = B[case_id, Assay],
         percentile_in_reference = map_dbl(Assay, \(a) round(pct(X[, a]))),
         ref_median = map_dbl(Assay, \(a) median(X[ref_ids, a], na.rm = TRUE)),
         ref_max = map_dbl(Assay, \(a) max(X[ref_ids, a], na.rm = TRUE)))

# ---- Excel ------------------------------------------------------------------------------------------------------------
readme <- tibble(sheet = c("profile", "sex_markers", "proteins", "hits", "only_in_case", "calibration", "checks", "pathways", "markers", "model_check", "how"),
  content = c(sprintf("what the proteome says about %s: estimate with 80%% range, or low/middle/high third of the reference, or 'not predictable'", case_name),
              "the sex-specific proteins: value of the case (z vs the reference) and which sex each one points to",
              sprintf("all %d tested proteins: Crawford-Howell single-case test against the %d other Leipzig samples", nrow(res), length(ref_ids)),
              sprintf("proteins that differ: FDR < %.2f, |difference| >= %.1f NPX (log2) and outside the range of the reference", fdr, min_diff),
              "proteins detected only in the case, or not detected only in the case (descriptive)",
              "number of hits when every reference sample is tested against the others in the same way",
              "sample quality, duplicate check, handling scores",
              "Hallmark pathways ranked by the case's t statistics (GSEA)",
              "selected single proteins as percentile among the reference",
              "leave-one-out predictions of each reference sample: how well each model works",
              paste("Estimates are calibrated on the Leipzig samples with clinical data, each model fixed in advance (reference/README.md).",
                    sprintf("A value is given as a number with an 80%% range only if the model explains >= %.0f%% of the variation in leave-one-out tests (r2_cv) and beats chance (permutation, BH q < 0.05);", 100 * r2_estimate),
                    sprintf("low/middle/high third if it explains %.0f-%.0f%%; otherwise 'not predictable'. Research use only.", 100 * r2_tertile, 100 * r2_estimate))))
writexl::write_xlsx(list(README = readme, profile = profile_out, sex_markers = sx$markers |> mutate(across(where(is.double), \(v) round(v, 2))),
                         proteins = res, hits = res |> filter(hit), only_in_case = only_case, calibration = calib_tbl, checks = checks,
                         pathways = pw, markers = markers_tbl, model_check = bind_rows(loo_tbl) |> relocate(parameter, SampleID)),
                    out_path(cfg, "leip_case", paste0(case_name, "_results.xlsx")))

# ---- figures ------------------------------------------------------------------------------------------------------------
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
violet <- colours$above; pink <- colours$below
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA))
fmt <- function(v) formatC(v, digits = 3, format = "fg", flag = "#") |> str_remove("\\.$") |> trimws()

prof_rows <- profile_out |> filter(parameter != "sex")
ref_long <- clin |> select(SampleID, any_of(prof_rows$parameter)) |> pivot_longer(-SampleID, names_to = "parameter") |>
  filter(!is.na(value)) |> left_join(prof_rows |> select(parameter, label), by = "parameter")
labs_df <- prof_rows |> mutate(
  text = case_when(verdict == "estimate" ~ sprintf("%s (80%%: %s-%s)", fmt(estimate), fmt(lo80), fmt(hi80)),
                   verdict == "low / middle / high only" ~ paste(third_of_reference, "third"),
                   verdict == "sex-specific average only" ~ if_else(is.na(estimate), "sex unclear", sprintf("about %s (sex average)", fmt(estimate))),
                   TRUE ~ "not predictable"),
  text = paste0(text, if_else(!is.na(r2_cv), sprintf("   r2_cv %.2f", r2_cv), "")),
  strip = paste0(label, ": ", text))
ref_long <- ref_long |> left_join(labs_df |> select(parameter, strip), by = "parameter")
case_pts <- labs_df |> filter(verdict %in% c("estimate", "sex-specific average only"), !is.na(estimate))
case_rng <- labs_df |> filter(verdict %in% c("estimate", "sex-specific average only", "low / middle / high only")) |>
  mutate(lo = if_else(verdict == "low / middle / high only", NA_real_, lo80), hi = if_else(verdict == "low / middle / high only", NA_real_, hi80))
terc_band <- ref_long |> group_by(parameter, strip) |> summarise(t1 = quantile(value, 1 / 3), t2 = quantile(value, 2 / 3), lo = min(value), hi = max(value), .groups = "drop") |>
  inner_join(labs_df |> filter(verdict == "low / middle / high only") |> select(parameter, third_of_reference), by = "parameter") |>
  mutate(xmin = case_when(third_of_reference == "low" ~ lo, third_of_reference == "middle" ~ t1, TRUE ~ t2),
         xmax = case_when(third_of_reference == "low" ~ t1, third_of_reference == "middle" ~ t2, TRUE ~ hi))
p_prof <- ggplot(ref_long, aes(value, 0)) +
  geom_rect(data = terc_band, aes(xmin = xmin, xmax = xmax, ymin = -0.45, ymax = 0.45), inherit.aes = FALSE, fill = violet, alpha = 0.18) +
  geom_jitter(height = 0.25, width = 0, colour = "grey55", size = 1.4) +
  geom_errorbar(data = case_rng |> filter(!is.na(lo)), aes(xmin = lo, xmax = hi, y = 0), inherit.aes = FALSE, orientation = "y", width = 0.35, colour = violet, linewidth = 0.8) +
  geom_point(data = case_pts |> mutate(value = estimate), colour = violet, size = 3.2, shape = 18) +
  facet_wrap(~strip, scales = "free_x", ncol = 2) + scale_y_continuous(NULL, breaks = NULL, limits = c(-0.5, 0.5)) +
  labs(x = NULL) + th + theme(strip.text = element_text(hjust = 0, size = 9))
sexd <- sx$ref_scores |> mutate(grp = recode(sex, M = "male", F = "female", case = case_name))
p_sex <- if (nrow(sexd)) ggplot(sexd, aes(score, grp, colour = grp == case_name)) + geom_jitter(height = 0.15, width = 0, size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = violet), guide = "none") +
  labs(x = "sex protein score (positive = male pattern)", y = NULL,
       title = sprintf("Sex: %s (P(male) = %.2f; %d of %d reference samples classified correctly)", sx$call, sx$p_male, sx$n_ref - sx$loo_errors, sx$n_ref)) + th else NULL
page1_note <- patchwork::plot_annotation(title = sprintf("%s: what the proteome says (calibrated on %d Leipzig samples)", case_name, nrow(clin)),
                             subtitle = "grey: reference samples; violet: estimate for the case with 80% range; shaded: estimated third of the reference")

top <- res |> filter(hit) |> slice_head(n = 12)
if (!nrow(top)) top <- res |> slice_head(n = 12)
lab <- res |> filter(hit) |> slice_head(n = 25)
p_volc <- ggplot(res, aes(diff, -log10(p))) + geom_point(aes(colour = hit), size = 1) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = violet), guide = "none") +
  geom_text(data = lab, aes(label = Assay), size = 2.6, vjust = -0.6, check_overlap = TRUE) +
  labs(x = "difference to the reference median (NPX, log2)", y = "-log10 p",
       title = sprintf("%d of %d proteins differ (FDR < %.2f, >= %.1f NPX, outside the reference range)", sum(res$hit), nrow(res), fdr, min_diff)) + th
strip_df <- map(top$Assay, \(a) tibble(Assay = a, value = X[c(ref_ids, case_id), a], is_case = c(ref_ids, case_id) == case_id)) |> bind_rows() |>
  mutate(Assay = factor(Assay, top$Assay))
p_strip <- ggplot(strip_df, aes(Assay, value)) + geom_jitter(data = \(d) filter(d, !is_case), width = 0.15, height = 0, colour = "grey55", size = 1.2) +
  geom_point(data = \(d) filter(d, is_case), colour = violet, size = 3, shape = 18) +
  labs(x = NULL, y = "NPX", title = if (any(res$hit)) "Strongest differences" else "No protein differs - lowest p values shown") + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cor <- ggplot(tibble(r = r_pairs), aes(r)) + geom_histogram(bins = 40, fill = "grey75") +
  geom_vline(xintercept = r_case[1:3], colour = violet) +
  labs(x = "correlation between two samples (centred protein profiles)", y = "reference pairs",
       title = sprintf("Duplicate check: closest samples to %s (violet)", case_name)) + th
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = is_case, size = is_case)) + geom_point() +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = violet), guide = "none") + scale_size_manual(values = c(1.5, 3.5), guide = "none") +
  labs(title = "Position among the Leipzig samples (PCA)") + th
p_cal <- ggplot(calib_tbl, aes(n_hits, fill = is_case)) + geom_histogram(binwidth = 1) +
  scale_fill_manual(values = c(`FALSE` = "grey75", `TRUE` = violet), guide = "none") +
  labs(x = "number of differing proteins when tested against the others", y = "samples", title = "Is the case more unusual than the others?") + th

fam <- fi$family
sf <- function(p) set_family(p, fam)
page1 <- patchwork::wrap_plots(compact(list(sf(p_sex), sf(p_prof))), ncol = 1, heights = if (is.null(p_sex)) 1 else c(1, 5)) + page1_note
page2 <- sf(p_volc) / sf(p_strip)
page3 <- (sf(p_cor) | sf(p_pca)) / sf(p_cal)
if (nzchar(fam)) { page1 <- page1 & theme(text = element_text(family = fam)); page2 <- page2 & theme(text = element_text(family = fam)); page3 <- page3 & theme(text = element_text(family = fam)) }
pdf_file <- out_path(cfg, "leip_case", paste0(case_name, "_profile.pdf"))
open_device(pdf_file, "pdf", fi)
for (pg in list(page1, page2, page3)) print(pg)
close_device(fi)
open_device(out_path(cfg, "leip_case", paste0(case_name, "_profile.png")), "png", fi)
print(page1)
close_device(fi)
msg("Leipzig case analysis: %s", file.path(cfg$paths$output, "leip_case"))
