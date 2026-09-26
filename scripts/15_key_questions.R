# 15 - Key questions of the project, answered with pre-specified tests
#   Q1  Are mast cell markers elevated in AD, or only in relapse vs non-relapse?
#   Q2  Is TNFRSF9 (CD137) or its ligand TNFSF9 a marker for mast cells in AD?
#   Q3  Do TNFRSF9 / TNFSF9 correlate with AD relapse?
#   Q4  Are TNFRSF9 / TNFSF9 a marker (change with lesion activity) or a predictor (values
#       before the relapse separate relapsers from non-relapsers) of AD relapse?
#   Q5  Is dISF superior to serum?
# Proteins and markers are set in config.yml (key_questions). All tests are single, pre-specified
# tests (unadjusted p-values); relapse analyses are exploratory (4 relapsers in MicroAD).
# Out: output/key_questions/  answers.csv, key_questions.xlsx, figures

source("R/utils.R")
source("R/models.R")
source("R/design.R")
set.seed(1)
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
clear_outputs(cfg, "key_questions")
kq    <- cfg$key_questions %||% list()
mast  <- unlist(kq$mast_cell_markers %||% c("KITLG", "CPA4", "FCER1A", "TPSAB1", "MS4A2", "TPSD1"))
kprot <- unlist(kq$proteins %||% c("TNFRSF9", "TNFSF9"))
B     <- kq$bootstrap %||% 2000
mgn   <- cfg$stats$min_group_n
fmt_p <- \(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
fmt_e <- \(e, p) ifelse(is.na(e), "n/a", sprintf("%+.2f (p = %s)", e, fmt_p(p)))

# ---- values: proteins and mast cell score ------------------------------------------------------------------
find_oid <- function(name) {
  N <- toupper(name)
  hit <- clean |> distinct(OlinkID, Assay) |>
    filter(toupper(Assay) == N | str_detect(toupper(Assay), paste0("(^|_)", N, "($|_)")) | OlinkID == name)
  if (nrow(hit)) hit$OlinkID[1] else NA_character_
}
vals_of <- \(oid) clean |> filter(OlinkID == oid) |> select(SampleID, matrix, value)
mast_oid <- setNames(vapply(mast, find_oid, ""), mast)
kp_oid   <- setNames(vapply(kprot, find_oid, ""), kprot)
msg("Mast cell markers found: %s; not in data: %s", paste(names(mast_oid)[!is.na(mast_oid)], collapse = ", "),
    if (any(is.na(mast_oid))) paste(names(mast_oid)[is.na(mast_oid)], collapse = ", ") else "-")
mast_oid <- mast_oid[!is.na(mast_oid)]; kp_oid <- kp_oid[!is.na(kp_oid)]
if (!length(kp_oid)) stop("None of the key proteins (", paste(kprot, collapse = ", "), ") is in the data.")

# mast cell score = mean of the z-scored markers (per matrix), at least 2 markers per sample
score <- if (length(mast_oid) >= 2) clean |> filter(OlinkID %in% mast_oid) |>
  group_by(matrix, OlinkID) |> mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)) |>
  group_by(SampleID, matrix) |>
  summarise(value = if (sum(!is.na(z)) >= 2) mean(z, na.rm = TRUE) else NA_real_, .groups = "drop") else NULL
entities <- c(map(mast_oid, vals_of), if (!is.null(score)) list(`Mast cell score` = score), map(kp_oid, vals_of))

isf_base <- isf_design(meta); ser_base <- serum_design(meta)
tests <- imap(entities, \(v, nm) {
  vv <- v |> select(SampleID, value)
  prespecified_tests(isf_base |> inner_join(vv, by = "SampleID"), ser_base |> inner_join(vv, by = "SampleID"), mgn) |>
    mutate(entity = nm, .before = 1)
}) |> bind_rows()
tt <- \(ent, mdl, ct) { r <- tests |> filter(entity == ent, model == mdl, contrast == ct); if (nrow(r)) r[1, ] else NULL }
est_p <- \(ent, mdl, ct) { r <- tt(ent, mdl, ct); c(e = r$estimate %||% NA_real_, p = r$p %||% NA_real_) }

# ---- Q1: mast cell markers - AD and/or relapse ----------------------------------------------------------------
q1_cmp <- tribble(
  ~label,                               ~model,                      ~contrast,        ~type,
  "dISF lesional vs healthy",           "states_all_visits",         "AD_L_vs_HC",     "AD",
  "dISF non-lesional vs healthy",       "states_all_visits",         "AD_NL_vs_HC",    "AD",
  "serum AD vs healthy",                "AD_vs_HC_in_study",         "AD_vs_HC",       "AD",
  "dISF relapse (ex-lesional)",         "relapse_ex_lesional",       "relapse_vs_non", "relapse",
  "dISF relapse (xL - NL)",             "relapse_delta_xL_minus_NL", "relapse_vs_non", "relapse",
  "serum relapse (MicroAD)",            "MicroAD_relapse",           "relapse_vs_non", "relapse",
  "serum relapse (RELAD/RELAD2)",       "RELAD_relapse",             "relapse_vs_non", "relapse")
q1_ent <- c(names(mast_oid), if (!is.null(score)) "Mast cell score")
q1 <- tidyr::expand_grid(entity = q1_ent, q1_cmp) |>
  mutate(ep = pmap(list(entity, model, contrast), est_p), estimate = map_dbl(ep, "e"), p = map_dbl(ep, "p")) |>
  select(-ep)
q1_verdict <- q1 |> group_by(entity) |>
  summarise(in_AD = any(type == "AD" & p < 0.05 & estimate > 0, na.rm = TRUE),
            where_AD = paste(label[type == "AD" & p < 0.05 & estimate > 0 & !is.na(p)], collapse = "; "),
            in_relapse = any(type == "relapse" & p < 0.05, na.rm = TRUE),
            where_relapse = paste(label[type == "relapse" & p < 0.05 & !is.na(p)], collapse = "; "), .groups = "drop") |>
  mutate(answer = case_when(in_AD & in_relapse ~ "elevated in AD AND associated with relapse",
                            in_AD ~ "elevated in AD, not associated with relapse",
                            in_relapse ~ "associated with relapse only, not elevated in AD",
                            TRUE ~ "neither elevated in AD nor associated with relapse"))
p <- ggplot(q1, aes(label, factor(entity, levels = rev(q1_ent)), fill = estimate)) + geom_tile() +
  geom_text(aes(label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")), size = 4) +
  scale_fill_gradient2(low = "steelblue", high = "firebrick") + facet_grid(~type, scales = "free_x", space = "free_x") +
  labs(title = "Q1  Mast cell markers: elevated in AD, or associated with relapse?",
       subtitle = "effect (log2 / score units); * p < 0.05, ** < 0.01, *** < 0.001 (single pre-specified tests)", x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, cfg, "key_questions", "Q1_mast_cell_markers.png", width = 11, height = 2.5 + 0.45 * length(q1_ent))

# ---- Q2: CD137 / CD137L vs mast cells in AD dISF -------------------------------------------------------------------
ad_isf <- meta |> filter(matrix == "ISF", group == "AD") |> select(SampleID, SubjectID, site, state)
corr_one <- function(x, y, subj) {
  ok <- !is.na(x) & !is.na(y)
  s <- if (sum(ok) >= 5) suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)) else NULL
  multi <- ok & subj %in% names(which(table(subj[ok]) >= 2))
  w <- if (sum(multi) >= 5 && n_distinct(subj[multi]) >= 2) tryCatch(suppressWarnings(
    rmcorr::rmcorr(participant = subj, measure1 = x, measure2 = y,
                   dataset = data.frame(subj = factor(subj[multi]), x = x[multi], y = y[multi]))), error = \(e) NULL) else NULL
  tibble(n = sum(ok), rho = unname(s$estimate %||% NA), p_rho = s$p.value %||% NA_real_,
         r_within = w$r %||% NA_real_, p_within = w$p %||% NA_real_)
}
partners <- c(if (!is.null(score)) "Mast cell score", names(mast_oid))
q2 <- tidyr::expand_grid(protein = names(kp_oid), partner = partners, subset = c("all AD dISF", "lesion site", "non-lesional")) |>
  pmap(\(protein, partner, subset) {
    a <- entities[[protein]] |> filter(matrix == "ISF") |> select(SampleID, x = value)
    b <- entities[[partner]] |> filter(matrix == "ISF") |> select(SampleID, y = value)
    d <- ad_isf |> inner_join(a, by = "SampleID") |> inner_join(b, by = "SampleID") |>
      filter(subset == "all AD dISF" | (subset == "lesion site" & site == "L") | (subset == "non-lesional" & site == "NL"))
    corr_one(d$x, d$y, d$SubjectID) |> mutate(protein = protein, partner = partner, subset = subset, .before = 1)
  }) |> bind_rows()
if (!is.null(score)) {
  sc <- tidyr::expand_grid(protein = names(kp_oid)) |> pmap(\(protein) {
    ad_isf |> inner_join(entities[[protein]] |> filter(matrix == "ISF") |> select(SampleID, x = value), by = "SampleID") |>
      inner_join(score |> filter(matrix == "ISF") |> select(SampleID, y = value), by = "SampleID") |> mutate(protein = protein)
  }) |> bind_rows() |>
    mutate(state = factor(state, levels = c("non-lesional", "ex-lesional", "lesional")))
  p <- ggplot(sc, aes(y, x, colour = state)) + geom_point(size = 1.4, alpha = 0.8) +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, colour = "grey30", linewidth = 0.5) +
    scale_colour_manual(values = c(`non-lesional` = "steelblue", `ex-lesional` = "orange3", lesional = "firebrick")) +
    facet_wrap(~protein, scales = "free_y") +
    labs(title = "Q2  CD137 / CD137L vs mast cell score in AD dISF", x = "mast cell score (mean z)", y = "NPX", colour = "skin state")
  save_plot(p, cfg, "key_questions", "Q2_cd137_vs_mast_score.png", width = 10, height = 5)
}

# ---- Q3: correlation with relapse -------------------------------------------------------------------------------------
rel_date <- meta |> filter(!is.na(relapse_visit), visit_num == relapse_visit) |> group_by(SubjectID) |>
  summarise(relapse_date = min(date, na.rm = TRUE))
q3_trend <- map(c(names(kp_oid), if (!is.null(score)) "Mast cell score"), \(ent) {
  v <- entities[[ent]] |> select(SampleID, value)
  isf <- isf_base |> inner_join(v, by = "SampleID") |> left_join(rel_date, by = "SubjectID") |>
    mutate(weeks_to_relapse = as.numeric(date - relapse_date) / 7)
  pr <- isf_delta_pairs(isf) |> left_join(isf |> select(AD_xL = SampleID, weeks_to_relapse), by = "AD_xL")
  pr$value <- isf$value[match(pr$AD_xL, isf$SampleID)] - isf$value[match(pr$AD_NL, isf$SampleID)]
  ser <- ser_base |> inner_join(v, by = "SampleID") |> left_join(rel_date, by = "SubjectID") |>
    mutate(weeks_to_relapse = as.numeric(date - relapse_date) / 7) |>
    filter(cohort == "MicroAD", lesion_state %in% "cleared", !is.na(weeks_to_relapse))
  bind_rows(
    test_single(pr |> filter(!is.na(weeks_to_relapse)), ~ weeks_to_relapse + (1 | SubjectID), c(per_week = "weeks_to_relapse"),
                "dISF xL - NL vs weeks to relapse (relapsers)", mgn),
    test_single(ser, ~ weeks_to_relapse + plate + (1 | SubjectID), c(per_week = "weeks_to_relapse"),
                "serum vs weeks to relapse (relapsers)", mgn)) |> mutate(entity = ent, .before = 1)
}) |> bind_rows()
q3 <- tests |> filter(entity %in% c(names(kp_oid), "Mast cell score"), contrast == "relapse_vs_non") |>
  select(entity, matrix, model, contrast, estimate, ci_low, ci_high, p, n_samples, n_subjects) |>
  bind_rows(q3_trend |> transmute(entity, matrix = if_else(str_detect(model, "^dISF"), "dISF", "serum"), model, contrast,
                                  estimate, ci_low, ci_high, p, n_samples, n_subjects))

# ---- Q4: marker vs predictor (AUC) -------------------------------------------------------------------------------------
auc_ci <- function(x, g) {
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- g[ok]
  a <- x[g]; b <- x[!g]
  if (length(a) < 2 || length(b) < 2) return(tibble(n_relapse = length(a), n_non = length(b), AUC = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  auc <- \(a, b) mean(outer(a, b, ">") + 0.5 * outer(a, b, "=="))
  bs <- replicate(B, auc(sample(a, replace = TRUE), sample(b, replace = TRUE)))
  tibble(n_relapse = length(a), n_non = length(b), AUC = auc(a, b), ci_low = unname(quantile(bs, 0.025)),
         ci_high = unname(quantile(bs, 0.975)), p = suppressWarnings(wilcox.test(a, b))$p.value)
}
q4_sets <- function(ent) {
  v <- entities[[ent]] |> select(SampleID, value)
  isf <- isf_base |> inner_join(v, by = "SampleID") |> filter(group == "AD", !is.na(relapse2))
  ser <- ser_base |> inner_join(v, by = "SampleID") |> filter(group == "AD", !is.na(relapse2))
  pr <- isf_delta_pairs(isf_base |> inner_join(v, by = "SampleID"))
  pr$value <- isf$value[match(pr$AD_xL, isf$SampleID)] - isf$value[match(pr$AD_NL, isf$SampleID)]
  per_patient <- \(d) d |> group_by(SubjectID, relapse2) |> summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
  sets <- list(
    `dISF lesion site at V1 (baseline)` = isf |> filter(site == "L", visit_num == 1) |> per_patient(),
    `dISF ex-lesional, before relapse (patient mean)` = isf |> filter(state == "ex-lesional", pre_relapse | is.na(relapse_visit)) |> per_patient(),
    `dISF xL - NL, before relapse (patient mean)` = pr |> filter(!is.na(value)) |> per_patient(),
    `serum MicroAD, cleared visits before relapse (patient mean)` = ser |> filter(cohort == "MicroAD", lesion_state %in% "cleared",
                                                                                 pre_relapse | is.na(relapse_visit)) |> per_patient(),
    `serum RELAD/RELAD2 (one sample per patient)` = ser |> filter(cohort %in% c("RELAD", "RELAD2")) |> per_patient())
  imap(sets, \(d, nm) auc_ci(d$value, d$relapse2 == "relapse") |> mutate(entity = ent, predictor = nm, .before = 1)) |> bind_rows()
}
q4_auc <- map(c(names(kp_oid), if (!is.null(score)) "Mast cell score"), q4_sets) |> bind_rows() |>
  mutate(matrix = if_else(str_detect(predictor, "^dISF"), "dISF", "serum"),
         direction = if_else(AUC > 0.5, "higher", "lower"),
         # with < 10 patients per group a single separation can occur by chance (4 vs 6: p = 0.01 for a
         # perfect split) and the bootstrap CI is unreliable -> only "possible", never a firm call
         verdict = case_when(is.na(AUC) ~ "not estimable",
                             p < 0.05 & pmin(n_relapse, n_non) >= 10 & (ci_low > 0.5 | ci_high < 0.5) ~
                               sprintf("predicts relapse (%s before relapse)", direction),
                             p < 0.05 ~ sprintf("possible predictor (%s before relapse; small groups, exploratory)", direction),
                             TRUE ~ "no predictive value shown"))
q4_marker <- tests |> filter(entity %in% c(names(kp_oid), "Mast cell score"),
                             paste(model, contrast) %in% c("states_all_visits AD_L_vs_xL", "states_all_visits AD_L_vs_NL",
                                                           "MicroAD_active_vs_cleared active_vs_cleared")) |>
  select(entity, matrix, model, contrast, estimate, ci_low, ci_high, p)
p <- ggplot(q4_auc |> filter(!is.na(AUC)), aes(AUC, predictor, colour = matrix)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high)) + facet_wrap(~entity) +
  scale_colour_manual(values = c(dISF = "firebrick", serum = "steelblue")) + coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Q4  Do values BEFORE the relapse separate relapsers from non-relapsers?",
       subtitle = "AUC with 95% bootstrap CI; 0.5 = no separation, >0.5 = higher in later relapsers", x = "AUC", y = NULL, colour = NULL)
save_plot(p, cfg, "key_questions", "Q4_relapse_prediction_auc.png", width = 12, height = 3 + 1.2 * ceiling(n_distinct(q4_auc$entity) / 3))

# ---- Q5: dISF vs serum ------------------------------------------------------------------------------------------------
det <- read_csv(file.path(cfg$paths$output, "qc", "assay_detection.csv"), show_col_types = FALSE)
q5_det <- det |> group_by(matrix) |> summarise(proteins_measured = n(), proteins_detected = sum(keep))
ov_path <- file.path(cfg$paths$output, "serum_vs_disf", "overlap_summary.csv")
q5_overlap <- if (file.exists(ov_path)) read_csv(ov_path, show_col_types = FALSE) |> filter(visit == "all visits") |>
  select(tier, question, isf_site, dISF_significant, serum_significant, both_same, both_opposite, dISF_only, dISF_only_not_in_serum, serum_only) else NULL
q5_effects <- map(c(names(kp_oid), q1_ent), \(ent) tibble(entity = ent,
    dISF_lesional_vs_healthy = fmt_e(est_p(ent, "states_all_visits", "AD_L_vs_HC")[["e"]], est_p(ent, "states_all_visits", "AD_L_vs_HC")[["p"]]),
    dISF_p = est_p(ent, "states_all_visits", "AD_L_vs_HC")[["p"]],
    serum_AD_vs_healthy = fmt_e(est_p(ent, "AD_vs_HC_in_study", "AD_vs_HC")[["e"]], est_p(ent, "AD_vs_HC_in_study", "AD_vs_HC")[["p"]]),
    serum_p = est_p(ent, "AD_vs_HC_in_study", "AD_vs_HC")[["p"]])) |> bind_rows() |> distinct(entity, .keep_all = TRUE)
q5_auc <- q4_auc |> filter(str_detect(predictor, "before relapse")) |>
  select(entity, predictor, matrix, AUC, ci_low, ci_high)

# ---- answers ------------------------------------------------------------------------------------------------------------
ans <- list()
add <- \(q, item, verdict, evidence) ans[[length(ans) + 1]] <<- tibble(question = q, item = item, verdict = verdict, evidence = evidence)
Q <- c(Q1 = "Q1 Are mast cell markers elevated in AD, or only in relapse vs non-relapse?",
       Q2 = "Q2 Is TNFRSF9 (CD137) or its ligand TNFSF9 a marker for mast cells in AD?",
       Q3 = "Q3 Do TNFRSF9 / TNFSF9 correlate with AD relapse?",
       Q4 = "Q4 Are TNFRSF9 / TNFSF9 a marker or predictor of AD relapse?",
       Q5 = "Q5 Is dISF superior to serum?")
for (i in seq_len(nrow(q1_verdict))) {
  r <- q1_verdict[i, ]; d <- q1 |> filter(entity == r$entity)
  add(Q[["Q1"]], r$entity, r$answer,
      paste(sprintf("%s %s", d$label, fmt_e(d$estimate, d$p)), collapse = "; "))
}
for (kp in names(kp_oid)) {
  s <- q2 |> filter(protein == kp, partner == "Mast cell score", subset == "all AD dISF")
  if (!nrow(s)) s <- q2 |> filter(protein == kp, subset == "all AD dISF") |> slice_min(p_rho, n = 1)
  pos <- nrow(s) && ((!is.na(s$p_within) && s$p_within < 0.05 && s$r_within > 0) || (!is.na(s$p_rho) && s$p_rho < 0.05 && s$rho > 0))
  neg <- nrow(s) && ((!is.na(s$p_within) && s$p_within < 0.05 && s$r_within < 0) || (!is.na(s$p_rho) && s$p_rho < 0.05 && s$rho < 0))
  mk <- q2 |> filter(protein == kp, partner != "Mast cell score", subset == "all AD dISF", p_rho < 0.05, rho > 0) |> pull(partner)
  add(Q[["Q2"]], kp, case_when(pos ~ "yes - co-varies with the mast cell signal in AD dISF",
                               neg ~ "inverse relation to the mast cell signal", TRUE ~ "no association with the mast cell signal shown"),
      sprintf("vs %s in AD dISF: within patients r = %.2f (p = %s); across samples rho = %.2f (p = %s); positively correlated single markers: %s",
              s$partner[1], s$r_within[1], fmt_p(s$p_within[1]), s$rho[1], fmt_p(s$p_rho[1]), if (length(mk)) paste(mk, collapse = ", ") else "none"))
}
for (kp in names(kp_oid)) {
  d <- q3 |> filter(entity == kp)
  sig <- d |> filter(p < 0.05)
  add(Q[["Q3"]], kp, if (nrow(sig)) sprintf("nominal association in %d of %d relapse tests (exploratory)", nrow(sig), nrow(d)) else
        sprintf("no association in %d relapse tests", nrow(d)),
      paste(sprintf("%s %s: %s", d$matrix, d$model, fmt_e(d$estimate, d$p)), collapse = "; "))
}
for (kp in names(kp_oid)) {
  m <- q4_marker |> filter(entity == kp)
  a <- q4_auc |> filter(entity == kp, !is.na(AUC))
  best <- a |> slice_max(abs(AUC - 0.5), n = 1, with_ties = FALSE)
  add(Q[["Q4"]], paste(kp, "- marker of lesion activity"),
      if (any(m$p < 0.05)) "yes - changes with lesion activity" else "no change with lesion activity shown",
      paste(sprintf("%s %s %s", m$matrix, m$contrast, fmt_e(m$estimate, m$p)), collapse = "; "))
  hits <- a$verdict[str_detect(a$verdict, "^predicts|^possible")]
  add(Q[["Q4"]], paste(kp, "- predictor of relapse"),
      if (length(hits)) sprintf("%s (in %d of %d predictor sets tested)", paste(unique(hits), collapse = "; "), length(hits), nrow(a))
      else sprintf("no predictive value shown (%d predictor sets tested)", nrow(a)),
      paste(sprintf("%s: AUC %.2f (95%% CI %.2f-%.2f, n = %d vs %d)", a$predictor, a$AUC, a$ci_low, a$ci_high, a$n_relapse, a$n_non), collapse = "; "))
}
det_txt <- paste(sprintf("%s %d of %d proteins detectable", q5_det$matrix, q5_det$proteins_detected, q5_det$proteins_measured), collapse = "; ")
add(Q[["Q5"]], "measurable proteome", det_txt, det_txt)
if (!is.null(q5_overlap)) for (i in which(q5_overlap$tier == "FDR")) {
  o <- q5_overlap[i, ]
  add(Q[["Q5"]], sprintf("%s (%s, all visits, FDR)", o$question, o$isf_site),
      if (o$dISF_significant > o$serum_significant) "more signal in dISF" else if (o$dISF_significant < o$serum_significant) "more signal in serum" else "equal",
      sprintf("dISF %d vs serum %d significant proteins; both %d; dISF only %d (+%d not measurable in serum); serum only %d",
              o$dISF_significant, o$serum_significant, o$both_same + o$both_opposite, o$dISF_only, o$dISF_only_not_in_serum, o$serum_only))
}
nd <- sum(q5_effects$dISF_p < 0.05, na.rm = TRUE); ns <- sum(q5_effects$serum_p < 0.05, na.rm = TRUE)
add(Q[["Q5"]], "key proteins and mast cell markers: AD vs healthy", if (nd > ns) "more signal in dISF" else if (nd < ns) "more signal in serum" else "equal",
    sprintf("p < 0.05 in dISF (lesional vs healthy skin) for %d of %d, in serum (AD vs healthy) for %d of %d", nd, nrow(q5_effects), ns, nrow(q5_effects)))
if (nrow(q5_auc)) {
  cmp <- q5_auc |> group_by(entity) |> summarise(best_dISF = suppressWarnings(max(abs(AUC[matrix == "dISF"] - 0.5), na.rm = TRUE)),
                                                 best_serum = suppressWarnings(max(abs(AUC[matrix == "serum"] - 0.5), na.rm = TRUE)))
  add(Q[["Q5"]], "relapse prediction, same MicroAD patients", sprintf("dISF separates better for %d of %d, serum for %d",
      sum(cmp$best_dISF > cmp$best_serum, na.rm = TRUE), nrow(cmp), sum(cmp$best_dISF < cmp$best_serum, na.rm = TRUE)),
      paste(sprintf("%s: |AUC - 0.5| dISF %.2f vs serum %.2f", cmp$entity, cmp$best_dISF, cmp$best_serum), collapse = "; "))
}
answers <- bind_rows(ans)
save_csv(answers, cfg, "key_questions", "answers.csv")
save_csv(q4_auc, cfg, "key_questions", "Q4_auc.csv")
save_csv(q1, cfg, "key_questions", "Q1_tests.csv")
writexl::write_xlsx(Filter(\(x) !is.null(x) && nrow(x), list(
  answers = answers, Q1_mast_markers = q1, Q1_verdicts = q1_verdict, Q2_correlations = q2, Q3_relapse = q3,
  Q4_marker = q4_marker, Q4_prediction_auc = q4_auc, Q5_detection = q5_det, Q5_overlap = q5_overlap,
  Q5_key_protein_effects = q5_effects, all_tests = tests)), out_path(cfg, "key_questions", "key_questions.xlsx"))
for (q in unique(answers$question)) {
  message("\n", q)
  a <- answers |> filter(question == q)
  for (i in seq_len(nrow(a))) message(sprintf("  - %s: %s", a$item[i], a$verdict[i]))
}
msg("Key questions: %s", file.path(cfg$paths$output, "key_questions"))
