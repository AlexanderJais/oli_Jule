# 20 - Do serum and dISF carry the same or different biological signatures?
# MicroAD only: RELAD / RELAD2 (and LEIP) serum is not used anywhere in this step.
#   1. Effect-size concordance: log2FC(serum) vs log2FC(dISF) per comparison (AD vs healthy with the
#      lesion site or the non-lesional site; relapse vs non-relapse before relapse), per visit and
#      pooled (models of step 14). Spearman rho, OLS and major-axis slope, % same sign; all proteins
#      measured in both, and restricted to proteins significant in either. Proteins below LOD in
#      MicroAD serum are left out and reported separately (dISF-only by detection).
#   2. Time-resolved model per compartment: limma + duplicateCorrelation (subject), cells group x
#      visit, moderated F-tests for group x visit. Healthy controls have one visit, so "AD vs
#      healthy x visit" = does the AD-vs-healthy difference change over visits. dISF-significant
#      proteins are classified by temporal profile (high at V1 and resolving / persistent /
#      late-rising) and compared with the serum profile of the same protein.
#   3. Compartment x group (x visit) interaction on paired samples (same subject and visit): dISF and
#      serum stacked, variance weights per compartment, subject blocking.
#   4. Within-subject paired correlation of serum vs dISF levels: tracks serum (systemic spill-over
#      likely) or not (local production more likely).
#   5. Signature sets (dISF-only, serum-only, shared-concordant, shared-discordant; FDR, nominal as
#      sensitivity) with protein lists, Reactome / GO:BP over-representation against the assayed
#      panel, and tissue origin from the Human Protein Atlas.
#   6. Relapse, predictive: dISF (and serum) at the visit before the relapse vs non-relapsers at the
#      same visits.
# Out: output/signatures/  signatures.xlsx, answers.csv, CSVs and figures

source("R/utils.R")
source("R/models.R")
source("R/design.R")
source("R/correlation.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
cats  <- read_step(cfg, "data", "serum_vs_disf_results.rds", step = "scripts/14_serum_vs_disf.R")
clear_outputs(cfg, "signatures")
fdr <- cfg$stats$fdr; min_det <- cfg$qc$min_detect_frac; vmin <- cfg$stats$visit_min_subjects %||% 5
sg <- cfg$signatures %||% list()
assay_map <- clean |> distinct(OlinkID, Assay)
fmt_p <- \(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
visit_levels <- c(paste0("V", 1:12), "all visits")

# ---- serum detectability in MicroAD serum ---------------------------------------------------------------------
ser_det <- clean |> filter(matrix == "Serum", cohort == "MicroAD", group %in% c("AD", "HC")) |>
  group_by(OlinkID, group) |> summarise(frac = mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = group, values_from = frac, names_prefix = "serum_frac_above_LOD_") |>
  mutate(serum_detected_MicroAD = pmax(coalesce(serum_frac_above_LOD_AD, 0), coalesce(serum_frac_above_LOD_HC, 0)) >= min_det)
serum_ok <- intersect(rownames(wide$Serum), ser_det$OlinkID[ser_det$serum_detected_MicroAD])
isf_ok   <- rownames(wide$ISF)
both     <- intersect(isf_ok, serum_ok)
serum_below <- setdiff(isf_ok, serum_ok)          # analysed in dISF, below LOD (or not measurable) in MicroAD serum
msg("dISF analysed: %d; measured in both (serum above LOD in MicroAD): %d; dISF only (serum below LOD): %d",
    length(isf_ok), length(both), length(serum_below))

cats <- cats |> mutate(visit = factor(as.character(visit), levels = visit_levels) |> droplevels(),
                       serum_measured = OlinkID %in% serum_ok, in_both = OlinkID %in% both)

# ---- 1. effect-size concordance -------------------------------------------------------------------------------
conc_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 5) return(tibble(n = sum(ok), spearman_rho = NA_real_, spearman_p = NA_real_, slope_ols = NA_real_,
                                 slope_major_axis = NA_real_, pct_same_sign = NA_real_))
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = sum(ok), spearman_rho = unname(s$estimate), spearman_p = s$p.value,
         slope_ols = unname(coef(lm(y[ok] ~ x[ok]))[2]), slope_major_axis = ma_slope(x[ok], y[ok]),
         pct_same_sign = round(100 * mean(sign(x[ok]) == sign(y[ok])), 1))
}
fdr_tier <- cats |> filter(tier == "FDR", in_both)
nom_tier <- cats |> filter(tier == "nominal", in_both)
conc <- bind_rows(
  fdr_tier |> group_by(question, isf_site, visit) |> group_modify(\(d, k) conc_stats(d$isf_logFC, d$ser_logFC)) |>
    mutate(proteins = "all measured in both"),
  fdr_tier |> filter(isf_sig %in% TRUE | ser_sig %in% TRUE) |> group_by(question, isf_site, visit) |>
    group_modify(\(d, k) conc_stats(d$isf_logFC, d$ser_logFC)) |> mutate(proteins = sprintf("FDR < %g in either", fdr)),
  nom_tier |> filter(isf_sig %in% TRUE | ser_sig %in% TRUE) |> group_by(question, isf_site, visit) |>
    group_modify(\(d, k) conc_stats(d$isf_logFC, d$ser_logFC)) |> mutate(proteins = "p < 0.05 in either (sensitivity)")) |>
  ungroup() |> relocate(proteins, .after = visit) |> arrange(question, isf_site, proteins, visit)
save_csv(conc, cfg, "signatures", "1_effect_concordance.csv")

below_sig <- cats |> filter(OlinkID %in% serum_below, isf_sig %in% TRUE) |>
  select(tier, question, isf_site, visit, Assay, OlinkID, isf_logFC, isf_p, isf_fdr) |>
  left_join(ser_det, by = "OlinkID")
below_summary <- cats |> filter(isf_sig %in% TRUE, !is.na(isf_logFC)) |>
  group_by(tier, question, isf_site, visit) |>
  summarise(dISF_significant = n(), of_which_serum_measured = sum(in_both),
            serum_measured_but_not_significant = sum(in_both & !(ser_sig %in% TRUE)),
            serum_below_LOD = sum(OlinkID %in% serum_below), .groups = "drop") |>
  mutate(pct_dISF_only_due_to_detection = round(100 * serum_below_LOD / pmax(dISF_significant - (of_which_serum_measured - serum_measured_but_not_significant), 1), 1))

for (qq in unique(fdr_tier$question)) {
  d <- fdr_tier |> filter(question == qq, !is.na(isf_logFC), !is.na(ser_logFC)) |>
    mutate(sig = case_when(isf_sig %in% TRUE & ser_sig %in% TRUE ~ "both", isf_sig %in% TRUE ~ "dISF", ser_sig %in% TRUE ~ "serum", TRUE ~ "neither"))
  lab <- conc |> filter(question == qq, proteins == "all measured in both") |>
    mutate(label = sprintf("rho %.2f\nslope %.2f (MA %.2f)\nsame sign %.0f%%", spearman_rho, slope_ols, slope_major_axis, pct_same_sign))
  p <- ggplot(d, aes(isf_logFC, ser_logFC)) + geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(data = \(x) filter(x, sig == "neither"), colour = "grey75", size = 0.5) +
    geom_point(data = \(x) filter(x, sig != "neither"), aes(colour = sig), size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.1, size = 2.3) +
    scale_colour_manual(values = c(both = "purple3", dISF = "firebrick", serum = "steelblue")) +
    facet_grid(isf_site ~ visit) +
    labs(title = sprintf("%s: effect in serum vs effect in dISF (proteins measured in both; MicroAD)", qq),
         x = "log2FC dISF", y = "log2FC serum", colour = sprintf("FDR < %g", fdr)) + theme(legend.position = "bottom")
  save_plot(p, cfg, "signatures", sprintf("1_concordance_%s.png", str_replace_all(qq, "[^A-Za-z]+", "_")),
            width = 3 + 2.4 * n_distinct(d$visit), height = 7)
}

# ---- 2. time-resolved model per compartment ------------------------------------------------------------------
isf_all <- isf_design(meta) |> filter(SampleID %in% colnames(wide$ISF))
visits_ok <- isf_all |> filter(group == "AD", !is.na(visit_num)) |> distinct(SubjectID, visit_num) |> count(visit_num) |>
  filter(n >= vmin) |> pull(visit_num) |> sort()
if (length(visits_ok) < 2) stop("Fewer than 2 visits with >= ", vmin, " patients - no time-resolved model possible.")
V <- paste0("V", visits_ok); v1 <- V[1]; later <- V[-1]

isf_t <- isf_all |> filter(group == "HC" | (group == "AD" & site %in% c("L", "NL") & visit_num %in% visits_ok)) |>
  mutate(cell = if_else(group == "HC", "HC", paste0(site, "_", visit)))
ct_isf <- c(setNames(sprintf("cellL_%s - cellHC", V), paste0("L_", V, "_vs_HC")),
            setNames(sprintf("cellNL_%s - cellHC", V), paste0("NL_", V, "_vs_HC")),
            setNames(sprintf("cellL_%s - cellL_%s", later, v1), paste0("L_", later, "_vs_", v1)),
            setNames(sprintf("cellNL_%s - cellNL_%s", later, v1), paste0("NL_", later, "_vs_", v1)),
            setNames(sprintf("(cellL_%s - cellNL_%s) - (cellL_%s - cellNL_%s)", later, later, v1, v1), paste0("siteXvisit_", later)))
ft_isf <- list(`dISF lesion site vs healthy, any visit` = paste0("L_", V, "_vs_HC"),
               `dISF non-lesional vs healthy, any visit` = paste0("NL_", V, "_vs_HC"),
               `dISF (AD lesion site vs healthy) x visit` = paste0("L_", later, "_vs_", v1),
               `dISF (AD non-lesional vs healthy) x visit` = paste0("NL_", later, "_vs_", v1),
               `dISF site (lesion - non-lesional) x visit` = paste0("siteXvisit_", later))
m_isf <- limma_block(wide$ISF, isf_t, ~ 0 + cell + plate, ct_isf, ft_isf, label = "dISF group x visit")

ser_all <- serum_design(meta) |> filter(cohort == "MicroAD", SampleID %in% colnames(wide$Serum), group %in% c("AD", "HC"))
ser_t <- ser_all |> filter(group == "HC" | visit_num %in% visits_ok) |>
  mutate(cell = if_else(group == "HC", "HC", paste0("AD_", visit)))
ct_ser <- c(setNames(sprintf("cellAD_%s - cellHC", V), paste0("AD_", V, "_vs_HC")),
            setNames(sprintf("cellAD_%s - cellAD_%s", later, v1), paste0("AD_", later, "_vs_", v1)))
ft_ser <- list(`serum AD vs healthy, any visit` = paste0("AD_", V, "_vs_HC"),
               `serum (AD vs healthy) x visit` = paste0("AD_", later, "_vs_", v1))
m_ser <- limma_block(wide$Serum[serum_ok, , drop = FALSE], ser_t, ~ 0 + cell + plate, ct_ser, ft_ser, label = "serum group x visit")

# relapse x visit (AD only; all visits - the relapse visit itself is lesional again, so exploratory)
relapse_visit_model <- function(expr, info, label) {
  info <- info |> filter(!is.na(relapse2), visit_num %in% visits_ok) |> mutate(cell = paste0(relapse2, "_", visit))
  vv <- info |> count(visit, relapse2) |> filter(n >= 2) |> count(visit) |> filter(n == 2) |> pull(visit)
  vv <- intersect(V, vv)
  if (length(vv) < 2) return(NULL)
  info <- info |> filter(visit %in% vv)
  ct <- c(setNames(sprintf("cellrelapse_%s - cellnon_relapse_%s", vv, vv), paste0("rel_", vv)),
          setNames(sprintf("(cellrelapse_%s - cellnon_relapse_%s) - (cellrelapse_%s - cellnon_relapse_%s)", vv[-1], vv[-1], vv[1], vv[1]),
                   paste0("relXvisit_", vv[-1])))
  ft <- setNames(list(paste0("rel_", vv), paste0("relXvisit_", vv[-1])),
                 c(paste(label, "relapse vs non-relapse, any visit"), paste(label, "relapse x visit")))
  limma_block(expr, info, ~ 0 + cell + plate, ct, ft, min_group_n = 2, label = paste(label, "relapse x visit"))
}
m_rel <- list(relapse_visit_model(wide$ISF, isf_all |> filter(group == "AD", site == "L"), "dISF lesion site"),
              relapse_visit_model(wide$ISF, isf_all |> filter(group == "AD", site == "NL"), "dISF non-lesional"),
              relapse_visit_model(wide$Serum[serum_ok, , drop = FALSE], ser_all |> filter(group == "AD"), "serum"))
time_models <- c(list(m_isf, m_ser), m_rel) |> Filter(f = Negate(is.null))
ftests <- map(time_models, "ftests") |> bind_rows() |> left_join(assay_map, by = "OlinkID") |> relocate(Assay, .after = OlinkID)
time_contrasts <- map(time_models, "contrasts") |> bind_rows() |> left_join(assay_map, by = "OlinkID") |> relocate(Assay, .after = OlinkID)
ftest_summary <- ftests |> group_by(model, ftest, n_samples, n_subjects, block_correlation) |>
  summarise(proteins = n(), significant_FDR = sum(adj.P.Val < fdr), nominal_p05 = sum(P.Value < 0.05),
            top_10 = paste(head(Assay[order(P.Value)], 10), collapse = ", "), .groups = "drop")
print(ftest_summary |> select(ftest, significant_FDR, nominal_p05) |> as.data.frame())
save_csv(ftests, cfg, "signatures", "2_time_models_Ftests.csv")

# temporal profiles of the dISF-significant proteins, and the serum profile of the same proteins
classify_profile <- function(e) {                      # e: effects vs healthy in visit order
  if (all(is.na(e))) return(c(direction = NA, profile = NA))
  s <- sign(e[which.max(abs(e))]); z <- s * e; peak <- max(z, na.rm = TRUE); first <- z[1]
  last <- mean(tail(z[!is.na(z)], 2))
  prof <- if (is.na(first) || first < 0.5 * peak) "late-rising" else if (last < 0.5 * first) "high at V1, resolving" else "persistent"
  c(direction = if (s > 0) "up in AD" else "down in AD", profile = prof)
}
prof_of <- function(prefix, model_ct) {
  time_contrasts |> filter(model == model_ct, str_detect(contrast, paste0("^", prefix, "_V[0-9]+_vs_HC$"))) |>
    mutate(visit = str_extract(contrast, "V[0-9]+")) |> select(OlinkID, Assay, visit, logFC)
}
profile_tab <- map(c("L", "NL"), \(sk) {
  site_lab <- if (sk == "L") "lesion site" else "non-lesional"
  ft_name <- if (sk == "L") "dISF lesion site vs healthy, any visit" else "dISF non-lesional vs healthy, any visit"
  sig <- ftests |> filter(ftest == ft_name, adj.P.Val < fdr) |> pull(OlinkID)
  if (!length(sig)) return(NULL)
  pi <- prof_of(sk, "dISF group x visit") |> filter(OlinkID %in% sig)
  ps <- prof_of("AD", "serum group x visit") |> filter(OlinkID %in% sig)
  map(sig, \(o) {
    ei <- pi |> filter(OlinkID == o) |> arrange(match(visit, V)); es <- ps |> filter(OlinkID == o) |> arrange(match(visit, V))
    ci <- classify_profile(ei$logFC)
    cs <- if (nrow(es)) classify_profile(es$logFC) else c(direction = NA, profile = NA)
    common <- intersect(ei$visit, es$visit)
    x <- ei$logFC[match(common, ei$visit)]; y <- es$logFC[match(common, es$visit)]
    tibble(OlinkID = o, site = site_lab, dISF_direction = ci[["direction"]], dISF_profile = ci[["profile"]],
           dISF_effects = paste(sprintf("%s %+.2f", ei$visit, ei$logFC), collapse = "; "),
           serum_measured = o %in% serum_ok,
           serum_profile = if (nrow(es)) paste(cs[["direction"]], cs[["profile"]], sep = ": ") else NA_character_,
           serum_effects = if (nrow(es)) paste(sprintf("%s %+.2f", es$visit, es$logFC), collapse = "; ") else NA_character_,
           profile_r = if (length(common) >= 3) suppressWarnings(cor(x, y)) else NA_real_,
           serum_amplitude = if (length(common) >= 2 && sum(x^2) > 0) sum(x * y) / sum(x^2) else NA_real_)
  }) |> bind_rows()
}) |> bind_rows()
if (nrow(profile_tab)) {
  serum_any <- ftests |> filter(ftest == "serum AD vs healthy, any visit") |> select(OlinkID, serum_any_visit_fdr = adj.P.Val)
  profile_tab <- profile_tab |> left_join(assay_map, by = "OlinkID") |> left_join(serum_any, by = "OlinkID") |>
    relocate(Assay, .after = OlinkID) |>
    mutate(serum_same_profile = serum_profile == paste(dISF_direction, dISF_profile, sep = ": "),
           serum_attenuated_copy = coalesce(profile_r > 0.5 & serum_amplitude > 0 & serum_amplitude < 1, FALSE))
  profile_summary <- profile_tab |> group_by(site, dISF_direction, dISF_profile) |>
    summarise(proteins = n(), serum_measured = sum(serum_measured),
              serum_significant_any_visit = sum(serum_any_visit_fdr < fdr, na.rm = TRUE),
              serum_same_profile = sum(serum_same_profile, na.rm = TRUE),
              median_profile_r = median(profile_r, na.rm = TRUE), median_serum_amplitude = median(serum_amplitude, na.rm = TRUE),
              serum_attenuated_copy = sum(serum_attenuated_copy), examples = paste(head(Assay, 12), collapse = ", "), .groups = "drop")
  print(as.data.frame(profile_summary |> select(-examples)), digits = 2)
  pl <- bind_rows(
    prof_of("L", "dISF group x visit") |> mutate(site = "lesion site"),
    prof_of("NL", "dISF group x visit") |> mutate(site = "non-lesional")) |>
    mutate(compartment = "dISF") |>
    inner_join(profile_tab |> select(OlinkID, site, dISF_direction, dISF_profile), by = c("OlinkID", "site"))
  pl <- bind_rows(pl, prof_of("AD", "serum group x visit") |> mutate(compartment = "serum") |>
                    inner_join(profile_tab |> select(OlinkID, site, dISF_direction, dISF_profile), by = "OlinkID", relationship = "many-to-many")) |>
    mutate(visit = factor(visit, levels = V), class = paste(dISF_direction, dISF_profile, sep = ": "))
  p <- ggplot(pl, aes(visit, logFC, group = OlinkID)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_line(alpha = 0.15, colour = "grey30") +
    stat_summary(aes(group = compartment, colour = compartment), fun = mean, geom = "line", linewidth = 1.2) +
    scale_colour_manual(values = c(dISF = "firebrick", serum = "steelblue")) +
    facet_grid(class ~ site + compartment, scales = "free_y") +
    labs(title = "Temporal profiles of dISF-significant proteins (AD vs healthy at each visit) and the same proteins in serum",
         subtitle = "thin lines = proteins, thick = mean; serum attenuated = same shape, smaller amplitude", x = NULL,
         y = "difference vs healthy (log2)", colour = NULL) + theme(legend.position = "bottom", strip.text.y = element_text(size = 7))
  save_plot(p, cfg, "signatures", "2_temporal_profiles.png", width = 12, height = 2.5 + 2 * n_distinct(pl$class))
} else profile_summary <- NULL

# ---- 3. compartment x group (x visit) on paired samples ------------------------------------------------------
paired_samples <- function(site_code) {
  i <- isf_all |> filter(group == "HC" | (group == "AD" & site == site_code)) |>
    select(SubjectID, visit, visit_num, grp = group, isf_id = SampleID, isf_plate = plate)
  s <- ser_all |> select(SubjectID, visit, ser_id = SampleID, ser_plate = plate)
  inner_join(i, s, by = c("SubjectID", "visit"))
}
stack_pairs <- function(pp) {
  expr <- cbind(wide$ISF[both, pp$isf_id, drop = FALSE], wide$Serum[both, pp$ser_id, drop = FALSE])
  colnames(expr) <- c(paste0("dISF:", pp$isf_id), paste0("serum:", pp$ser_id))
  info <- bind_rows(pp |> transmute(SampleID = paste0("dISF:", isf_id), SubjectID, visit, visit_num, grp, compartment = "dISF", plate = isf_plate),
                    pp |> transmute(SampleID = paste0("serum:", ser_id), SubjectID, visit, visit_num, grp, compartment = "serum", plate = ser_plate))
  list(expr = expr, info = info)
}
cxg <- map(c(L = "lesion site", NL = "non-lesional"), \(site_lab) {
  sk <- if (site_lab == "lesion site") "L" else "NL"
  pp <- paired_samples(sk) |> filter(grp == "HC" | visit_num %in% visits_ok)
  if (nrow(pp) < 8 || length(both) < 10) return(NULL)
  st <- stack_pairs(pp)
  info <- st$info |> mutate(cell = paste(compartment, grp, sep = "_"))
  pooled <- limma_block(st$expr, info, ~ 0 + cell + compartment:plate,
                        c(dISF_AD_vs_HC = "celldISF_AD - celldISF_HC", serum_AD_vs_HC = "cellserum_AD - cellserum_HC",
                          compartment_x_group = "(celldISF_AD - celldISF_HC) - (cellserum_AD - cellserum_HC)"),
                        weights_by = "compartment", label = paste("compartment x group,", site_lab))
  info_v <- st$info |> mutate(cell = if_else(grp == "HC", paste(compartment, "HC", sep = "_"), paste(compartment, "AD", visit, sep = "_")))
  vv <- intersect(V, unique(info_v$visit[info_v$grp == "AD"]))
  ctv <- setNames(sprintf("(celldISF_AD_%s - celldISF_HC) - (cellserum_AD_%s - cellserum_HC)", vv, vv), paste0("compXgroup_", vv))
  ctv3 <- setNames(sprintf("(celldISF_AD_%s - cellserum_AD_%s) - (celldISF_AD_%s - cellserum_AD_%s)", vv[-1], vv[-1], vv[1], vv[1]),
                   paste0("compXgroupXvisit_", vv[-1]))
  by_visit <- limma_block(st$expr, info_v, ~ 0 + cell + compartment:plate, c(ctv, ctv3),
                          list(`compartment x group, any visit` = names(ctv), `compartment x group x visit` = names(ctv3)),
                          weights_by = "compartment", label = paste("compartment x group x visit,", site_lab))
  list(pooled = pooled, by_visit = by_visit, site = site_lab, n_pairs = nrow(pp))
}) |> Filter(f = Negate(is.null))
cxg_tab <- map(cxg, \(r) {
  w <- r$pooled$contrasts |> select(OlinkID, contrast, logFC, P.Value, adj.P.Val) |>
    pivot_wider(names_from = contrast, values_from = c(logFC, P.Value, adj.P.Val))
  f <- r$by_visit$ftests |> select(OlinkID, ftest, P.Value, adj.P.Val) |>
    pivot_wider(names_from = ftest, values_from = c(P.Value, adj.P.Val), names_glue = "{ftest} {.value}")
  pv <- r$by_visit$contrasts |> filter(str_detect(contrast, "^compXgroup_V")) |>
    transmute(OlinkID, col = paste0("interaction_logFC_", str_remove(contrast, "compXgroup_")), logFC) |>
    pivot_wider(names_from = col, values_from = logFC)
  w |> left_join(f, by = "OlinkID") |> left_join(pv, by = "OlinkID") |> mutate(site = r$site, n_pairs = r$n_pairs, .before = 1)
}) |> bind_rows()
if (nrow(cxg_tab)) {
  cxg_tab <- cxg_tab |> left_join(assay_map, by = "OlinkID") |> relocate(Assay, .after = OlinkID) |>
    mutate(pattern = case_when(
      !coalesce(adj.P.Val_compartment_x_group < fdr, FALSE) ~ "same disease effect in both (no interaction)",
      sign(logFC_dISF_AD_vs_HC) != sign(logFC_serum_AD_vs_HC) & adj.P.Val_dISF_AD_vs_HC < fdr & adj.P.Val_serum_AD_vs_HC < fdr ~ "opposite direction",
      abs(logFC_dISF_AD_vs_HC) > abs(logFC_serum_AD_vs_HC) ~ "stronger in dISF",
      TRUE ~ "stronger in serum")) |>
    arrange(site, P.Value_compartment_x_group)
  save_csv(cxg_tab, cfg, "signatures", "3_compartment_x_group.csv")
  print(count(cxg_tab, site, pattern))
}

# ---- 4. within-subject paired correlation ------------------------------------------------------------------------
fast_rm <- function(x, y, subj) {         # repeated-measures correlation (= rmcorr): subject-centred values
  ok <- !is.na(x) & !is.na(y); keep <- ok & subj %in% names(which(table(subj[ok]) >= 2))
  if (sum(keep) < 5 || n_distinct(subj[keep]) < 2) return(c(NA, NA, NA))
  xc <- x[keep] - ave(x[keep], subj[keep]); yc <- y[keep] - ave(y[keep], subj[keep])
  r <- sum(xc * yc) / sqrt(sum(xc^2) * sum(yc^2)); dfree <- sum(keep) - n_distinct(subj[keep]) - 1
  c(r, 2 * pt(-abs(r * sqrt(dfree / max(1 - r^2, 1e-12))), dfree), dfree)
}
paired_corr <- map(c("L", "NL"), \(sk) {
  pp <- paired_samples(sk)
  if (nrow(pp) < 8) return(NULL)
  X <- wide$ISF[both, pp$isf_id, drop = FALSE]; Y <- wide$Serum[both, pp$ser_id, drop = FALSE]
  ad <- pp$grp == "AD"
  map(both, \(o) {
    x <- X[o, ]; y <- Y[o, ]
    s <- spearman_row(x, y)
    w <- fast_rm(x[ad], y[ad], pp$SubjectID[ad])
    m <- tibble(s = pp$SubjectID, x, y) |> group_by(s) |> summarise(x = mean(x, na.rm = TRUE), y = mean(y, na.rm = TRUE))
    b <- spearman_row(m$x, m$y)
    tibble(OlinkID = o, n_pairs = s$n, rho_paired = s$rho, p_paired = s$p, r_within_subject = w[1], p_within_subject = w[2],
           df_within = w[3], rho_between_subjects = b$rho, p_between_subjects = b$p)
  }) |> bind_rows() |> mutate(site = if (sk == "L") "lesion site" else "non-lesional", .before = 1)
}) |> bind_rows()
if (nrow(paired_corr)) {
  paired_corr <- paired_corr |> group_by(site) |>
    mutate(fdr_paired = p.adjust(p_paired, "BH"), fdr_within_subject = p.adjust(p_within_subject, "BH")) |> ungroup() |>
    mutate(tracks_serum = case_when(coalesce(fdr_paired < fdr & rho_paired > 0, FALSE) | coalesce(fdr_within_subject < fdr & r_within_subject > 0, FALSE) ~
                                      "tracks serum: systemic spill-over likely",
                                    TRUE ~ "does not track serum: local production more likely")) |>
    left_join(assay_map, by = "OlinkID") |> relocate(Assay, .after = OlinkID)
  enr_path <- file.path(cfg$paths$output, "matrix_comparison", "relative_enrichment.csv")
  if (file.exists(enr_path)) {
    enr <- read_csv(enr_path, show_col_types = FALSE) |>
      transmute(OlinkID, site = if_else(str_detect(model, "non-lesional"), "non-lesional", "lesion site"),
                relative_dISF_serum_log2 = rel_log2_isf_vs_serum, relative_enrichment = direction)
    paired_corr <- paired_corr |> left_join(enr, by = c("OlinkID", "site"))
  }
  save_csv(paired_corr, cfg, "signatures", "4_paired_correlation.csv")
  print(count(paired_corr, site, tracks_serum))
}
not_assessable <- tibble(OlinkID = serum_below) |> left_join(assay_map, by = "OlinkID") |> left_join(ser_det, by = "OlinkID") |>
  mutate(note = "below LOD (or not measurable) in MicroAD serum - serum tracking not assessable")

# ---- 5. signature sets, enrichment, tissue origin ------------------------------------------------------------------
sets <- cats |> filter(visit == "all visits", !is.na(isf_logFC) | !is.na(ser_logFC)) |>
  mutate(ser_sig_m = ser_sig %in% TRUE & serum_measured,
         set = case_when(isf_sig %in% TRUE & ser_sig_m & sign(isf_logFC) == sign(ser_logFC) ~ "shared-concordant",
                         isf_sig %in% TRUE & ser_sig_m ~ "shared-discordant",
                         isf_sig %in% TRUE & !serum_measured ~ "dISF-only (serum below LOD)",
                         isf_sig %in% TRUE ~ "dISF-only",
                         ser_sig_m & !(OlinkID %in% isf_ok) ~ "serum-only (not analysed in dISF)",
                         ser_sig_m ~ "serum-only",
                         TRUE ~ NA_character_))
set_lists <- sets |> filter(!is.na(set)) |>
  select(tier, question, isf_site, set, Assay, OlinkID, isf_logFC, isf_p, isf_fdr, ser_logFC, ser_p, ser_fdr) |>
  arrange(tier, question, isf_site, set, pmin(isf_p, ser_p, na.rm = TRUE))
set_counts <- sets |> filter(!is.na(set)) |> count(tier, question, isf_site, set) |>
  pivot_wider(names_from = set, values_from = n, values_fill = 0)
print(as.data.frame(set_counts))

# Human Protein Atlas: tissue / cell-type specificity (cached in paths$hpa; downloaded if missing)
hpa_path <- cfg$paths$hpa %||% "data/hpa_annotation.tsv"
hpa_url <- sg$hpa_url %||% "https://www.proteinatlas.org/api/search_download.php?search=&format=tsv&columns=g,up,rnats,rnatsm,rnascs,rnascsm,rnabcs,rnabcsm,secl&compress=no"
if (!file.exists(hpa_path)) {
  msg("Downloading Human Protein Atlas annotation to %s", hpa_path)
  dir.create(dirname(hpa_path), showWarnings = FALSE, recursive = TRUE)
  ok <- tryCatch(utils::download.file(hpa_url, hpa_path, quiet = TRUE, mode = "wb") == 0, error = \(e) FALSE, warning = \(w) FALSE)
  if (!ok) { unlink(hpa_path); msg("HPA download failed - tissue origin not annotated (place the file at %s)", hpa_path) }
}
hpa <- if (file.exists(hpa_path)) read_tsv(hpa_path, show_col_types = FALSE, progress = FALSE) else NULL
origin <- NULL
if (!is.null(hpa)) {
  names(hpa) <- str_replace_all(tolower(names(hpa)), "[^a-z]+", "_") |> str_remove("_$")
  col <- \(pat) { hit <- names(hpa)[str_detect(names(hpa), pat)]; if (length(hit)) hpa[[hit[1]]] else rep(NA_character_, nrow(hpa)) }
  hp <- tibble(gene = hpa$gene, tissue_specificity = col("^rna_tissue_specificity$"), tissue_nTPM = col("^rna_tissue_specific_ntpm"),
               celltype_specificity = col("^rna_single_cell_type_specificity$"), celltype_nCPM = col("^rna_single_cell_type_specific_n"),
               blood_cell_specificity = col("^rna_blood_cell_specificity$"), secretome = col("^secretome_location")) |>
    distinct(gene, .keep_all = TRUE)
  immune_cells <- "T-cells|B-cells|NK-cells|macrophage|monocyte|dendritic|cDC|pDC|granulocyte|neutrophil|eosinophil|basophil|mast cell|plasma cell|Langerhans|Kupffer|microglia|Hofbauer"
  hp <- hp |> mutate(
    skin = str_detect(coalesce(tissue_nTPM, ""), regex("skin", TRUE)) | str_detect(coalesce(celltype_nCPM, ""), regex("keratinocyte|melanocyte", TRUE)),
    immune = str_detect(coalesce(blood_cell_specificity, ""), "enriched|enhanced") & !str_detect(coalesce(blood_cell_specificity, ""), "Low|Not detected") |
      str_detect(coalesce(tissue_nTPM, ""), regex("lymphoid|bone marrow|spleen|thymus|tonsil|lymph node", TRUE)) |
      str_detect(coalesce(celltype_nCPM, ""), regex(immune_cells, TRUE)),
    liver = str_detect(coalesce(tissue_nTPM, ""), regex("liver", TRUE)) | str_detect(coalesce(celltype_nCPM, ""), regex("hepatocyte", TRUE)))
  origin <- assay_map |> mutate(gene = str_split(Assay, "_")) |> unnest(gene) |> inner_join(hp, by = "gene") |>
    group_by(OlinkID) |>
    summarise(skin = any(skin), immune = any(immune), liver = any(liver),
              tissue_specificity = paste(unique(na.omit(tissue_specificity)), collapse = "; "),
              tissue_nTPM = paste(unique(na.omit(tissue_nTPM)), collapse = "; "),
              celltype_nCPM = str_trunc(paste(unique(na.omit(celltype_nCPM)), collapse = "; "), 250),
              blood_cell_specificity = paste(unique(na.omit(blood_cell_specificity)), collapse = "; "),
              secretome = paste(unique(na.omit(secretome)), collapse = "; "), .groups = "drop") |>
    mutate(tissue_origin = pmap_chr(list(skin, immune, liver), \(s, i, l) {
      x <- c(if (s) "skin/keratinocyte", if (i) "immune", if (l) "liver"); if (length(x)) paste(x, collapse = " + ") else "other / not specific" }))
  set_lists <- set_lists |> left_join(origin |> select(OlinkID, tissue_origin, skin, immune, liver, tissue_nTPM, celltype_nCPM,
                                                       blood_cell_specificity, secretome), by = "OlinkID")
  bg <- sets |> distinct(tier, question, isf_site, OlinkID) |> left_join(origin |> select(OlinkID, skin, immune, liver), by = "OlinkID")
  origin_summary <- set_lists |> group_by(tier, question, isf_site, set) |>
    summarise(proteins = n(), annotated = sum(!is.na(tissue_origin)), skin = sum(skin, na.rm = TRUE),
              immune = sum(immune, na.rm = TRUE), liver = sum(liver, na.rm = TRUE), .groups = "drop") |>
    left_join(bg |> group_by(tier, question, isf_site) |>
                summarise(bg_n = n(), bg_skin = sum(skin, na.rm = TRUE), bg_immune = sum(immune, na.rm = TRUE),
                          bg_liver = sum(liver, na.rm = TRUE), .groups = "drop"), by = c("tier", "question", "isf_site")) |>
    mutate(across(c(skin, immune, liver), \(x) round(100 * x / pmax(proteins, 1), 1), .names = "pct_{.col}"),
           across(c(bg_skin, bg_immune, bg_liver), \(x) round(100 * x / pmax(bg_n, 1), 1), .names = "pct_{.col}"),
           p_skin = pmap_dbl(list(skin, proteins, bg_skin, bg_n), \(a, n, b, N) fisher.test(matrix(c(a, n - a, b - a, N - n - b + a), 2))$p.value),
           p_immune = pmap_dbl(list(immune, proteins, bg_immune, bg_n), \(a, n, b, N) fisher.test(matrix(c(a, n - a, b - a, N - n - b + a), 2))$p.value),
           p_liver = pmap_dbl(list(liver, proteins, bg_liver, bg_n), \(a, n, b, N) fisher.test(matrix(c(a, n - a, b - a, N - n - b + a), 2))$p.value))
  po <- set_lists |> filter(tier == "FDR", !is.na(tissue_origin)) |> count(question, isf_site, set, tissue_origin)
  if (nrow(po)) {
    p <- ggplot(po, aes(set, n, fill = tissue_origin)) + geom_col(position = "fill") +
      facet_grid(question ~ isf_site) + scale_y_continuous(labels = scales::percent) +
      labs(title = sprintf("Tissue origin (Human Protein Atlas) of the signature sets (FDR < %g, all visits)", fdr),
           x = NULL, y = "share of proteins", fill = NULL) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
    save_plot(p, cfg, "signatures", "5_tissue_origin.png", width = 12, height = 8)
  }
} else origin_summary <- NULL

# over-representation (Reactome, GO:BP) with the assayed panel as background
colls <- unlist(sg$collections %||% c("C2:CP:REACTOME", "C5:GO:BP"))
gsets <- tryCatch(map(colls, \(cl) {
  parts <- str_split_fixed(cl, ":", 2)
  g <- if (parts[2] == "") msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1])
       else msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = parts[2])
  split(g$gene_symbol, g$gs_name)
}) |> unlist(recursive = FALSE), error = \(e) { msg("Gene sets not available (%s) - enrichment skipped", conditionMessage(e)); NULL })
ora <- NULL
if (!is.null(gsets)) {
  genes_of <- \(a) unique(unlist(str_split(a, "_")))
  ora <- sets |> filter(!is.na(set)) |> group_by(tier, question, isf_site, set) |>
    group_modify(\(d, k) {
      universe <- genes_of(sets$Assay[sets$tier == k$tier & sets$question == k$question & sets$isf_site == k$isf_site])
      g <- genes_of(d$Assay)
      if (length(g) < 3) return(tibble())
      r <- suppressWarnings(fgsea::fora(gsets, genes = g, universe = universe, minSize = cfg$enrichment$min_size %||% 10,
                                        maxSize = cfg$enrichment$max_size %||% 500))
      as_tibble(r) |> filter(overlap >= 2) |> mutate(overlapGenes = map_chr(overlapGenes, paste, collapse = ";"),
                                                     set_genes = length(g), universe = length(universe))
    }) |> ungroup() |> arrange(tier, question, isf_site, set, pval)
  save_csv(ora, cfg, "signatures", "5_enrichment.csv")
}
save_csv(set_lists, cfg, "signatures", "5_signature_protein_lists.csv")

# ---- 6. relapse, predictive: the visit before the relapse ---------------------------------------------------------
pre_relapse_model <- function(expr, info, label) {
  info <- info |> filter(group == "AD", !is.na(relapse2), !is.na(visit_num))
  pre <- info |> filter(relapse2 == "relapse", !is.na(relapse_visit), visit_num < relapse_visit) |>
    group_by(SubjectID) |> filter(visit_num == max(visit_num)) |> ungroup()
  if (n_distinct(pre$SubjectID) < 2) return(NULL)
  non <- info |> filter(relapse2 == "non_relapse", visit_num %in% pre$visit_num)
  d <- bind_rows(pre, non) |> mutate(visitf = factor(visit))
  r <- limma_block(expr, d, ~ 0 + relapse2 + visitf + plate, c(pre_relapse_vs_non = "relapse2relapse - relapse2non_relapse"),
                   min_group_n = 2, label = label)
  if (is.null(r)) return(NULL)
  r$contrasts |> mutate(n_relapse = n_distinct(pre$SubjectID), n_non_relapse_subjects = n_distinct(non$SubjectID),
                        n_non_relapse_samples = nrow(non),
                        visits = paste(sort(unique(pre$visit)), collapse = ","))
}
isf_ad <- isf_all |> filter(group == "AD")
pre_res <- list(
  pre_relapse_model(wide$ISF, isf_ad |> filter(site == "L"), "dISF lesion site (ex-lesional), visit before relapse"),
  pre_relapse_model(wide$ISF, isf_ad |> filter(site == "NL"), "dISF non-lesional, visit before relapse"),
  pre_relapse_model(wide$Serum[serum_ok, , drop = FALSE], ser_all, "serum (MicroAD), visit before relapse"))
# xL - NL difference at the same visits (removes plate and day-to-day variation)
dp <- isf_ad |> filter(site %in% c("L", "NL")) |> select(SubjectID, visit, visit_num, site, SampleID, relapse2, relapse_visit, plate) |>
  pivot_wider(id_cols = c(SubjectID, visit, visit_num, relapse2, relapse_visit), names_from = site, values_from = c(SampleID, plate)) |>
  filter(!is.na(SampleID_L), !is.na(SampleID_NL)) |> mutate(SampleID = paste(SubjectID, visit, sep = "_"), group = "AD", plate = plate_L)
if (nrow(dp) >= 6) {
  de <- wide$ISF[, dp$SampleID_L, drop = FALSE] - wide$ISF[, dp$SampleID_NL, drop = FALSE]
  colnames(de) <- dp$SampleID
  pre_res <- c(pre_res, list(pre_relapse_model(de, dp, "dISF lesion site minus non-lesional, visit before relapse")))
}
pre_res <- bind_rows(Filter(Negate(is.null), pre_res))
if (nrow(pre_res)) {
  pre_res <- pre_res |> left_join(assay_map, by = "OlinkID") |> relocate(Assay, .after = OlinkID) |>
    mutate(significant = adj.P.Val < fdr) |> arrange(model, P.Value)
  save_csv(pre_res, cfg, "signatures", "6_pre_relapse.csv")
  cross_sec <- cats |> filter(tier == "FDR", question == "relapse vs non-relapse", visit != "all visits") |>
    group_by(visit, isf_site) |> summarise(dISF_FDR = sum(isf_sig %in% TRUE), serum_FDR = sum(ser_sig %in% TRUE), .groups = "drop")
  pre_summary <- pre_res |> group_by(model, visits, n_relapse, n_non_relapse_subjects, n_non_relapse_samples) |>
    summarise(significant_FDR = sum(significant), nominal_p05 = sum(P.Value < 0.05),
              top_10 = paste(head(Assay, 10), collapse = ", "), .groups = "drop")
  print(as.data.frame(pre_summary |> select(model, n_relapse, significant_FDR, nominal_p05)))
  p <- ggplot(pre_res, aes(logFC, -log10(P.Value), colour = significant)) + geom_point(size = 0.7) +
    geom_text(data = \(d) d |> group_by(model) |> slice_min(P.Value, n = 8), aes(label = Assay), size = 2.4, vjust = -0.6,
              colour = "black", check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) + facet_wrap(~str_wrap(model, 40)) +
    labs(title = "Visit before relapse: later relapsers vs non-relapsers at the same visits (exploratory)",
         x = "difference (log2)", y = "-log10 p") + theme(legend.position = "bottom")
  save_plot(p, cfg, "signatures", "6_pre_relapse_volcano.png", width = 11, height = 8)
} else { pre_summary <- NULL; cross_sec <- NULL }

# ---- answers --------------------------------------------------------------------------------------------------------
ans <- list()
add <- \(q, item, verdict, evidence) ans[[length(ans) + 1]] <<- tibble(question = q, item = item, verdict = verdict, evidence = evidence)
Q <- c("1 Effect-size concordance", "2 Time-resolved models", "3 Compartment x group", "4 Paired correlation",
       "5 Signature sets", "6 Relapse, predictive")
add(Q[1], "scope", "MicroAD only", sprintf("RELAD / RELAD2 / LEIP serum not used. Measured in both: %d proteins; dISF-analysed but below LOD in MicroAD serum: %d (reported separately).",
                                           length(both), length(serum_below)))
for (i in which(conc$visit == "all visits")) {
  r <- conc[i, ]
  add(Q[1], sprintf("%s, %s (%s)", r$question, r$isf_site, r$proteins),
      if (is.na(r$spearman_rho)) "not estimable" else if (r$spearman_p < 0.05 && r$spearman_rho > 0)
        sprintf("concordant (serum slope %.2f of dISF)", r$slope_major_axis) else "no concordance shown",
      sprintf("n %d; rho %.2f (p %s); OLS slope %.2f; major-axis slope %.2f; same sign %.0f%%", r$n, r$spearman_rho, fmt_p(r$spearman_p),
              r$slope_ols, r$slope_major_axis, r$pct_same_sign))
}
bs <- below_summary |> filter(tier == "FDR", visit == "all visits")
for (i in seq_len(nrow(bs)))
  add(Q[1], sprintf("dISF-only: absence of effect or of detection? %s, %s", bs$question[i], bs$isf_site[i]),
      sprintf("%d of %d dISF-significant proteins are below LOD in serum", bs$serum_below_LOD[i], bs$dISF_significant[i]),
      sprintf("measured in serum but not significant: %d; below LOD in serum: %d", bs$serum_measured_but_not_significant[i], bs$serum_below_LOD[i]))
for (i in seq_len(nrow(ftest_summary)))
  add(Q[2], ftest_summary$ftest[i], sprintf("%d proteins at FDR < %g", ftest_summary$significant_FDR[i], fdr),
      sprintf("nominal p < 0.05: %d; top: %s", ftest_summary$nominal_p05[i], ftest_summary$top_10[i]))
if (!is.null(profile_summary)) for (i in seq_len(nrow(profile_summary))) {
  r <- profile_summary[i, ]
  add(Q[2], sprintf("%s: %s, %s", r$site, r$dISF_direction, r$dISF_profile), sprintf("%d proteins; serum attenuated copy for %d of %d measured", r$proteins, r$serum_attenuated_copy, r$serum_measured),
      sprintf("median profile correlation dISF-serum %.2f; median serum amplitude %.2f x dISF; same serum profile class %d; e.g. %s",
              r$median_profile_r, r$median_serum_amplitude, r$serum_same_profile, r$examples))
}
if (nrow(cxg_tab)) for (st in unique(cxg_tab$site)) {
  d <- cxg_tab |> filter(site == st)
  add(Q[3], st, sprintf("%d proteins change with disease differently in dISF than in serum (FDR < %g)", sum(d$adj.P.Val_compartment_x_group < fdr, na.rm = TRUE), fdr),
      paste(sprintf("%s: %d", names(table(d$pattern)), table(d$pattern)), collapse = "; "))
}
if (nrow(paired_corr)) for (st in unique(paired_corr$site)) {
  d <- paired_corr |> filter(site == st)
  add(Q[4], st, sprintf("%d of %d proteins track serum", sum(str_detect(d$tracks_serum, "^tracks")), nrow(d)),
      sprintf("tracking (systemic): %s ...; %d dISF proteins below LOD in serum are not assessable",
              paste(head(d$Assay[order(d$p_paired)][str_detect(d$tracks_serum[order(d$p_paired)], "^tracks")], 15), collapse = ", "),
              length(serum_below)))
}
sc <- set_counts |> filter(tier == "FDR")
for (i in seq_len(nrow(sc))) {
  r <- sc[i, ]; cols <- setdiff(names(r), c("tier", "question", "isf_site"))
  add(Q[5], sprintf("%s, %s (FDR, all visits)", r$question, r$isf_site), "protein counts per set",
      paste(sprintf("%s %d", cols, unlist(r[cols])), collapse = "; "))
}
if (!is.null(ora)) {
  top_ora <- ora |> filter(tier == "FDR", padj < fdr) |> group_by(question, isf_site, set) |> slice_min(padj, n = 3, with_ties = FALSE) |> ungroup()
  add(Q[5], "enrichment (background = assayed panel)", sprintf("%d set x term hits at FDR < %g", sum(ora$tier == "FDR" & ora$padj < fdr, na.rm = TRUE), fdr),
      if (nrow(top_ora)) paste(sprintf("%s / %s / %s: %s", top_ora$question, top_ora$isf_site, top_ora$set, top_ora$pathway), collapse = "; ") else "none")
}
if (!is.null(pre_summary)) for (i in seq_len(nrow(pre_summary))) {
  r <- pre_summary[i, ]
  add(Q[6], r$model, sprintf("%d proteins at FDR < %g", r$significant_FDR, fdr),
      sprintf("%d relapsers (visits %s) vs %d non-relapsers (%d samples, same visits); nominal p < 0.05: %d; top: %s",
              r$n_relapse, r$visits, r$n_non_relapse_subjects, r$n_non_relapse_samples, r$nominal_p05, r$top_10))
}
answers <- bind_rows(ans)
save_csv(answers, cfg, "signatures", "answers.csv")

methods <- tibble(item = Q, method = c(
  "Per-visit and pooled models of step 14 (MicroAD; dISF lesion site or non-lesional vs healthy skin, serum AD vs healthy; relapse vs non-relapse on visits before the relapse). Proteins analysed in dISF and above LOD in >= 50% of AD or healthy MicroAD serum. OLS slope = serum on dISF; major-axis slope is symmetric (not attenuated by noise in the dISF effect).",
  "limma + duplicateCorrelation (subject), plate adjusted. dISF cells: healthy, lesion site x visit, non-lesional x visit (visits with >= stats$visit_min_subjects patients). Healthy controls have one visit, so group x visit = does the AD-vs-healthy difference change over visits (F-test of Vk - V1). Profile classes use the AD-vs-healthy effects of the F-significant proteins: peak at V1 and last two visits < 50% of V1 = resolving; >= 50% = persistent; V1 < 50% of the peak = late-rising. Serum: same classification, profile correlation and amplitude (least-squares scale of the serum profile to the dISF profile; 0 < amplitude < 1 with r > 0.5 = attenuated copy).",
  "Paired subject-visits with dISF and serum; values stacked, cells compartment x group (x visit), plate within compartment, variance weights per compartment (arrayWeights), subject blocking. Interaction = (dISF AD - dISF healthy) - (serum AD - serum healthy).",
  "Per protein measured in both: Spearman over all paired samples, repeated-measures correlation within AD patients (subject-centred), Spearman of subject means; BH within site.",
  "Pooled (all visits) results of step 14. Sets at FDR < 0.05 (tier FDR) and p < 0.05 (tier nominal, sensitivity). Over-representation: fgsea::fora, MSigDB Reactome and GO:BP, background = genes of all proteins analysed for that comparison. Tissue origin: Human Protein Atlas RNA tissue / single-cell / blood-cell specificity (skin incl. keratinocytes, immune, liver); Fisher test against the same background.",
  "Relapsers: the last visit before the relapse visit (first return of lesional skin). Non-relapsers: all their samples at the same visit numbers. limma, relapse + visit + plate, subject blocking; also on the lesion-site minus non-lesional difference and in MicroAD serum. 4 relapsers - exploratory."))
writexl::write_xlsx(Filter(\(x) !is.null(x) && nrow(x), list(
  answers = answers, methods = methods,
  `1_concordance` = conc, `1_dISF_sig_serum_below_LOD` = below_sig, `1_detection_summary` = below_summary,
  `2_Ftest_summary` = ftest_summary, `2_Ftests` = ftests, `2_profiles` = profile_tab, `2_profile_summary` = profile_summary,
  `3_compartment_x_group` = cxg_tab, `4_paired_correlation` = paired_corr, `4_not_assessable` = not_assessable,
  `5_set_counts` = set_counts, `5_protein_lists` = set_lists, `5_enrichment` = ora, `5_tissue_origin` = origin_summary,
  `6_pre_relapse_summary` = pre_summary, `6_pre_relapse` = pre_res, `6_cross_sectional_per_visit` = cross_sec,
  serum_detection_MicroAD = ser_det |> left_join(assay_map, by = "OlinkID"))),
  out_path(cfg, "signatures", "signatures.xlsx"))
for (q in unique(answers$question)) {
  message("\n", q)
  a <- answers |> filter(question == q)
  for (i in seq_len(nrow(a))) message(sprintf("  - %s: %s", a$item[i], a$verdict[i]))
}
msg("Serum vs dISF signatures: %s", file.path(cfg$paths$output, "signatures"))
