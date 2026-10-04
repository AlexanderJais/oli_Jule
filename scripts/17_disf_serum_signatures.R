# 17 - Do serum and dISF carry the same or different biological signatures?
# MicroAD only: RELAD, RELAD2 and LEIP serum are NOT used in this step.
# "Measured in both" = passes the dISF detection filter AND is above LOD in >= qc$min_detect_frac of at least one
# MicroAD serum group (AD or healthy). dISF proteins not detected in MicroAD serum are reported separately.
#   1  Effect-size concordance: log2FC(serum) vs log2FC(dISF) per comparison, per visit and pooled:
#      Spearman rho, slope, % same sign - all shared proteins, and those FDR-significant (or p < 0.05) in either.
#      Comparisons: AD vs healthy (dISF lesion site / non-lesional skin vs healthy skin; serum AD vs healthy) and
#      relapse vs non-relapse before the relapse (dISF lesion site; serum; V1 lesional for all, later visits only
#      while the lesion is cleared; pooled = cleared visits, adjusted for weeks since V1).
#   2  Time-resolved model per compartment: group x visit (limma + duplicateCorrelation). Healthy controls have one
#      visit, so the group x visit test asks whether the AD - healthy difference changes over V1..V6. Relapse x visit
#      likewise within AD (before relapse). dISF-significant proteins are classified by temporal profile
#      (resolving / persistent / late-rising / reversing / fluctuating) and the same proteins' serum profiles compared.
#   3  Compartment x group: per subject-visit, dISF minus serum (centred per pair, as in step 10); AD vs healthy on
#      this difference = proteins whose disease effect differs between dISF and serum.
#   4  Within-subject paired correlation of serum and dISF levels (from step 06): tracks serum = systemic
#      spill-over candidate; does not track (and relatively enriched in dISF, step 10) = local production candidate.
#   5  Signature sets (dISF-only, serum-only, shared-concordant, shared-discordant; FDR and nominal), Reactome / GO
#      over-representation against the assayed panel, tissue origin from the Human Protein Atlas (paths$hpa).
#   6  Relapse, predictive framing: dISF (and serum) at the visit before the relapse vs non-relapsers' cleared visits,
#      one value per patient.
# Out: output/signatures/  (signatures.xlsx + figures)

source("R/utils.R")
source("R/models.R")
source("R/design.R")
source("R/correlate.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clear_outputs(cfg, "signatures")
sg   <- cfg$signatures %||% list()
fdr  <- cfg$stats$fdr; mgn <- cfg$stats$min_group_n
out  <- \(...) file.path("signatures", ...)
assay_map <- clean |> distinct(OlinkID, Assay)
aname <- setNames(assay_map$Assay, assay_map$OlinkID)

# ---- proteins measured in both compartments (MicroAD serum) -------------------------------------------------------
ser_det <- clean |> filter(matrix == "Serum", cohort == "MicroAD", group %in% c("AD", "HC")) |>
  group_by(OlinkID, group) |> summarise(f = mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  group_by(OlinkID) |> summarise(serum_pct_above_LOD_best_group = round(100 * max(f, na.rm = TRUE), 1), .groups = "drop") |>
  # no LOD available (NaN / -Inf) = cannot be judged = kept, as in step 02
  mutate(serum_pct_above_LOD_best_group = if_else(is.finite(serum_pct_above_LOD_best_group), serum_pct_above_LOD_best_group, NA_real_),
         serum_detected = is.na(serum_pct_above_LOD_best_group) | serum_pct_above_LOD_best_group >= 100 * cfg$qc$min_detect_frac)
ser_ok <- intersect(rownames(wide$Serum), ser_det$OlinkID[ser_det$serum_detected])
isf_ok <- rownames(wide$ISF)
shared <- intersect(isf_ok, ser_ok)
isf_not_in_serum <- setdiff(isf_ok, ser_ok)
msg("dISF proteins: %d; detected in MicroAD serum as well: %d; dISF only (not detected in MicroAD serum): %d",
    length(isf_ok), length(shared), length(isf_not_in_serum))

# ---- samples ---------------------------------------------------------------------------------------------------------
isf_info <- visit_specs(isf_design(meta) |> filter(cohort == "MicroAD", SampleID %in% colnames(wide$ISF)))$info |>
  mutate(before_relapse = is.na(relapse_visit) | visit_num < relapse_visit)
ser_info <- serum_design(meta) |> filter(cohort == "MicroAD", SampleID %in% colnames(wide$Serum)) |>
  mutate(before_relapse = is.na(relapse_visit) | visit_num < relapse_visit)
visits <- sort(unique(isf_info$visit_num[isf_info$group == "AD"]))
hc_isf <- isf_info$SampleID[isf_info$site_cond %in% "HC"]
hc_ser <- ser_info$SampleID[ser_info$status %in% "HC"]

fit1 <- \(expr, info, ids, form, ct, name) {
  r <- fit_contrasts(expr, info |> filter(SampleID %in% ids), form, ct, name, mgn)
  if (is.null(r)) NULL else r |> select(OlinkID, logFC, t, P.Value, adj.P.Val, n_samples, n_subjects)
}
# relapse sets: V1 (lesional for everyone) or cleared visits, always before the relapse
rel_isf <- \(d) d |> filter(group == "AD", !is.na(relapse2), before_relapse, visit_num == 1 | state %in% "ex-lesional" | site_cond %in% "NL")
rel_ser <- \(d) d |> filter(group == "AD", !is.na(relapse2), before_relapse, visit_num == 1 | lesion_state %in% "cleared")

# ---- 1. effects per comparison, visit and compartment --------------------------------------------------------------
effects <- list()
add_eff <- function(comparison, visit, compartment, r) if (!is.null(r))
  effects[[length(effects) + 1]] <<- r |> mutate(comparison = comparison, visit = visit, compartment = compartment, .before = 1)
for (v in c(as.list(visits), list("all"))) {
  pooled <- identical(v, "all"); vlab <- if (pooled) "all visits" else paste0("V", v)
  inv <- \(d) if (pooled) rep(TRUE, nrow(d)) else d$visit_num %in% v
  re <- \(f) if (pooled) update(f, ~ . + (1 | SubjectID)) else f
  add_eff("AD vs healthy", vlab, "serum",
          fit1(wide$Serum[ser_ok, , drop = FALSE], ser_info, c(ser_info$SampleID[ser_info$status %in% "AD" & inv(ser_info)], hc_ser),
               re(~ 0 + status + plate), c(AD_vs_HC = "statusAD - statusHC"), paste("serum AD vs HC", vlab)))
  for (st in c("Lsite", "NL")) {
    lab <- c(Lsite = "lesion site", NL = "non-lesional skin")[[st]]
    add_eff(paste("AD vs healthy -", lab), vlab, "dISF",
            fit1(wide$ISF, isf_info, c(isf_info$SampleID[isf_info$site_cond %in% st & inv(isf_info)], hc_isf),
                 re(~ 0 + site_cond + plate), setNames(sprintf("site_cond%s - site_condHC", st), "AD_vs_HC"), paste("dISF", st, vlab)))
  }
  # relapse: pooled = cleared visits only (V1 is lesional for all), adjusted for weeks since V1
  ri <- rel_isf(isf_info); ri <- ri[ri$site_cond %in% "Lsite" & inv(ri) & (!pooled | ri$visit_num > 1), ] |> mutate(weeks = days_since_v1 / 7)
  rs <- rel_ser(ser_info); rs <- rs[inv(rs) & (!pooled | rs$visit_num > 1), ] |> mutate(weeks = days_since_v1 / 7)
  rf <- if (pooled) ~ 0 + relapse2 + weeks + plate + (1 | SubjectID) else ~ 0 + relapse2 + plate
  add_eff("relapse vs non-relapse (before relapse) - lesion site", vlab, "dISF",
          fit1(wide$ISF, ri, ri$SampleID, rf, c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"), paste("dISF relapse", vlab)))
  add_eff("relapse vs non-relapse (before relapse)", vlab, "serum",
          fit1(wide$Serum[ser_ok, , drop = FALSE], rs, rs$SampleID, rf, c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"), paste("serum relapse", vlab)))
}
effects <- bind_rows(effects) |> mutate(Assay = aname[OlinkID], .after = OlinkID)
save_csv(effects, cfg, out("effects_per_comparison.csv"))

# pair each dISF comparison with its serum counterpart
pairs_def <- tribble(~comparison, ~isf_comp, ~ser_comp,
  "AD vs healthy - lesion site", "AD vs healthy - lesion site", "AD vs healthy",
  "AD vs healthy - non-lesional skin", "AD vs healthy - non-lesional skin", "AD vs healthy",
  "relapse vs non-relapse (before relapse) - lesion site", "relapse vs non-relapse (before relapse) - lesion site", "relapse vs non-relapse (before relapse)")
vis_lv <- c(paste0("V", 1:12), "all visits")
joined <- pmap(pairs_def, \(comparison, isf_comp, ser_comp) {
  i <- effects |> filter(compartment == "dISF", comparison == isf_comp) |>
    transmute(visit, OlinkID, Assay, isf_logFC = logFC, isf_p = P.Value, isf_fdr = adj.P.Val)
  s <- effects |> filter(compartment == "serum", comparison == ser_comp) |>
    transmute(visit, OlinkID, ser_logFC = logFC, ser_p = P.Value, ser_fdr = adj.P.Val)
  full_join(i, s, by = c("visit", "OlinkID")) |> mutate(comparison = comparison, .before = 1)
}) |> bind_rows() |>
  mutate(Assay = coalesce(Assay, aname[OlinkID]), visit = factor(visit, levels = vis_lv) |> droplevels(),
         measured = case_when(OlinkID %in% shared ~ "both", OlinkID %in% isf_not_in_serum ~ "dISF only (not detected in serum)",
                              TRUE ~ "serum only (not detected in dISF)"))

conc_stats <- function(d) {
  d <- d |> filter(!is.na(isf_logFC), !is.na(ser_logFC))
  if (nrow(d) < 5) return(tibble(n = nrow(d), spearman_rho = NA_real_, p = NA_real_, slope_serum_on_dISF = NA_real_, pct_same_sign = NA_real_))
  ct <- suppressWarnings(cor.test(d$isf_logFC, d$ser_logFC, method = "spearman", exact = FALSE))
  tibble(n = nrow(d), spearman_rho = unname(ct$estimate), p = ct$p.value,
         slope_serum_on_dISF = unname(coef(lm(ser_logFC ~ isf_logFC, d))[2]),
         pct_same_sign = round(100 * mean(sign(d$isf_logFC) == sign(d$ser_logFC)), 1))
}
concordance <- joined |> filter(measured == "both") |> group_by(comparison, visit) |>
  group_modify(\(d, k) bind_rows(
    conc_stats(d) |> mutate(proteins = "all measured in both"),
    conc_stats(d |> filter(isf_fdr < fdr | ser_fdr < fdr)) |> mutate(proteins = sprintf("FDR < %g in either", fdr)),
    conc_stats(d |> filter(isf_p < 0.05 | ser_p < 0.05)) |> mutate(proteins = "p < 0.05 in either (sensitivity)"))) |>
  ungroup() |> relocate(comparison, visit, proteins)
save_csv(concordance, cfg, out("effect_concordance.csv"))
print(concordance |> filter(visit == "all visits") |> as.data.frame(), digits = 2)
not_in_serum <- joined |> filter(measured == "dISF only (not detected in serum)", !is.na(isf_logFC)) |>
  left_join(ser_det |> select(OlinkID, serum_pct_above_LOD_best_group), by = "OlinkID") |>
  select(comparison, visit, Assay, OlinkID, isf_logFC, isf_p, isf_fdr, serum_pct_above_LOD_best_group) |> arrange(comparison, visit, isf_p)
nis_sum <- not_in_serum |> group_by(comparison, visit) |>
  summarise(dISF_proteins_not_detected_in_serum = n(), of_which_dISF_FDR_sig = sum(isf_fdr < fdr, na.rm = TRUE),
            of_which_dISF_p05 = sum(isf_p < 0.05, na.rm = TRUE), .groups = "drop")

p <- joined |> filter(measured == "both", !is.na(isf_logFC), !is.na(ser_logFC)) |>
  mutate(sig = case_when(isf_fdr < fdr & ser_fdr < fdr ~ "FDR both", isf_fdr < fdr ~ "FDR dISF", ser_fdr < fdr ~ "FDR serum", TRUE ~ "n.s."),
         sig = factor(sig, levels = c("n.s.", "FDR dISF", "FDR serum", "FDR both"))) |> arrange(sig) |>
  ggplot(aes(isf_logFC, ser_logFC, colour = sig)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 0.6, alpha = 0.7) + geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
  geom_text(data = concordance |> filter(proteins == "all measured in both"),
            aes(x = -Inf, y = Inf, label = sprintf("rho %.2f\nslope %.2f\n%.0f%% same sign", spearman_rho, slope_serum_on_dISF, pct_same_sign)),
            hjust = -0.05, vjust = 1.1, size = 2.3, inherit.aes = FALSE) +
  scale_colour_manual(values = c(`n.s.` = "grey75", `FDR dISF` = "firebrick", `FDR serum` = "steelblue", `FDR both` = "purple3")) +
  facet_grid(str_wrap(comparison, 22) ~ visit, scales = "free") +
  labs(title = "Effect-size concordance: serum vs dISF log2 fold change (MicroAD; proteins measured in both)",
       x = "dISF log2FC", y = "serum log2FC", colour = NULL) + theme(legend.position = "bottom", strip.text = element_text(size = 7))
save_plot(p, cfg, out("effect_concordance.png"), width = 3 + 2.1 * n_distinct(joined$visit), height = 9)
report_plots <- list(concordance = p)
p <- concordance |> filter(visit != "all visits", !is.na(spearman_rho)) |>
  ggplot(aes(visit, spearman_rho, colour = proteins, group = proteins)) + geom_hline(yintercept = 0, colour = "grey60") +
  geom_line() + geom_point() + facet_wrap(~str_wrap(comparison, 40)) +
  labs(title = "Concordance of serum and dISF effects per visit (Spearman rho of log2FC)", x = NULL, y = "rho", colour = NULL) +
  theme(legend.position = "bottom")
save_plot(p, cfg, out("effect_concordance_by_visit.png"), width = 12, height = 4.5)
report_plots$concordance_by_visit <- p

# ---- 2. time-resolved models: group x visit ------------------------------------------------------------------------------
vk <- paste0("V", visits)
time_model <- function(expr, info, idcol_group, label) {
  d <- info |> mutate(gv = if_else(.data[[idcol_group]] == "HC", "HC", paste0("AD_", visit)))
  vv <- intersect(vk, unique(d$visit[d$gv != "HC"]))
  ct <- c(setNames(sprintf("gvAD_%s - gvHC", vv), paste0(vv, "_AD_vs_HC")),
          setNames(sprintf("gvAD_%s - gvAD_%s", vv[-1], vv[1]), paste0(vv[-1], "_minus_", vv[1], "_AD")))
  r <- fit_limma_f(expr, d, ~ 0 + gv + plate, "SubjectID", ct,
                   list(`group x visit (AD - healthy difference changes over visits)` = names(ct)[str_detect(names(ct), "_minus_")],
                        `AD vs healthy at any visit` = names(ct)[str_detect(names(ct), "_AD_vs_HC$")]), mgn)
  if (is.null(r)) return(NULL)
  list(effects = r$effects |> mutate(model = label), f = r$f |> mutate(model = label, n_samples = r$n_samples, n_subjects = r$n_subjects))
}
relapse_time_model <- function(expr, info, label) {
  d <- info |> mutate(rv = paste0(relapse2, "_", visit))
  vv <- vk[vapply(vk, \(v) all(c(sum(d$rv == paste0("relapse_", v)), sum(d$rv == paste0("non_relapse_", v))) >= 2), logical(1))]
  if (length(vv) < 2) return(NULL)
  d <- d |> filter(visit %in% vv)
  ct <- c(setNames(sprintf("rvrelapse_%s - rvnon_relapse_%s", vv, vv), paste0(vv, "_relapse_vs_non")),
          setNames(sprintf("(rvrelapse_%s - rvnon_relapse_%s) - (rvrelapse_%s - rvnon_relapse_%s)", vv[-1], vv[-1], vv[1], vv[1]),
                   paste0(vv[-1], "_vs_", vv[1], "_interaction")))
  r <- fit_limma_f(expr, d, ~ 0 + rv + plate, "SubjectID", ct,
                   list(`relapse x visit` = names(ct)[str_detect(names(ct), "_interaction$")]), 2)
  if (is.null(r)) return(NULL)
  list(effects = r$effects |> mutate(model = label), f = r$f |> mutate(model = label, n_samples = r$n_samples, n_subjects = r$n_subjects))
}
tm <- compact(list(
  time_model(wide$ISF, isf_info |> filter(site_cond %in% c("Lsite", "HC")), "site_cond", "dISF lesion site"),
  time_model(wide$ISF, isf_info |> filter(site_cond %in% c("NL", "HC")), "site_cond", "dISF non-lesional skin"),
  time_model(wide$Serum[ser_ok, , drop = FALSE], ser_info |> filter(status %in% c("AD", "HC")) |> mutate(grp = status), "grp", "serum"),
  relapse_time_model(wide$ISF, rel_isf(isf_info) |> filter(site_cond %in% "Lsite"), "relapse: dISF lesion site"),
  relapse_time_model(wide$Serum[ser_ok, , drop = FALSE], rel_ser(ser_info), "relapse: serum")))
tm_eff <- map(tm, "effects") |> bind_rows() |> mutate(Assay = aname[OlinkID], .after = OlinkID)
tm_f <- map(tm, "f") |> bind_rows() |> mutate(Assay = aname[OlinkID], .after = OlinkID)
tm_sum <- tm_f |> group_by(model, test, n_samples, n_subjects) |>
  summarise(proteins = n(), FDR_sig = sum(adj.P.Val < fdr), p05 = sum(P.Value < 0.05), .groups = "drop")
print(as.data.frame(tm_sum))

# temporal profiles of dISF-significant proteins (AD - healthy difference per visit)
classify <- function(e) {                          # e: effects V1..Vk (AD - healthy)
  e <- e[!is.na(e)]; if (length(e) < 3) return("not classifiable")
  s <- sign(e[which.max(abs(e))]); a <- e * s       # sign-aligned to the strongest effect
  early <- a[1]; late <- mean(tail(a, 2))
  if (early > 0.25 * max(a) && late < -0.25 * max(a)) return("reversing")
  if (early >= max(a) * 0.75 && late < 0.5 * early) return("high at V1, resolving")
  if (late >= 1.5 * max(early, 0) && late >= 0.75 * max(a)) return("late-rising")
  if (min(a) >= 0.4 * max(a)) return("persistent")
  "fluctuating"
}
profiles <- map(c("dISF lesion site", "dISF non-lesional skin"), \(md) {
  sig <- tm_f |> filter(model == md, adj.P.Val < fdr) |> distinct(OlinkID)
  if (!nrow(sig)) return(NULL)
  e <- tm_eff |> filter(model == md, str_detect(contrast, "_AD_vs_HC$"), OlinkID %in% sig$OlinkID) |>
    mutate(visit = str_extract(contrast, "^V\\d+"))
  cls <- e |> arrange(OlinkID, match(visit, vk)) |> group_by(OlinkID) |>
    summarise(profile = classify(logFC), direction = if_else(logFC[which.max(abs(logFC))] > 0, "up in AD", "down in AD"), .groups = "drop")
  ser_e <- tm_eff |> filter(model == "serum", str_detect(contrast, "_AD_vs_HC$")) |> mutate(visit = str_extract(contrast, "^V\\d+")) |>
    select(OlinkID, visit, serum_logFC = logFC, serum_fdr = adj.P.Val)
  e |> select(OlinkID, Assay, visit, dISF_logFC = logFC, dISF_fdr = adj.P.Val) |>
    left_join(cls, by = "OlinkID") |> left_join(ser_e, by = c("OlinkID", "visit")) |> mutate(site = md, .before = 1)
}) |> bind_rows()
if (nrow(profiles)) {
  prof_sum <- profiles |> mutate(sgn = if_else(direction == "up in AD", 1, -1)) |>
    group_by(site, profile) |>
    summarise(proteins = n_distinct(OlinkID), in_serum = n_distinct(OlinkID[!is.na(serum_logFC)]),
              attenuation_slope = if (sum(!is.na(serum_logFC)) >= 5) unname(coef(lm(I(sgn * serum_logFC) ~ 0 + I(sgn * dISF_logFC)))[1]) else NA_real_,
              .groups = "drop")
  pm <- profiles |> mutate(sgn = if_else(direction == "up in AD", 1, -1)) |>
    pivot_longer(c(dISF_logFC, serum_logFC), names_to = "compartment", values_to = "logFC") |>
    filter(!is.na(logFC)) |> mutate(compartment = sub("_logFC", "", compartment), aligned = sgn * logFC) |>
    group_by(site, profile, compartment, visit) |>
    summarise(mean = mean(aligned), se = sd(aligned) / sqrt(n()), n = n(), .groups = "drop") |>
    left_join(prof_sum |> select(site, profile, proteins), by = c("site", "profile")) |>
    mutate(visit = factor(visit, levels = vk), panel = sprintf("%s (%d)", profile, proteins))
  p <- ggplot(pm, aes(visit, mean, colour = compartment, group = compartment)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_line() + geom_pointrange(aes(ymin = mean - se, ymax = mean + se), size = 0.3) +
    scale_colour_manual(values = c(dISF = "firebrick", serum = "steelblue")) + facet_grid(site ~ panel) +
    labs(title = "Temporal profiles of dISF-significant proteins (AD - healthy, sign-aligned) and the same proteins in serum",
         subtitle = "Classes: high at V1 and resolving / persistent / late-rising / reversing / fluctuating; (n proteins). Serum line lower = attenuated.",
         x = NULL, y = "mean log2FC (aligned to the dISF direction)", colour = NULL) + theme(legend.position = "bottom")
  save_plot(p, cfg, out("temporal_profiles.png"), width = 14, height = 7)
  report_plots$profiles <- p
} else prof_sum <- tibble()

# ---- 3. compartment x group (paired subject-visit samples) ---------------------------------------------------------------
ser_pairs <- ser_info |> select(SubjectID, visit, serum_id = SampleID)
cxg <- map(c("Lsite", "NL"), \(st) {
  pp <- isf_info |> filter(site_cond %in% c(st, "HC")) |> select(SubjectID, visit, group, isf_id = SampleID, plate) |>
    inner_join(ser_pairs, by = c("SubjectID", "visit")) |>
    mutate(SampleID = paste(isf_id, serum_id, sep = "_"), grp = if_else(group == "HC", "HC", "AD"))
  if (nrow(pp) < 6 || length(shared) < 10) return(NULL)
  delta <- wide$ISF[shared, pp$isf_id, drop = FALSE] - wide$Serum[shared, pp$serum_id, drop = FALSE]
  delta <- sweep(delta, 2, apply(delta, 2, median, na.rm = TRUE)); colnames(delta) <- pp$SampleID
  lab <- c(Lsite = "lesion site", NL = "non-lesional skin")[[st]]
  bind_rows(
    fit_contrasts(delta, pp, ~ 0 + grp + (1 | SubjectID), c(compartment_x_group = "grpAD - grpHC"), paste("all visits,", lab), mgn),
    fit_contrasts(delta, pp |> filter(visit == "V1"), ~ 0 + grp, c(compartment_x_group = "grpAD - grpHC"), paste("V1,", lab), mgn))
}) |> bind_rows()
if (nrow(cxg)) {
  pooled_eff <- joined |> filter(visit == "all visits") |>
    mutate(site = if_else(str_detect(comparison, "lesion site"), "lesion site", "non-lesional skin")) |>
    filter(str_detect(comparison, "^AD vs healthy")) |> select(site, OlinkID, isf_logFC, ser_logFC)
  cxg <- cxg |> mutate(site = if_else(str_detect(model, "lesion site"), "lesion site", "non-lesional skin"), Assay = aname[OlinkID]) |>
    left_join(pooled_eff, by = c("site", "OlinkID")) |>
    transmute(model, Assay, OlinkID, interaction_log2 = logFC, t, P.Value, FDR = adj.P.Val, n_samples, n_subjects,
              dISF_AD_vs_HC_all_visits = isf_logFC, serum_AD_vs_HC_all_visits = ser_logFC,
              pattern = case_when(FDR >= fdr ~ "same disease effect in both (no interaction)",
                                  sign(dISF_AD_vs_HC_all_visits) != sign(serum_AD_vs_HC_all_visits) & abs(dISF_AD_vs_HC_all_visits) > 0.1 & abs(serum_AD_vs_HC_all_visits) > 0.1 ~ "opposite direction",
                                  abs(dISF_AD_vs_HC_all_visits) > abs(serum_AD_vs_HC_all_visits) ~ "stronger in dISF",
                                  TRUE ~ "stronger in serum")) |> arrange(model, P.Value)
  save_csv(cxg, cfg, out("compartment_x_group.csv"))
  print(count(cxg, model, pattern))
}

# ---- 4. within-subject paired correlation (step 06) ------------------------------------------------------------------------
corr_path <- file.path(cfg$paths$output, "isf_serum", "isf_serum_correlation.csv")
enr_path <- file.path(cfg$paths$output, "matrix_comparison", "relative_enrichment.csv")
tracking <- if (file.exists(corr_path)) {
  enr <- if (file.exists(enr_path)) read_csv(enr_path, show_col_types = FALSE) |>
    mutate(site = case_when(str_detect(model, "AD lesional") ~ "L", str_detect(model, "AD non-lesional") ~ "NL")) |>
    filter(!is.na(site)) |> select(site, OlinkID, rel_log2_isf_vs_serum, enrichment = direction) else NULL
  read_csv(corr_path, show_col_types = FALSE) |>
    (\(d) if (is.null(enr)) d else left_join(d, enr, by = c("site", "OlinkID")))() |>
    mutate(tracks_serum = coalesce(fdr_within < fdr & r_within > 0, FALSE) | coalesce(fdr_between < fdr & r_between > 0, FALSE),
           origin_hint = case_when(tracks_serum ~ "tracks serum - systemic spill-over candidate",
                                   enrichment %in% "enriched in dISF" ~ "does not track serum, enriched in dISF - local production candidate",
                                   TRUE ~ "does not track serum"),
           site = recode(site, L = "lesion site", NL = "non-lesional / healthy skin")) |>
    select(site, Assay, OlinkID, n_pairs, n_subjects, r_within, fdr_within, r_between, fdr_between, any_of(c("rel_log2_isf_vs_serum", "enrichment")),
           origin_hint) |> arrange(site, fdr_within)
} else { msg("Step 06 results not found - paired correlation skipped."); NULL }

# ---- 5. signature sets, enrichment, tissue origin ---------------------------------------------------------------------------
hpa_path <- cfg$paths$hpa %||% "data/reference/proteinatlas.tsv.zip"
hpa <- hpa_annotate(unique(str_extract(assay_map$Assay, "^[^_]+")), hpa_path)
if (!nrow(hpa)) msg("No Human Protein Atlas file at %s - tissue origin not annotated (see tools/download_hpa.R).", hpa_path)
sets_of <- function(d, tier) {
  sig_i <- if (tier == "FDR") d$isf_fdr < fdr else d$isf_p < 0.05
  sig_s <- if (tier == "FDR") d$ser_fdr < fdr else d$ser_p < 0.05
  d |> mutate(tier = tier, set = case_when(
    sig_i %in% TRUE & sig_s %in% TRUE & sign(isf_logFC) == sign(ser_logFC) ~ "shared-concordant",
    sig_i %in% TRUE & sig_s %in% TRUE ~ "shared-discordant",
    sig_i %in% TRUE & measured == "both" ~ "dISF-only",
    sig_i %in% TRUE ~ "dISF-only (not detected in serum)",
    sig_s %in% TRUE & measured == "both" ~ "serum-only",
    sig_s %in% TRUE ~ "serum-only (not detected in dISF)",
    TRUE ~ NA_character_)) |> filter(!is.na(set))
}
sig_sets <- joined |> filter(visit == "all visits") |>
  (\(d) bind_rows(sets_of(d, "FDR"), sets_of(d, "nominal")))() |>
  add_hpa(hpa) |>
  left_join(if (is.null(tracking)) tibble(OlinkID = character()) else
              tracking |> filter(site == "lesion site") |> select(OlinkID, origin_hint_lesion_site = origin_hint), by = "OlinkID") |>
  select(tier, comparison, set, Assay, OlinkID, isf_logFC, isf_p, isf_fdr, ser_logFC, ser_p, ser_fdr, any_of(c("hpa_origin", "hpa_secretome", "origin_hint_lesion_site"))) |>
  arrange(tier, comparison, set, pmin(isf_p, ser_p, na.rm = TRUE))
set_sum <- sig_sets |> count(tier, comparison, set) |> pivot_wider(names_from = set, values_from = n, values_fill = 0)
origin_sum <- if ("hpa_origin" %in% names(sig_sets)) sig_sets |> filter(!is.na(hpa_origin)) |> count(tier, comparison, set, hpa_origin) else tibble()
print(as.data.frame(set_sum))

colls <- unlist(sg$collections %||% c("C2:CP:REACTOME", "C5:GO:BP"))
gsets <- map(colls, \(cl) {
  parts <- str_split_fixed(cl, ":", 2)
  g <- if (parts[2] == "") msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1])
       else msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = parts[2])
  split(g$gene_symbol, g$gs_name)
}) |> unlist(recursive = FALSE)
genes_of <- \(a) unique(unlist(str_split(a, "_")))
ora <- sig_sets |> group_by(tier, comparison, set) |> group_modify(\(d, k) {
  universe <- genes_of(joined |> filter(comparison == k$comparison, visit == "all visits", !is.na(isf_p) | !is.na(ser_p)) |> pull(Assay))
  g <- genes_of(d$Assay)
  if (length(g) < 3) return(tibble())
  r <- fgsea::fora(gsets, genes = g, universe = universe, minSize = sg$min_size %||% 5, maxSize = cfg$enrichment$max_size %||% 500)
  if (!nrow(r)) return(tibble())
  as_tibble(r) |> mutate(overlapGenes = map_chr(overlapGenes, paste, collapse = ";"), universe_genes = length(universe)) |> arrange(pval)
}) |> ungroup()
msg("Over-representation: %d set x gene-set tests, %d with padj < %g", nrow(ora), sum(ora$padj < fdr, na.rm = TRUE), fdr)

# ---- 6. relapse, predictive framing ----------------------------------------------------------------------------------------
pre <- function(info, cleared) {               # one value set per patient: relapsers = visit before relapse
  info <- info |> filter(group == "AD", !is.na(relapse2))
  keep <- (info$relapse2 == "relapse" & !is.na(info$relapse_visit) & info$visit_num %in% (info$relapse_visit - 1) & cleared(info)) |
    (info$relapse2 == "non_relapse" & cleared(info))
  info[keep, ] |> select(SubjectID, relapse2, SampleID)
}
patient_matrix <- function(expr, sel) {
  ids <- split(sel$SampleID, sel$SubjectID)
  m <- sapply(ids, \(s) rowMeans(expr[, s, drop = FALSE], na.rm = TRUE))
  list(m = m, info = sel |> distinct(SubjectID, relapse2) |> mutate(SampleID = SubjectID))
}
pred_sets <- compact(list(
  `dISF lesion site` = { s <- pre(isf_info |> filter(site_cond %in% "Lsite"), \(d) d$state %in% "ex-lesional"); if (nrow(s)) patient_matrix(wide$ISF, s) },
  `dISF non-lesional skin` = { s <- pre(isf_info |> filter(site_cond %in% "NL"), \(d) d$lesion_state %in% "cleared"); if (nrow(s)) patient_matrix(wide$ISF, s) },
  `serum` = { s <- pre(ser_info, \(d) d$lesion_state %in% "cleared"); if (nrow(s)) patient_matrix(wide$Serum[ser_ok, , drop = FALSE], s) }))
if (all(c("dISF lesion site", "dISF non-lesional skin") %in% names(pred_sets))) {
  a <- pred_sets$`dISF lesion site`; b <- pred_sets$`dISF non-lesional skin`
  common <- intersect(colnames(a$m), colnames(b$m)); rows <- intersect(rownames(a$m), rownames(b$m))
  pred_sets$`dISF ex-lesional minus non-lesional` <- list(m = a$m[rows, common, drop = FALSE] - b$m[rows, common, drop = FALSE],
                                                          info = a$info |> filter(SubjectID %in% common))
}
predictive <- imap(pred_sets, \(ps, nm) {
  r <- fit_contrasts(ps$m, ps$info, ~ 0 + relapse2, c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"), nm, mgn)
  if (is.null(r)) return(NULL)
  r |> mutate(set = nm, n_relapse = sum(ps$info$relapse2 == "relapse"), n_non = sum(ps$info$relapse2 == "non_relapse"))
}) |> bind_rows()
if (nrow(predictive)) {
  predictive <- predictive |> mutate(Assay = aname[OlinkID]) |>
    transmute(set, Assay, OlinkID, logFC, t, P.Value, FDR = adj.P.Val, n_relapse, n_non) |> arrange(set, P.Value)
  pred_sum <- predictive |> group_by(set, n_relapse, n_non) |>
    summarise(proteins = n(), FDR_sig = sum(FDR < fdr), p05 = sum(P.Value < 0.05), expected_p05_by_chance = round(0.05 * n()), .groups = "drop")
  print(as.data.frame(pred_sum))
  p <- ggplot(predictive, aes(logFC, -log10(P.Value), colour = FDR < fdr)) + geom_point(size = 0.7, alpha = 0.7) +
    geom_text(data = \(d) d |> group_by(set) |> slice_min(P.Value, n = 8), aes(label = Assay), size = 2.4, vjust = -0.6, colour = "black", check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"), labels = c(`FALSE` = "n.s.", `TRUE` = sprintf("FDR < %g", fdr))) +
    facet_wrap(~set) + labs(title = "Relapse, predictive framing: the visit before the relapse vs non-relapsers' cleared visits (one value per patient)",
                            x = "log2FC (relapse - non-relapse)", y = "-log10 p", colour = NULL) + theme(legend.position = "bottom")
  save_plot(p, cfg, out("relapse_predictive_volcano.png"), width = 12, height = 7)
  report_plots$predictive <- p
} else pred_sum <- tibble()

# ---- workbook ------------------------------------------------------------------------------------------------------------------
readme <- tibble(sheet = c("README", "concordance", "not_detected_in_serum", "not_detected_in_serum_list", "time_models_F", "time_models_summary",
                           "temporal_profiles", "temporal_profiles_summary", "compartment_x_group", "paired_correlation",
                           "signature_sets", "signature_set_counts", "signature_origin_counts", "enrichment", "relapse_predictive",
                           "relapse_predictive_summary", "effects_all"),
                 content = c(sprintf("MicroAD only (RELAD, RELAD2, LEIP excluded). %d proteins in dISF, %d measured in both (dISF filter + detected in MicroAD serum), %d dISF-only. FDR < %g.",
                                     length(isf_ok), length(shared), length(isf_not_in_serum), fdr),
                             "1: Spearman rho, slope (serum on dISF) and % same sign of log2FC, per comparison and visit; all shared, FDR in either, p < 0.05 in either",
                             "1: per comparison, dISF proteins not detected in serum and how many of them are significant in dISF",
                             "1: the list of those proteins with their dISF effects",
                             "2: F-tests per protein: group x visit interaction and 'AD vs healthy at any visit' (dISF per site, serum); relapse x visit",
                             "2: number of significant proteins per model and test",
                             "2: dISF-significant proteins with their temporal profile class, dISF and serum effects per visit",
                             "2: proteins per profile class; attenuation_slope < 1 = serum shows a weaker version of the dISF profile",
                             "3: interaction = (AD - healthy in dISF) - (AD - healthy in serum) on paired, pair-centred differences; pattern from the pooled effects",
                             "4: within-subject (rmcorr) and between-subject correlation of dISF and serum levels (step 06; MicroAD incl. the 3 CPUO patients, healthy skin counted as non-lesional) with origin hint",
                             "5: protein lists: dISF-only, serum-only, shared-concordant, shared-discordant (FDR and nominal), with HPA tissue origin",
                             "5: number of proteins per set", "5: HPA origin per set (empty if no HPA file)",
                             "5: Reactome / GO over-representation per set, background = assayed panel of that comparison",
                             "6: relapsers at the visit before relapse vs non-relapsers (cleared visits, patient mean); all proteins",
                             "6: number of significant proteins (and expected at p < 0.05 by chance)",
                             "all fitted effects per comparison, visit and compartment"))
writexl::write_xlsx(Filter(\(x) !is.null(x) && nrow(x), list(
  README = readme, concordance = concordance, not_detected_in_serum = nis_sum, not_detected_in_serum_list = not_in_serum,
  time_models_F = tm_f, time_models_summary = tm_sum, temporal_profiles = profiles, temporal_profiles_summary = prof_sum,
  compartment_x_group = cxg, paired_correlation = tracking, signature_sets = sig_sets, signature_set_counts = set_sum,
  signature_origin_counts = origin_sum, enrichment = ora, relapse_predictive = predictive, relapse_predictive_summary = pred_sum,
  effects_all = effects)), out_path(cfg, out("signatures.xlsx")))
saveRDS(list(concordance = concordance, nis_sum = nis_sum, tm_sum = tm_sum, prof_sum = prof_sum, set_sum = set_sum,
             cxg_sum = if (nrow(cxg)) count(cxg, model, pattern) else tibble(), pred_sum = pred_sum,
             n = c(isf = length(isf_ok), shared = length(shared), isf_only = length(isf_not_in_serum)), hpa = nrow(hpa) > 0,
             plots = report_plots, sets_fdr = sig_sets |> filter(tier == "FDR"), predictive = predictive),
        out_path(cfg, out("summary.rds")))
msg("Serum vs dISF signatures: %s", file.path(cfg$paths$output, "signatures"))
