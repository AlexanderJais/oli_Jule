# 19 - Which measured proteins correlate with TNFRSF9 (CD137 / 4-1BB) in dISF?
# Protein and partners are set in config.yml (tnfrsf9). Skin sites: the tracked lesion site
# (lesional / ex-lesional) and the non-lesional site; healthy skin is the healthy group of both.
#   1. Targeted: TNFRSF9 vs each partner (IL33, IL4, ...), Spearman rho / p / n pooled (AD + healthy),
#      within each group (AD, healthy; relapse, non-relapse) and per visit. A correlation that is
#      seen only when pooling is driven by group differences, not by a protein-protein relation.
#   2. Proteome-wide: TNFRSF9 vs every measured dISF protein, adjusted for group (skin state),
#      visit and plate: mixed model (dream, subject random effect; ranking and FDR) and partial
#      Spearman correlation (ranks, same covariates; effect size).
#   3. Longitudinal: do changes between visits track within the same patient? Delta-delta Spearman
#      and repeated-measures correlation (= correlation of subject-centred values), plus a mixed
#      model on subject-centred TNFRSF9 adjusted for skin state.
#   4. Detectability: % of samples above LOD per site and group. Below-LOD values are used as
#      measured (Olink recommendation; primary). Sensitivity: below-LOD values set to the LOD (they
#      tie at the lowest rank) and only samples with both proteins above LOD.
#   5. Cross-compartment (MicroAD only, no RELAD / RELAD2): dISF vs serum TNFRSF9 in the same
#      subject and visit, and TNFRSF9 vs the partners in serum.
# Out: output/tnfrsf9/  TNFRSF9_correlations.xlsx, answers.csv, proteome_wide_*.csv, figures

source("R/utils.R")
source("R/models.R")
source("R/design.R")
source("R/correlation.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
clear_outputs(cfg, "tnfrsf9")
mgn <- cfg$stats$min_group_n; fdr <- cfg$stats$fdr
tn  <- cfg$tnfrsf9 %||% list()
tp_name   <- tn$protein %||% "TNFRSF9"
partners  <- unlist(tn$partners %||% c("IL33", "IL4", "CSF2", "IL6", "IL18", "CXCL8", "IL1RL1", "KIT", "KITLG",
                                       "TPSAB1", "TPSB2", "FCER1A"))
highlight <- unlist(tn$highlight %||% c("IL33", "IL4"))
min_det   <- cfg$qc$min_detect_frac
fmt_p <- \(p) ifelse(is.na(p), "n/a", ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))

assay_map <- clean |> distinct(OlinkID, Assay)
tp <- find_assay(assay_map, tp_name)
if (is.na(tp)) stop(tp_name, " is not in the data.")
p_oid <- setNames(map_chr(partners, \(n) find_assay(assay_map, n)), partners)
not_measured <- names(p_oid)[is.na(p_oid)]
p_oid <- p_oid[!is.na(p_oid) & p_oid != tp]
p_oid <- p_oid[!duplicated(p_oid)]                     # e.g. TPSAB1 and TPSB2 on one assay
lab <- \(o) assay_map$Assay[match(o, assay_map$OlinkID)]
msg("%s = %s; partners measured: %s; not measured: %s", tp_name, tp, paste(names(p_oid), collapse = ", "),
    if (length(not_measured)) paste(not_measured, collapse = ", ") else "-")
hl_oid <- p_oid[names(p_oid) %in% highlight]

# ---- samples: dISF per site set; MicroAD serum ------------------------------------------------------------------
isf <- isf_design(meta) |>
  filter(group %in% c("AD", "HC"), SampleID %in% clean$SampleID[clean$matrix == "ISF"]) |>
  mutate(grp = if_else(group == "AD", "AD", "healthy"),
         rel = case_when(group == "AD" & relapse %in% "relapse" ~ "relapse",
                         group == "AD" & relapse %in% "non-relapse" ~ "non-relapse"),
         colour_grp = case_when(grp == "healthy" ~ "healthy", !is.na(rel) ~ paste("AD", rel), TRUE ~ "AD other (dropout)"))
site_df <- bind_rows(isf |> filter(site == "L" | group == "HC") |> mutate(site_set = "lesional site"),
                     isf |> filter(site == "NL" | group == "HC") |> mutate(site_set = "non-lesional site"))
ser <- serum_design(meta) |>
  filter(cohort == "MicroAD", group %in% c("AD", "HC"), SampleID %in% clean$SampleID[clean$matrix == "Serum"]) |>
  mutate(grp = if_else(group == "AD", "AD", "healthy"),
         rel = case_when(group == "AD" & relapse %in% "relapse" ~ "relapse",
                         group == "AD" & relapse %in% "non-relapse" ~ "non-relapse"),
         colour_grp = case_when(grp == "healthy" ~ "healthy", !is.na(rel) ~ paste("AD", rel), TRUE ~ "AD other (dropout)"),
         site_set = "serum")

lv <- clean |> filter(OlinkID %in% c(tp, p_oid)) |> select(SampleID, matrix, OlinkID, value, LOD, below_lod)
getv <- \(o, mx, nm) lv |> filter(OlinkID == o, matrix == mx) |>
  transmute(SampleID, "{nm}" := value, "{nm}_lod" := LOD, "{nm}_below" := below_lod)
pair_data <- \(samples, o, mx = "ISF") samples |> inner_join(getv(tp, mx, "x"), by = "SampleID") |>
  inner_join(getv(o, mx, "y"), by = "SampleID")

# ---- 4. detectability ---------------------------------------------------------------------------------------------
det_summ <- \(d) d |> summarise(n = sum(!is.na(value)), n_above_LOD = sum(!below_lod, na.rm = TRUE),
                                pct_above_LOD = round(100 * mean(!below_lod, na.rm = TRUE), 1),
                                LOD_median = median(LOD, na.rm = TRUE), .groups = "drop")
detection <- map(c(tp, p_oid), \(o) {
  d <- bind_rows(site_df |> inner_join(lv |> filter(OlinkID == o, matrix == "ISF"), by = "SampleID"),
                 ser |> inner_join(lv |> filter(OlinkID == o, matrix == "Serum"), by = "SampleID"))
  bind_rows(d |> group_by(site_set) |> det_summ() |> mutate(group = "all"),
            d |> group_by(site_set, group = grp) |> det_summ(),
            d |> filter(!is.na(rel)) |> group_by(site_set, group = paste("AD", rel)) |> det_summ(),
            d |> filter(grp == "AD", !is.na(visit)) |> group_by(site_set, group = paste("AD", visit)) |> det_summ()) |>
    mutate(protein = lab(o), OlinkID = o, .before = 1)
}) |> bind_rows() |>
  mutate(compartment = if_else(site_set == "serum", "serum (MicroAD)", "dISF"), .after = OlinkID)
det_isf_overall <- clean |> filter(matrix == "ISF", OlinkID %in% c(tp, p_oid), SampleID %in% site_df$SampleID) |>
  group_by(OlinkID) |> summarise(pct_above_dISF = 100 * mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  mutate(protein = lab(OlinkID))
low_det <- det_isf_overall |> filter(pct_above_dISF < 100 * min_det)
lod_note <- if (nrow(low_det)) sprintf("Largely below LOD in dISF (< %d%% of samples above LOD): %s. Correlations involving them rest mostly on values in the noise range and should not be interpreted as a protein relationship.",
                                       round(100 * min_det), paste(sprintf("%s (%.0f%%)", low_det$protein, low_det$pct_above_dISF), collapse = ", ")) else
  "All analysed proteins are above LOD in at least half of the dISF samples."
msg(lod_note)

# ---- 1. targeted correlations ----------------------------------------------------------------------------------------
corr_block <- function(d) {
  cx <- censor_at_lod(d$x, d$x_lod); cy <- censor_at_lod(d$y, d$y_lod)
  above <- !(d$x_below %in% TRUE) & !(d$y_below %in% TRUE)
  bind_cols(spearman_row(d$x, d$y),
            spearman_row(cx, cy) |> select(rho_censored = rho, p_censored = p),
            spearman_row(d$x[above], d$y[above]) |> select(n_both_above_LOD = n, rho_above_LOD = rho, p_above_LOD = p),
            rmcorr_row(d$x, d$y, d$SubjectID)) |>
    mutate(n_subjects = n_distinct(d$SubjectID), .after = n)
}
strata_of <- function(d) {
  vis <- sort(unique(d$visit_num[d$grp == "AD"]))
  states <- if ("state" %in% names(d)) sort(unique(d$state[d$grp == "AD" & d$state %in% c("lesional", "ex-lesional")])) else character()
  c(list(`pooled (AD + healthy)` = d, AD = filter(d, grp == "AD"), healthy = filter(d, grp == "healthy"),
         `AD relapse` = filter(d, rel %in% "relapse"), `AD non-relapse` = filter(d, rel %in% "non-relapse")),
    # at the lesion site the skin state differs between samples of the AD group, so it is a group too
    if (length(states) > 1) setNames(map(states, \(s) filter(d, grp == "AD", state == s)), paste("AD", states)),
    list(`V1 pooled (AD + healthy)` = filter(d, visit_num %in% 1)),
    setNames(map(vis, \(v) filter(d, grp == "AD", visit_num == v)), paste0("AD V", vis)))
}
stratum_type <- \(s) case_when(str_detect(s, "pooled") & !str_detect(s, "^V") ~ "a) pooled",
                               s %in% c("AD", "healthy", "AD relapse", "AD non-relapse", "AD lesional", "AD ex-lesional") ~ "b) within group",
                               TRUE ~ "c) per visit")
corr_table <- function(samples, compartment, mx) {
  map(p_oid, \(o) {
    d <- pair_data(samples, o, mx)
    map(unique(d$site_set), \(st) {
      imap(strata_of(d |> filter(site_set == st)), \(dd, nm) corr_block(dd) |> mutate(stratum = nm, .before = 1)) |>
        bind_rows() |> mutate(site = st, .before = 1)
    }) |> bind_rows() |> mutate(partner = lab(o), OlinkID = o, .before = 1)
  }) |> bind_rows() |>
    mutate(compartment = compartment, level = stratum_type(stratum), .after = OlinkID)
}
targeted <- corr_table(site_df, "dISF", "ISF")

finest_hits <- function(stratum, p) {
  fine <- c("healthy", "AD lesional", "AD ex-lesional")
  if (!any(stratum %in% c("AD lesional", "AD ex-lesional"))) fine <- c(fine, "AD")
  paste(stratum[stratum %in% fine & coalesce(p < 0.05, FALSE)], collapse = ", ")
}
verdict_of <- function(tab) {
  tab |> group_by(partner, OlinkID, compartment, site) |>
    summarise(pooled_rho = rho[stratum == "pooled (AD + healthy)"], pooled_p = p[stratum == "pooled (AD + healthy)"],
              pooled_n = n[stratum == "pooled (AD + healthy)"],
              AD_rho = rho[stratum == "AD"], AD_p = p[stratum == "AD"], AD_n = n[stratum == "AD"],
              healthy_rho = rho[stratum == "healthy"], healthy_p = p[stratum == "healthy"], healthy_n = n[stratum == "healthy"],
              within_subject_r = r_within[stratum == "AD"], within_subject_p = p_within[stratum == "AD"],
              # finest groups: healthy, and AD within one skin state (lesion site: lesional / ex-lesional separately)
              finest_p05 = finest_hits(stratum, p),
              visits_p05 = sum(p[level == "c) per visit"] < 0.05, na.rm = TRUE), visits_tested = sum(!is.na(p[level == "c) per visit"])),
              pooled_rho_above_LOD = rho_above_LOD[stratum == "pooled (AD + healthy)"], .groups = "drop") |>
    mutate(within_group = finest_p05 != "",
           verdict = case_when(within_group ~ sprintf("correlated within groups: a relationship between the proteins (p < 0.05 in %s)", finest_p05),
                               coalesce(pooled_p < 0.05, FALSE) | coalesce(AD_p < 0.05, FALSE) ~
                                 "pooled only: driven by group / skin-state differences, not within groups",
                               TRUE ~ "no correlation shown"))
}
verdicts <- verdict_of(targeted)

# ---- 2. proteome-wide ----------------------------------------------------------------------------------------------
ew <- clean |> filter(matrix == "ISF", SampleID %in% site_df$SampleID) |> select(OlinkID, SampleID, value) |>
  pivot_wider(names_from = SampleID, values_from = value)
E <- as.matrix(ew[, -1]); rownames(E) <- ew$OlinkID
E <- E[rowSums(!is.na(E)) >= 8, , drop = FALSE]
keep_isf <- clean |> filter(matrix == "ISF") |> distinct(OlinkID, keep)
pct_isf <- clean |> filter(matrix == "ISF", SampleID %in% site_df$SampleID) |> group_by(OlinkID) |>
  summarise(pct_above_LOD_dISF = round(100 * mean(!below_lod, na.rm = TRUE), 1), .groups = "drop")
tnf_val <- clean |> filter(matrix == "ISF", OlinkID == tp) |> select(SampleID, tnf = value)

partial_spearman <- function(E, info, X) {
  map(rownames(E), \(o) {
    y <- E[o, info$SampleID]; ok <- !is.na(y) & !is.na(info$tnf)
    if (sum(ok) < 8) return(NULL)
    Xo <- estimable_design(X[ok, , drop = FALSE]); q <- qr(Xo)
    rx <- qr.resid(q, rank(info$tnf[ok])); ry <- qr.resid(q, rank(y[ok]))
    r <- suppressWarnings(cor(rx, ry)); dfree <- sum(ok) - ncol(Xo) - 1
    tt <- r * sqrt(dfree / max(1 - r^2, 1e-12))
    tibble(OlinkID = o, n = sum(ok), partial_rho = r, partial_p = 2 * pt(-abs(tt), dfree),
           rho_unadjusted = suppressWarnings(cor(rank(info$tnf[ok]), rank(y[ok]))))
  }) |> bind_rows()
}
proteome <- map(c("lesional site", "non-lesional site"), \(st) {
  info <- site_df |> filter(site_set == st) |> inner_join(tnf_val, by = "SampleID") |>
    filter(!is.na(tnf), !is.na(cond)) |> mutate(visit = factor(visit))
  msg("Proteome-wide, %s: %d samples, %d proteins", st, nrow(info), nrow(E) - 1)
  Eo <- E[setdiff(rownames(E), tp), info$SampleID, drop = FALSE]
  mm <- fit_contrasts(Eo, info, ~ tnf + cond + visit + plate + (1 | SubjectID), c(TNFRSF9 = "tnf"),
                      paste(tp_name, st), mgn)
  X <- model.matrix(~ cond + visit + plate, droplevels(as.data.frame(info)))
  ps <- partial_spearman(Eo, info, X)
  ps |> mutate(partial_fdr = p.adjust(partial_p, "BH")) |>
    left_join(if (is.null(mm)) tibble(OlinkID = character()) else
      mm |> transmute(OlinkID, mixed_slope = logFC, mixed_t = t, mixed_p = P.Value, mixed_fdr = adj.P.Val, mixed_method = method),
      by = "OlinkID") |>
    mutate(site = st, .before = 1)
}) |> bind_rows() |>
  left_join(assay_map, by = "OlinkID") |> left_join(keep_isf |> rename(detected_dISF = keep), by = "OlinkID") |>
  left_join(pct_isf, by = "OlinkID") |>
  mutate(partner = OlinkID %in% p_oid, highlight = OlinkID %in% hl_oid,
         rank_p = coalesce(mixed_p, partial_p)) |>
  group_by(site) |> arrange(rank_p, .by_group = TRUE) |> mutate(rank = row_number(), of = n()) |> ungroup() |>
  select(site, rank, of, Assay, OlinkID, n, partial_rho, partial_p, partial_fdr, mixed_slope, mixed_t, mixed_p, mixed_fdr,
         rho_unadjusted, detected_dISF, pct_above_LOD_dISF, partner, highlight, any_of("mixed_method"))
for (st in unique(proteome$site))
  save_csv(proteome |> filter(site == st), cfg, "tnfrsf9", sprintf("proteome_wide_%s.csv", str_replace_all(st, " ", "_")))
hl_rank <- proteome |> filter(partner) |> select(site, Assay, rank, of, partial_rho, mixed_p, mixed_fdr)
print(as.data.frame(hl_rank), digits = 3)

# ---- 3. longitudinal (within-subject) --------------------------------------------------------------------------------
longit <- map(p_oid, \(o) {
  d <- pair_data(site_df |> filter(grp == "AD", !is.na(visit_num)), o)
  map(unique(d$site_set), \(st) {
    dd <- d |> filter(site_set == st) |> arrange(SubjectID, visit_num)
    dl <- dd |> group_by(SubjectID) |>
      mutate(dx = x - lag(x), dy = y - lag(y), from = lag(visit_num)) |> ungroup() |> filter(!is.na(dx), !is.na(dy))
    cen <- dd |> group_by(SubjectID) |> filter(n() >= 2) |> mutate(x_c = x - mean(x)) |> ungroup() |> mutate(value = y)
    mm <- test_single(cen, ~ x_c + cond + (1 | SubjectID), c(within_slope = "x_c"), "subject-centred", mgn)
    bind_cols(spearman_row(dl$dx, dl$dy) |> rename(n_deltas = n, delta_rho = rho, delta_p = p),
              tibble(n_subjects_deltas = n_distinct(dl$SubjectID)),
              rmcorr_row(dd$x, dd$y, dd$SubjectID),
              tibble(mixed_within_slope = if (is.null(mm) || !nrow(mm)) NA_real_ else mm$estimate[1],
                     mixed_within_p = if (is.null(mm) || !nrow(mm)) NA_real_ else mm$p[1])) |>
      mutate(site = st, .before = 1)
  }) |> bind_rows() |> mutate(partner = lab(o), OlinkID = o, .before = 1)
}) |> bind_rows()

# ---- 5. cross-compartment (MicroAD only) ------------------------------------------------------------------------------
x_isf <- getv(tp, "ISF", "x"); x_ser <- getv(tp, "Serum", "y")
cc_pairs <- site_df |> filter(cohort == "MicroAD") |> inner_join(x_isf, by = "SampleID") |>
  inner_join(ser |> select(SubjectID, visit, serum_id = SampleID), by = c("SubjectID", "visit")) |>
  inner_join(x_ser |> rename(serum_id = SampleID), by = "serum_id")
cross <- map(unique(cc_pairs$site_set), \(st) {
  imap(strata_of(cc_pairs |> filter(site_set == st)), \(dd, nm) corr_block(dd) |> mutate(stratum = nm, .before = 1)) |>
    bind_rows() |> mutate(site = st, .before = 1)
}) |> bind_rows() |> mutate(comparison = sprintf("dISF %s vs serum %s (same subject and visit)", tp_name, tp_name), .before = 1) |>
  mutate(level = stratum_type(stratum), .after = stratum)
serum_corr <- corr_table(ser, "serum (MicroAD)", "Serum")
serum_verdicts <- verdict_of(serum_corr)

# ---- figures ----------------------------------------------------------------------------------------------------------
grp_cols <- c(healthy = "grey45", `AD non-relapse` = "steelblue", `AD relapse` = "firebrick", `AD other (dropout)` = "orange3")
scatter <- function(tab_samples, oids, mx, title) {
  d <- map(oids, \(o) pair_data(tab_samples, o, mx) |> mutate(partner = lab(o))) |> bind_rows() |>
    mutate(below = if_else(x_below %in% TRUE | y_below %in% TRUE, "below LOD (either)", "above LOD"))
  if (!nrow(d)) return(NULL)
  stats <- d |> group_by(partner, site_set) |>
    summarise(label = sprintf("pooled rho %.2f (p %s, n %d)\nAD %.2f (p %s) | healthy %.2f (p %s)",
                              spearman_row(x, y)$rho, fmt_p(spearman_row(x, y)$p), sum(!is.na(x) & !is.na(y)),
                              spearman_row(x[grp == "AD"], y[grp == "AD"])$rho, fmt_p(spearman_row(x[grp == "AD"], y[grp == "AD"])$p),
                              spearman_row(x[grp == "healthy"], y[grp == "healthy"])$rho,
                              fmt_p(spearman_row(x[grp == "healthy"], y[grp == "healthy"])$p)), .groups = "drop")
  ggplot(d, aes(y, x, colour = colour_grp)) +
    geom_point(aes(shape = below), size = 1.6, alpha = 0.8) +
    geom_smooth(aes(group = grp, linetype = grp), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5, colour = "grey20") +
    geom_text(data = stats, aes(x = -Inf, y = Inf, label = label), hjust = -0.03, vjust = 1.1, size = 2.5, inherit.aes = FALSE) +
    scale_colour_manual(values = grp_cols) + scale_shape_manual(values = c(`above LOD` = 16, `below LOD (either)` = 1)) +
    facet_grid(partner ~ site_set, scales = "free") +
    labs(title = title, x = "partner NPX", y = paste(tp_name, "NPX"), colour = "group", shape = NULL, linetype = "fit") +
    theme(legend.position = "bottom", legend.box = "vertical")
}
if (length(hl_oid)) {
  save_plot(scatter(site_df, hl_oid, "ISF", sprintf("%s vs %s in dISF, per site", tp_name, paste(names(hl_oid), collapse = " / "))),
            cfg, "tnfrsf9", "scatter_TNFRSF9_vs_IL33_IL4_dISF.png", width = 10, height = 3 + 3.2 * length(hl_oid))
  save_plot(scatter(ser, hl_oid, "Serum", sprintf("%s vs %s in serum (MicroAD)", tp_name, paste(names(hl_oid), collapse = " / "))),
            cfg, "tnfrsf9", "scatter_TNFRSF9_vs_IL33_IL4_serum.png", width = 6.5, height = 3 + 3.2 * length(hl_oid))
}
save_plot(scatter(site_df, p_oid, "ISF", sprintf("%s vs all partner proteins in dISF", tp_name)),
          cfg, "tnfrsf9", "scatter_TNFRSF9_vs_partners_dISF.png", width = 10, height = 3 + 2.8 * length(p_oid))

p <- proteome |> ggplot(aes(partial_rho, -log10(coalesce(mixed_p, partial_p)))) +
  geom_point(aes(colour = coalesce(mixed_fdr, partial_fdr) < fdr), size = 0.7, alpha = 0.7) +
  geom_point(data = \(d) filter(d, partner), colour = "black", shape = 1, size = 2.2) +
  geom_text(data = \(d) d |> filter(partner | rank <= 10), aes(label = Assay), size = 2.5, vjust = -0.7, check_overlap = TRUE) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), labels = c("FDR >= 0.05", "FDR < 0.05")) +
  facet_wrap(~site) +
  labs(title = sprintf("Proteins correlated with %s in dISF (adjusted for skin state, visit, plate; subject random effect)", tp_name),
       subtitle = "circles = partner proteins (config tnfrsf9$partners); labels = partners and top 10",
       x = "partial Spearman rho", y = "-log10 p (mixed model)", colour = NULL) + theme(legend.position = "bottom")
save_plot(p, cfg, "tnfrsf9", "proteome_wide_volcano.png", width = 12, height = 6)

dd_plot <- map(hl_oid, \(o) pair_data(site_df |> filter(grp == "AD", !is.na(visit_num)), o) |> arrange(SubjectID, visit_num) |>
                 group_by(SubjectID, site_set) |> mutate(dx = x - lag(x), dy = y - lag(y)) |> ungroup() |>
                 filter(!is.na(dx), !is.na(dy)) |> mutate(partner = lab(o))) |> bind_rows()
if (nrow(dd_plot)) {
  p <- ggplot(dd_plot, aes(dy, dx, colour = colour_grp)) + geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") + geom_point() +
    geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linewidth = 0.5) +
    scale_colour_manual(values = grp_cols) + facet_grid(partner ~ site_set, scales = "free") +
    labs(title = sprintf("Change between consecutive visits: %s vs partner (AD, same patient)", tp_name),
         x = "change in partner NPX", y = sprintf("change in %s NPX", tp_name), colour = NULL) + theme(legend.position = "bottom")
  save_plot(p, cfg, "tnfrsf9", "delta_delta_TNFRSF9.png", width = 9, height = 2.5 + 3 * length(hl_oid))
}
if (nrow(cc_pairs)) {
  p <- ggplot(cc_pairs, aes(y, x, colour = colour_grp)) + geom_point() +
    geom_line(aes(group = SubjectID), alpha = 0.3) + scale_colour_manual(values = grp_cols) + facet_wrap(~site_set) +
    labs(title = sprintf("%s: dISF vs serum, same subject and visit (MicroAD); lines join visits of one patient", tp_name),
         x = "serum NPX", y = "dISF NPX", colour = NULL) + theme(legend.position = "bottom")
  save_plot(p, cfg, "tnfrsf9", "TNFRSF9_dISF_vs_serum.png", width = 10, height = 5)
}

# ---- answers ------------------------------------------------------------------------------------------------------------
ans <- list()
add <- \(q, item, verdict, evidence) ans[[length(ans) + 1]] <<- tibble(question = q, item = item, verdict = verdict, evidence = evidence)
Q <- c(det = "4 Detectability (read first)", tgt = "1 Targeted correlations in dISF", pw = "2 Proteome-wide (adjusted for group, visit, plate)",
       lon = "3 Longitudinal (within patients)", cc = "5 Cross-compartment (MicroAD serum)")
for (o in c(tp, p_oid)) {
  dd <- detection |> filter(OlinkID == o, group %in% c("all", "AD", "healthy"))
  add(Q[["det"]], lab(o),
      { pa <- det_isf_overall$pct_above_dISF[det_isf_overall$OlinkID == o]
        if (!length(pa)) "not measured in dISF" else if (pa < 100 * min_det) sprintf("largely below LOD in dISF (%.0f%% above)", pa) else sprintf("detectable in dISF (%.0f%% above LOD)", pa) },
      paste(sprintf("%s %s: %.0f%% (n %d)", dd$site_set, dd$group, dd$pct_above_LOD, dd$n), collapse = "; "))
}
if (length(not_measured)) add(Q[["det"]], paste(not_measured, collapse = ", "), "not measured on the panel", "no assay with this name, OlinkID or UniProt")
add(Q[["det"]], "below-LOD handling", "values used as measured (primary)",
    paste("Olink recommends keeping below-LOD NPX values; Spearman uses ranks only. Sensitivity columns: rho_censored (below-LOD values set to the LOD, so they tie at the lowest rank) and rho_above_LOD (only samples with both proteins above LOD).", lod_note))
for (i in seq_len(nrow(verdicts))) {
  v <- verdicts[i, ]
  add(Q[["tgt"]], sprintf("%s vs %s (%s)", tp_name, v$partner, v$site),
      paste0(if (v$OlinkID %in% low_det$OlinkID) "CAUTION largely below LOD - " else "", v$verdict),
      sprintf("pooled rho %.2f (p %s, n %d); AD rho %.2f (p %s, n %d); healthy rho %.2f (p %s, n %d); within patients r %.2f (p %s); p < 0.05 at %d of %d visits",
              v$pooled_rho, fmt_p(v$pooled_p), v$pooled_n, v$AD_rho, fmt_p(v$AD_p), v$AD_n, v$healthy_rho, fmt_p(v$healthy_p),
              v$healthy_n, v$within_subject_r, fmt_p(v$within_subject_p), v$visits_p05, v$visits_tested))
}
for (st in unique(proteome$site)) {
  pw <- proteome |> filter(site == st)
  top <- pw |> filter(coalesce(mixed_fdr, partial_fdr) < fdr)
  add(Q[["pw"]], sprintf("%s: proteins correlated with %s", st, tp_name),
      sprintf("%d of %d proteins at FDR < %g", nrow(top), nrow(pw), fdr),
      if (nrow(top)) paste(sprintf("%s (rho %.2f)", head(top$Assay, 25), head(top$partial_rho, 25)), collapse = ", ") else "none")
  for (h in names(hl_oid)) {
    r <- pw |> filter(OlinkID == hl_oid[[h]])
    if (nrow(r)) add(Q[["pw"]], sprintf("%s: where %s falls", st, h), sprintf("rank %d of %d", r$rank, r$of),
                     sprintf("partial rho %.2f, p %s, FDR %s", r$partial_rho, fmt_p(coalesce(r$mixed_p, r$partial_p)),
                             fmt_p(coalesce(r$mixed_fdr, r$partial_fdr))))
  }
}
for (i in seq_len(nrow(longit))) {
  l <- longit[i, ]
  add(Q[["lon"]], sprintf("%s vs %s (%s)", tp_name, l$partner, l$site),
      if (coalesce(l$delta_p < 0.05, FALSE) || coalesce(l$p_within < 0.05, FALSE))
        sprintf("changes track within patients (%s)", if (coalesce(l$r_within, l$delta_rho, 0) > 0) "same direction" else "opposite direction") else "no within-patient tracking shown",
      sprintf("delta-delta rho %.2f (p %s, %d changes in %d patients); repeated-measures r %.2f (p %s); mixed model on subject-centred values, adjusted for skin state: slope %.2f (p %s)",
              l$delta_rho, fmt_p(l$delta_p), l$n_deltas, l$n_subjects_deltas, l$r_within, fmt_p(l$p_within),
              l$mixed_within_slope, fmt_p(l$mixed_within_p)))
}
for (st in unique(cross$site)) {
  c1 <- cross |> filter(site == st, stratum == "pooled (AD + healthy)"); c2 <- cross |> filter(site == st, stratum == "AD")
  add(Q[["cc"]], sprintf("dISF (%s) vs serum %s", st, tp_name),
      if (coalesce(c1$p < 0.05, FALSE) || coalesce(c2$p_within < 0.05, FALSE)) "dISF and serum TNFRSF9 correlate" else "no correlation between dISF and serum shown",
      sprintf("pooled rho %.2f (p %s, n %d pairs); AD rho %.2f (p %s); within patients r %.2f (p %s)",
              c1$rho, fmt_p(c1$p), c1$n, c2$rho, fmt_p(c2$p), c2$r_within, fmt_p(c2$p_within)))
}
for (i in seq_len(nrow(serum_verdicts))) {
  v <- serum_verdicts[i, ]
  add(Q[["cc"]], sprintf("serum: %s vs %s", tp_name, v$partner), v$verdict,
      sprintf("pooled rho %.2f (p %s, n %d); AD rho %.2f (p %s); healthy rho %.2f (p %s); within patients r %.2f (p %s)",
              v$pooled_rho, fmt_p(v$pooled_p), v$pooled_n, v$AD_rho, fmt_p(v$AD_p), v$healthy_rho, fmt_p(v$healthy_p),
              v$within_subject_r, fmt_p(v$within_subject_p)))
}
answers <- bind_rows(ans)
save_csv(answers, cfg, "tnfrsf9", "answers.csv")
save_csv(detection, cfg, "tnfrsf9", "LOD_summary.csv")
save_csv(targeted, cfg, "tnfrsf9", "targeted_correlations.csv")

methods <- tibble(topic = c("samples", "groups", "targeted", "proteome-wide", "longitudinal", "below LOD", "cross-compartment", "caveat"),
                  text = c(
  "dISF of MicroAD. Lesional site = tracked lesion site (lesional or ex-lesional) of AD patients; non-lesional site = non-lesional skin of AD patients. Healthy skin of healthy volunteers is the healthy group at both sites. CPUO is not included.",
  "AD vs healthy; within AD relapse vs non-relapse (manifest Relapse column; dropouts only in the AD and pooled strata).",
  "Spearman rho, p, n: a) pooled AD + healthy, b) within each group, c) per visit (AD; V1 also pooled with healthy). Samples of one patient at several visits are not independent: the pooled and within-group p-values over several visits are therefore optimistic; r_within (repeated-measures correlation) is the valid within-patient test.",
  "Every measured dISF protein (also those below the detection filter; see detected_dISF and pct_above_LOD_dISF) vs TNFRSF9, adjusted for skin state (group), visit and plate. Ranking and FDR: mixed model (dream: protein ~ TNFRSF9 + state + visit + plate + (1 | subject)). Effect size: partial Spearman rho (ranks residualised on the same covariates, without the random effect).",
  "AD patients: change between consecutive available visits in TNFRSF9 vs in the partner (delta-delta Spearman); repeated-measures correlation (correlation of subject-centred values); mixed model partner ~ subject-centred TNFRSF9 + skin state + (1 | subject).",
  "Below-LOD values are used as measured (primary). Sensitivity: set to the LOD (ties) and restricted to samples with both proteins above LOD.",
  "MicroAD only (RELAD / RELAD2 serum is not used): dISF TNFRSF9 vs serum TNFRSF9 at the same subject and visit, and TNFRSF9 vs the partners in serum.",
  "Relapse strata: 4 relapsing patients in MicroAD - exploratory."))
writexl::write_xlsx(Filter(\(x) !is.null(x) && nrow(x), list(
  answers = answers, methods = methods, LOD_summary = detection, targeted_dISF = targeted, targeted_summary = verdicts,
  proteome_wide = proteome, partners_in_ranking = hl_rank, longitudinal = longit, cross_compartment_TNFRSF9 = cross,
  serum_partners = serum_corr, serum_summary = serum_verdicts)),
  out_path(cfg, "tnfrsf9", "TNFRSF9_correlations.xlsx"))
for (q in unique(answers$question)) {
  message("\n", q)
  a <- answers |> filter(question == q)
  for (i in seq_len(nrow(a))) message(sprintf("  - %s: %s", a$item[i], a$verdict[i]))
}
msg("TNFRSF9 correlations: %s", file.path(cfg$paths$output, "tnfrsf9"))
