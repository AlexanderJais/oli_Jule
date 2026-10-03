# 13 - dISF comparisons visit by visit, and proteins regulated at ALL visits ("time course")
# At each visit (V1..V6, if >= stats$visit_min_subjects AD patients):
#   Lsite_vs_HC  tracked lesion site (lesional at V1/relapse, ex-lesional after clearing) vs healthy skin
#   Lsite_vs_NL  lesion site vs non-lesional skin of the same patients (paired)
#   NL_vs_HC     non-lesional AD skin vs healthy skin
# Healthy controls (one visit) are the reference at every visit. Same model as step 04
# (dream: subject random effect, plate fixed), with a volcano plot per visit.
# "Regulated at all visits": significant at every analysed visit, same direction
#   - strict: FDR < stats$fdr at every visit;  - nominal: p < 0.05 at every visit.
# Out: output/models/ISF_by_visit_*, output/visit_course/*

source("R/utils.R")
source("R/models.R")
source("R/design.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
clear_outputs(cfg, "models", "^ISF_by_visit"); clear_outputs(cfg, "visit_course")
expr <- wide$ISF
fdr  <- cfg$stats$fdr

vs <- visit_specs(isf_design(meta) |> filter(SampleID %in% colnames(expr)), cfg$stats$visit_min_subjects %||% 5)
info <- vs$info
comp <- info |> filter(group == "AD", site == "L") |> count(visit, state) |>
  pivot_wider(names_from = state, values_from = n, values_fill = 0)
save_csv(comp, cfg, "visit_course", "lesion_site_state_per_visit.csv")
msg("Lesion-site state per visit:"); print(as.data.frame(comp))

res <- run_model_specs(vs$specs, expr, info, cfg, "ISF_by_visit", assay_map)
if (!nrow(res)) stop("No per-visit model could be fitted.")
res <- res |> mutate(visit = factor(model, levels = paste0("V", 1:12)) |> droplevels(),
                     se = abs(logFC / t), ci_low = logFC - 1.96 * se, ci_high = logFC + 1.96 * se)

# ---- regulated at all visits -------------------------------------------------------------------------
consistency <- res |>
  group_by(contrast, OlinkID, Assay) |>
  summarise(n_visits = n(), n_fdr = sum(adj.P.Val < fdr), n_nominal = sum(P.Value < 0.05),
            n_up = sum(logFC > 0), n_down = sum(logFC < 0), mean_logFC = mean(logFC),
            max_p = max(P.Value), .groups = "drop") |>
  # a protein not modelled at some visit (e.g. too many missing values there) cannot be "regulated at all visits"
  group_by(contrast) |> mutate(visits_analysed = max(n_visits)) |> ungroup() |>
  mutate(same_direction = n_up == n_visits | n_down == n_visits,
         all_visits_fdr = same_direction & n_fdr == n_visits & n_visits == visits_analysed,
         all_visits_nominal = same_direction & n_nominal == n_visits & n_visits == visits_analysed,
         direction = if_else(mean_logFC > 0, "up", "down")) |>
  arrange(contrast, desc(all_visits_fdr), desc(all_visits_nominal), max_p)
save_csv(consistency, cfg, "visit_course", "consistency_across_visits.csv")
cons_sum <- consistency |> group_by(contrast) |>
  summarise(visits = max(n_visits),
            up_fdr = sum(all_visits_fdr & direction == "up"), down_fdr = sum(all_visits_fdr & direction == "down"),
            up_nominal = sum(all_visits_nominal & direction == "up"), down_nominal = sum(all_visits_nominal & direction == "down"),
            all_visits_fdr = sum(all_visits_fdr), all_visits_nominal = sum(all_visits_nominal))
save_csv(cons_sum, cfg, "visit_course", "consistency_summary.csv")
print(as.data.frame(cons_sum))

# ---- figures -----------------------------------------------------------------------------------------------
counts <- res |> group_by(visit, contrast, n_subjects) |>
  summarise(up = sum(significant & logFC > 0), down = sum(significant & logFC < 0), .groups = "drop") |>
  pivot_longer(c(up, down), names_to = "dir", values_to = "n")
p <- ggplot(counts, aes(visit, if_else(dir == "up", n, -n), fill = dir)) + geom_col() +
  geom_hline(yintercept = 0) + facet_wrap(~contrast) +
  scale_fill_manual(values = c(up = "firebrick", down = "steelblue")) +
  labs(title = sprintf("dISF: significant proteins per visit (FDR < %g)", fdr), x = NULL, y = "proteins (down < 0 < up)", fill = NULL)
save_plot(p, cfg, "visit_course", "significant_per_visit.png", width = 11, height = 4.5)

for (ct in unique(res$contrast)) {
  cs <- consistency |> filter(contrast == ct)
  sel <- cs |> filter(all_visits_fdr)
  tier <- "FDR at all visits"
  if (!nrow(sel)) { sel <- cs |> filter(all_visits_nominal); tier <- "p < 0.05 at all visits" }
  if (!nrow(sel)) next
  top <- sel |> slice_min(max_p, n = 24, with_ties = FALSE)
  d <- res |> filter(contrast == ct, OlinkID %in% top$OlinkID) |>
    mutate(Assay = factor(Assay, levels = unique(top$Assay)))
  p <- ggplot(d, aes(visit, logFC, group = 1)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, colour = "grey50") +
    geom_line(colour = "grey40") + geom_point(aes(colour = significant), size = 2) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"),
                        labels = c(`FALSE` = "n.s.", `TRUE` = sprintf("FDR < %g", fdr))) +
    facet_wrap(~Assay, scales = "free_y") +
    labs(title = sprintf("Time course of %s - proteins regulated at all visits (%s; top %d of %d)",
                         ct, tier, nrow(top), nrow(sel)), x = NULL, y = "difference in NPX (log2, 95% CI)", colour = NULL)
  save_plot(p, cfg, "visit_course", sprintf("time_course_effects_%s.png", ct), width = 13, height = 9)

  hm_ids <- sel |> slice_min(max_p, n = 60, with_ties = FALSE) |> arrange(mean_logFC)
  hm <- res |> filter(contrast == ct, OlinkID %in% hm_ids$OlinkID) |>
    mutate(Assay = factor(Assay, levels = unique(hm_ids$Assay)))
  p <- ggplot(hm, aes(visit, Assay, fill = logFC)) + geom_tile() +
    geom_text(aes(label = if_else(significant, "*", "")), size = 3) +
    scale_fill_gradient2(low = "steelblue", high = "firebrick") +
    labs(title = sprintf("%s across visits (%s; * FDR < %g)", ct, tier, fdr), x = NULL, y = NULL)
  save_plot(p, cfg, "visit_course", sprintf("time_course_heatmap_%s.png", ct), width = 7, height = 2 + 0.18 * nrow(hm_ids))
}

# NPX levels over visits for the top proteins regulated at all visits (lesion site vs healthy skin)
lv_ids <- consistency |> filter(contrast == "Lsite_vs_HC", all_visits_fdr | all_visits_nominal) |>
  slice_min(max_p, n = 12, with_ties = FALSE)
if (nrow(lv_ids)) {
  lv <- clean |> filter(matrix == "ISF", OlinkID %in% lv_ids$OlinkID) |>
    inner_join(info |> select(SampleID, site_cond), by = "SampleID") |> filter(!is.na(site_cond))
  hc_band <- lv |> filter(site_cond == "HC") |> group_by(Assay) |>
    summarise(m = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(n()))
  ad <- lv |> filter(site_cond != "HC") |> group_by(Assay, visit, site_cond) |>
    summarise(m = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop") |>
    mutate(site_cond = recode(site_cond, Lsite = "lesion site", NL = "non-lesional"))
  p <- ggplot(ad, aes(visit, m, colour = site_cond, group = site_cond)) +
    geom_rect(data = hc_band, aes(xmin = -Inf, xmax = Inf, ymin = m - se, ymax = m + se), inherit.aes = FALSE,
              fill = "grey85") +
    geom_hline(data = hc_band, aes(yintercept = m), colour = "grey50", linetype = 2) +
    geom_line() + geom_pointrange(aes(ymin = m - se, ymax = m + se), size = 0.3) +
    scale_colour_manual(values = c(`lesion site` = "firebrick", `non-lesional` = "steelblue")) +
    facet_wrap(~Assay, scales = "free_y") +
    labs(title = "dISF levels over visits (mean +/- SE; grey band: healthy skin)", x = NULL, y = "NPX", colour = NULL) +
    theme(legend.position = "bottom")
  save_plot(p, cfg, "visit_course", "time_course_levels.png", width = 12, height = 8)
}

writexl::write_xlsx(list(summary = cons_sum, consistency = consistency, lesion_site_state = comp,
                         per_visit_results = res |> select(-formula)),
                    out_path(cfg, "visit_course", "visit_course.xlsx"))
msg("Per-visit volcano plots: models/volcano/ISF_by_visit_V*.png; time course: visit_course/")
