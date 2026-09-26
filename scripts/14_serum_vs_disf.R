# 14 - Serum vs dISF, visit by visit: overlap and the additional information dISF provides
# The same question is asked in both matrices (MicroAD, matched visits):
#   AD vs healthy          dISF lesion site vs healthy skin  |  dISF non-lesional vs healthy skin
#                          serum AD vs healthy (MicroAD volunteers; one sample per person per visit)
#   relapse vs non-relapse dISF lesion site  |  dISF non-lesional  |  serum
#                          (only visits BEFORE the relapse, so it is a question about prediction)
# per visit (V1..V6, where group sizes allow) and pooled over all visits (subject as random effect).
# Significant protein sets are compared at two levels: FDR < stats$fdr (strict) and p < 0.05
# (nominal, exploratory - per-visit groups are small):
#   both (same / opposite direction), dISF only, serum only, and dISF-significant proteins that
#   are not measurable in serum at all ("additional via detectability").
# Out: output/serum_vs_disf/*  (overlap summary, protein lists, Venn diagrams, volcano/effect plots)

source("R/utils.R")
source("R/models.R")
source("R/design.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
clear_outputs(cfg, "serum_vs_disf")
fdr <- cfg$stats$fdr; mgn <- cfg$stats$min_group_n
shared <- intersect(rownames(wide$ISF), rownames(wide$Serum))

isf_info <- visit_specs(isf_design(meta) |> filter(SampleID %in% colnames(wide$ISF)))$info |>
  mutate(before_relapse = is.na(relapse_visit) | visit_num < relapse_visit)
ser_info <- serum_design(meta) |> filter(SampleID %in% colnames(wide$Serum), cohort == "MicroAD") |>
  mutate(before_relapse = is.na(relapse_visit) | visit_num < relapse_visit)
visits <- sort(unique(isf_info$visit_num[isf_info$group == "AD"]))

fit <- \(expr, info, ids, form, ct, name) {
  r <- fit_contrasts(expr, info |> filter(SampleID %in% ids), form, ct, name, mgn)
  if (is.null(r)) NULL else r |> left_join(assay_map, by = "OlinkID")
}
hc_isf <- isf_info$SampleID[isf_info$site_cond %in% "HC"]
hc_ser <- ser_info$SampleID[ser_info$status %in% "HC"]

# ---- the questions ----------------------------------------------------------------------------------------
questions <- list()
add_q <- function(question, visit, isf_site, isf_res, ser_res) {
  if (is.null(isf_res) || is.null(ser_res)) return(invisible())
  questions[[length(questions) + 1]] <<- list(question = question, visit = visit, isf_site = isf_site,
                                              isf = isf_res, serum = ser_res)
}
for (v in c(as.list(visits), list("all"))) {
  vlab <- if (identical(v, "all")) "all visits" else paste0("V", v)
  in_v <- \(info) if (identical(v, "all")) rep(TRUE, nrow(info)) else info$visit_num %in% v
  rand <- identical(v, "all")                 # pooled: several samples per person -> subject random effect
  msg("---- %s ----", vlab)

  # AD vs healthy
  ser_ids <- c(ser_info$SampleID[ser_info$status %in% "AD" & in_v(ser_info)], hc_ser)
  ser_ad <- fit(wide$Serum, ser_info, ser_ids,
                if (rand) ~ 0 + status + plate + (1 | SubjectID) else ~ 0 + status + plate,
                c(AD_vs_HC = "statusAD - statusHC"), paste("serum AD vs HC", vlab))
  for (st in c("Lsite", "NL")) {
    ids <- c(isf_info$SampleID[isf_info$site_cond %in% st & in_v(isf_info)], hc_isf)
    isf_ad <- fit(wide$ISF, isf_info, ids,
                  if (rand) ~ 0 + site_cond + plate + (1 | SubjectID) else ~ 0 + site_cond + plate,
                  setNames(sprintf("site_cond%s - site_condHC", st), "AD_vs_HC"), paste("dISF", st, "vs HC", vlab))
    add_q("AD vs healthy", vlab, if (st == "Lsite") "lesion site" else "non-lesional", isf_ad, ser_ad)
  }

  # relapse vs non-relapse (visits before relapse)
  ser_ids <- ser_info$SampleID[ser_info$group %in% "AD" & !is.na(ser_info$relapse2) & ser_info$before_relapse & in_v(ser_info)]
  ser_rl <- fit(wide$Serum, ser_info, ser_ids,
                if (rand) ~ 0 + relapse2 + plate + (1 | SubjectID) else ~ 0 + relapse2 + plate,
                c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"), paste("serum relapse", vlab))
  for (st in c("Lsite", "NL")) {
    ids <- isf_info$SampleID[isf_info$site_cond %in% st & !is.na(isf_info$relapse2) & isf_info$before_relapse & in_v(isf_info)]
    isf_rl <- fit(wide$ISF, isf_info, ids,
                  if (rand) ~ 0 + relapse2 + plate + (1 | SubjectID) else ~ 0 + relapse2 + plate,
                  c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"), paste("dISF", st, "relapse", vlab))
    add_q("relapse vs non-relapse", vlab, if (st == "Lsite") "lesion site" else "non-lesional", isf_rl, ser_rl)
  }
}
if (!length(questions)) stop("No question could be fitted in both matrices.")

# ---- overlap ------------------------------------------------------------------------------------------------------
categorise <- function(q, tier) {
  sig <- \(r) if (tier == "FDR") r$adj.P.Val < fdr else r$P.Value < 0.05
  i <- q$isf |> transmute(OlinkID, Assay, isf_logFC = logFC, isf_p = P.Value, isf_fdr = adj.P.Val, isf_sig = sig(q$isf))
  s <- q$serum |> transmute(OlinkID, ser_logFC = logFC, ser_p = P.Value, ser_fdr = adj.P.Val, ser_sig = sig(q$serum))
  full_join(i, s, by = "OlinkID") |>
    mutate(Assay = coalesce(Assay, assay_map$Assay[match(OlinkID, assay_map$OlinkID)]),
           in_both_panels = OlinkID %in% shared,
           category = case_when(
             isf_sig %in% TRUE & ser_sig %in% TRUE & sign(isf_logFC) == sign(ser_logFC) ~ "both, same direction",
             isf_sig %in% TRUE & ser_sig %in% TRUE ~ "both, opposite direction",
             isf_sig %in% TRUE & !in_both_panels ~ "dISF only - not measurable in serum",
             isf_sig %in% TRUE ~ "dISF only",
             ser_sig %in% TRUE & !in_both_panels ~ "serum only - not measurable in dISF",
             ser_sig %in% TRUE ~ "serum only",
             TRUE ~ "not significant"),
           question = q$question, visit = q$visit, isf_site = q$isf_site, tier = tier)
}
cats <- map(questions, \(q) bind_rows(categorise(q, "FDR"), categorise(q, "nominal"))) |> bind_rows() |>
  mutate(visit = factor(visit, levels = c(paste0("V", 1:12), "all visits")) |> droplevels())
design_n <- map(questions, \(q) tibble(question = q$question, visit = q$visit, isf_site = q$isf_site,
                                       isf_samples = q$isf$n_samples[1], serum_samples = q$serum$n_samples[1])) |> bind_rows()

summ <- cats |>
  group_by(question, visit, isf_site, tier) |>
  summarise(dISF_significant = sum(isf_sig %in% TRUE), serum_significant = sum(ser_sig %in% TRUE),
            both_same = sum(category == "both, same direction"), both_opposite = sum(category == "both, opposite direction"),
            dISF_only = sum(category == "dISF only"), dISF_only_not_in_serum = sum(category == "dISF only - not measurable in serum"),
            serum_only = sum(category == "serum only"), serum_only_not_in_dISF = sum(category == "serum only - not measurable in dISF"),
            .groups = "drop") |>
  mutate(overlap_pct_of_serum = round(100 * (both_same + both_opposite) / pmax(serum_significant, 1), 1),
         additional_from_dISF = dISF_only + dISF_only_not_in_serum, visit = as.character(visit)) |>
  left_join(design_n, by = c("question", "visit", "isf_site")) |>
  mutate(visit = factor(visit, levels = c(paste0("V", 1:12), "all visits")) |> droplevels()) |>
  arrange(tier, question, isf_site, visit)
save_csv(summ, cfg, "serum_vs_disf", "overlap_summary.csv")
save_csv(cats |> filter(category != "not significant"), cfg, "serum_vs_disf", "overlap_protein_lists.csv")
print(summ |> filter(tier == "nominal") |> select(question, visit, isf_site, dISF_significant, serum_significant,
                                                   both_same, dISF_only, dISF_only_not_in_serum, serum_only) |> as.data.frame())

# ---- figures --------------------------------------------------------------------------------------------------------
cat_cols <- c(`both, same direction` = "purple3", `both, opposite direction` = "orange3", `dISF only` = "firebrick",
              `dISF only - not measurable in serum` = "darkred", `serum only` = "steelblue",
              `serum only - not measurable in dISF` = "navy", `not significant` = "grey80")

# stacked bars: what dISF adds, per visit
bars <- summ |> select(question, visit, isf_site, tier, both_same, both_opposite, dISF_only, dISF_only_not_in_serum,
                       serum_only, serum_only_not_in_dISF) |>
  pivot_longer(c(both_same, both_opposite, dISF_only, dISF_only_not_in_serum, serum_only, serum_only_not_in_dISF),
               names_to = "cat", values_to = "n") |>
  mutate(cat = recode(cat, both_same = "both, same direction", both_opposite = "both, opposite direction", dISF_only = "dISF only",
                      dISF_only_not_in_serum = "dISF only - not measurable in serum", serum_only = "serum only",
                      serum_only_not_in_dISF = "serum only - not measurable in dISF"),
         cat = factor(cat, levels = setdiff(names(cat_cols), "not significant")))
for (tr in c("FDR", "nominal")) {
  p <- ggplot(bars |> filter(tier == tr), aes(visit, n, fill = cat)) + geom_col() +
    facet_grid(question ~ isf_site, scales = "free_y") + scale_fill_manual(values = cat_cols) +
    labs(title = sprintf("Serum vs dISF: significant proteins per visit (%s)", if (tr == "FDR") sprintf("FDR < %g", fdr) else "p < 0.05, exploratory"),
         x = NULL, y = "proteins", fill = NULL) + theme(legend.position = "bottom") + guides(fill = guide_legend(nrow = 2))
  save_plot(p, cfg, "serum_vs_disf", sprintf("overlap_bars_%s.png", tr), width = 12, height = 7)
}

# Venn diagrams (drawn with ggplot): left = dISF, right = serum
circ <- \(cx) tibble(x = cx + cos(seq(0, 2 * pi, length.out = 120)), y = sin(seq(0, 2 * pi, length.out = 120)))
circles <- bind_rows(circ(-0.55) |> mutate(set = "dISF"), circ(0.55) |> mutate(set = "serum"))
for (tr in c("FDR", "nominal")) for (qq in unique(summ$question)) {
  d <- summ |> filter(tier == tr, question == qq)
  if (!nrow(d)) next
  lab <- d |> transmute(visit, isf_site,
                        left = sprintf("%d", dISF_only + dISF_only_not_in_serum),
                        left2 = if_else(dISF_only_not_in_serum > 0, sprintf("(%d not in\nserum panel)", dISF_only_not_in_serum), ""),
                        mid = sprintf("%d", both_same + both_opposite),
                        mid2 = if_else(both_opposite > 0, sprintf("(%d opposite)", both_opposite), ""),
                        right = sprintf("%d", serum_only + serum_only_not_in_dISF))
  circ_all <- tidyr::expand_grid(d |> distinct(visit, isf_site), circles)   # expand_grid keeps the point order
  lv <- c(paste0("V", 1:12), "all visits")
  lab <- lab |> mutate(visit = factor(as.character(visit), levels = lv) |> droplevels())
  circ_all <- circ_all |> mutate(visit = factor(as.character(visit), levels = levels(lab$visit)))
  p <- ggplot() +
    geom_polygon(data = circ_all, aes(x, y, group = set, fill = set), alpha = 0.3, colour = "grey30") +
    geom_text(data = lab, aes(x = -1.05, y = 0.1, label = left), size = 4.5) +
    geom_text(data = lab, aes(x = -1.05, y = -0.35, label = left2), size = 2.3) +
    geom_text(data = lab, aes(x = 0, y = 0.1, label = mid), size = 4.5) +
    geom_text(data = lab, aes(x = 0, y = -0.3, label = mid2), size = 2.3) +
    geom_text(data = lab, aes(x = 1.05, y = 0.1, label = right), size = 4.5) +
    scale_fill_manual(values = c(dISF = "firebrick", serum = "steelblue")) +
    coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.1, 1.1)) + facet_grid(isf_site ~ visit) +
    labs(title = sprintf("%s: significant proteins in dISF and serum (%s)", qq,
                         if (tr == "FDR") sprintf("FDR < %g", fdr) else "p < 0.05, exploratory"),
         fill = NULL, x = NULL, y = NULL) +
    theme_void() + theme(legend.position = "bottom", strip.text = element_text(size = 9), plot.title = element_text(size = 12))
  save_plot(p, cfg, "serum_vs_disf", sprintf("venn_%s_%s.png", str_replace_all(qq, "[^A-Za-z]+", "_"), tr),
            width = 2 + 2.2 * n_distinct(d$visit), height = 5.5)
}

# dISF volcano coloured by overlap with serum, and dISF vs serum effects (nominal tier)
for (qq in unique(cats$question)) {
  d <- cats |> filter(tier == "nominal", question == qq, !is.na(isf_p))
  p <- ggplot(d, aes(isf_logFC, -log10(isf_p), colour = category)) +
    geom_point(data = \(x) filter(x, category == "not significant"), size = 0.5) +
    geom_point(data = \(x) filter(x, category != "not significant"), size = 0.9) +
    scale_colour_manual(values = cat_cols) + facet_grid(isf_site ~ visit) +
    labs(title = sprintf("%s - dISF volcano, coloured by what serum shows (p < 0.05)", qq),
         x = "dISF effect (log2)", y = "-log10 p (dISF)", colour = NULL) +
    theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2)))
  save_plot(p, cfg, "serum_vs_disf", sprintf("volcano_dISF_coloured_%s.png", str_replace_all(qq, "[^A-Za-z]+", "_")),
            width = 3 + 2.3 * n_distinct(d$visit), height = 7)
  e <- cats |> filter(tier == "nominal", question == qq, in_both_panels, !is.na(isf_logFC), !is.na(ser_logFC))
  p <- ggplot(e, aes(isf_logFC, ser_logFC, colour = category)) +
    geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(data = \(x) filter(x, category == "not significant"), size = 0.4) +
    geom_point(data = \(x) filter(x, category != "not significant"), size = 0.9) +
    scale_colour_manual(values = cat_cols) + facet_grid(isf_site ~ visit) +
    labs(title = sprintf("%s - effect in dISF vs serum (proteins measured in both; p < 0.05)", qq),
         x = "dISF effect (log2)", y = "serum effect (log2)", colour = NULL) +
    theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2)))
  save_plot(p, cfg, "serum_vs_disf", sprintf("effects_dISF_vs_serum_%s.png", str_replace_all(qq, "[^A-Za-z]+", "_")),
            width = 3 + 2.3 * n_distinct(e$visit), height = 7)
}

writexl::write_xlsx(list(summary = summ, protein_lists = cats |> filter(category != "not significant"),
                         all_results = cats),
                    out_path(cfg, "serum_vs_disf", "serum_vs_disf.xlsx"))
msg("Serum vs dISF comparison: %s", file.path(cfg$paths$output, "serum_vs_disf"))
