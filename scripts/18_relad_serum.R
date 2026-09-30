# 18 - RELAD and RELAD2 serum: all results in one Excel workbook
# One serum sample per person (AD patients in remission with / without later relapse, active AD,
# healthy controls). Models per protein (limma, adjusted for plate; pooled models also for cohort):
#   relapse vs non-relapse   RELAD + RELAD2 pooled, without samples with conflicting labels, each cohort alone
#   AD vs healthy            pooled and per cohort
#   time to relapse          RELAD2 relapse < 1 week vs > 1 week; RELAD TimeToRelapse as a number
#   active AD vs remission   RELAD2 (ClinicalStateSkin)
# plus detection per cohort and group, the key proteins, TNFRSF9 correlations, and the NPX values.
# Out: output/relad/RELAD_RELAD2_serum_results.xlsx (+ volcano plots, CSV of all results)

source("R/utils.R")
source("R/models.R")
source("R/design.R")
source("R/correlation.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
raw   <- read_step(cfg, "data", "npx_all_samples.rds", step = "scripts/02_import_qc.R")
clear_outputs(cfg, "relad")
fdr <- cfg$stats$fdr; mgn <- cfg$stats$min_group_n
assay_map <- clean |> distinct(OlinkID, Assay, UniProt)
cohorts <- c("RELAD", "RELAD2")

info <- serum_design(meta) |>
  filter(cohort %in% cohorts, SampleID %in% colnames(wide$Serum)) |>
  mutate(ttr_number = suppressWarnings(as.numeric(time_to_relapse)),
         ttr_class = case_when(time_to_relapse == "<1w" ~ "lt1w", time_to_relapse == ">1w" ~ "gt1w"),
         activity = case_when(clinical_state == "active AD" ~ "active", clinical_state == "remission" ~ "remission"))
if (!nrow(info)) stop("No RELAD / RELAD2 serum samples in the data.")
msg("RELAD / RELAD2 serum samples: %s", paste(names(table(info$cohort)), table(info$cohort), collapse = ", "))

# ---- models ------------------------------------------------------------------------------------------------
ids <- \(...) info |> filter(...) |> pull(SampleID)
rl <- c(relapse_vs_non = "relapse2relapse - relapse2non_relapse")
adhc <- c(AD_vs_HC = "statusAD - statusHC")
specs <- list(
  list(name = "relapse RELAD+RELAD2", samples = ids(group == "AD", !is.na(relapse2)),
       formula = ~ 0 + relapse2 + cohort + plate, contrasts = rl, groups = "relapse2"),
  list(name = "relapse RELAD+RELAD2 unflagged", samples = ids(group == "AD", !is.na(relapse2), !flagged),
       formula = ~ 0 + relapse2 + cohort + plate, contrasts = rl, groups = "relapse2"),
  list(name = "relapse RELAD", samples = ids(cohort == "RELAD", group == "AD", !is.na(relapse2)),
       formula = ~ 0 + relapse2 + plate, contrasts = rl, groups = "relapse2"),
  list(name = "relapse RELAD2", samples = ids(cohort == "RELAD2", group == "AD", !is.na(relapse2)),
       formula = ~ 0 + relapse2 + plate, contrasts = rl, groups = "relapse2"),
  list(name = "AD vs healthy RELAD+RELAD2", samples = ids(status %in% c("AD", "HC"), !dropout),
       formula = ~ 0 + status + cohort + plate, contrasts = adhc, groups = "status"),
  list(name = "AD vs healthy RELAD", samples = ids(cohort == "RELAD", status %in% c("AD", "HC"), !dropout),
       formula = ~ 0 + status + plate, contrasts = adhc, groups = "status"),
  list(name = "AD vs healthy RELAD2", samples = ids(cohort == "RELAD2", status %in% c("AD", "HC"), !dropout),
       formula = ~ 0 + status + plate, contrasts = adhc, groups = "status"),
  list(name = "RELAD2 relapse <1w vs >1w", samples = ids(cohort == "RELAD2", relapse %in% "relapse", !is.na(ttr_class)),
       formula = ~ 0 + ttr_class + plate, contrasts = c(lt1w_vs_gt1w = "ttr_classlt1w - ttr_classgt1w"), groups = "ttr_class"),
  list(name = "RELAD time to relapse", samples = ids(cohort == "RELAD", relapse %in% "relapse", !is.na(ttr_number)),
       formula = ~ ttr_number + plate, contrasts = c(per_unit_TimeToRelapse = "ttr_number"), groups = NULL),
  list(name = "RELAD2 active AD vs remission", samples = ids(cohort == "RELAD2", !is.na(activity)),
       formula = ~ 0 + activity + plate, contrasts = c(active_vs_remission = "activityactive - activityremission"),
       groups = "activity")
)

res <- map(specs, \(sp) {
  msg("Model %s: %d samples", sp$name, length(sp$samples))
  fit_contrasts(wide$Serum, info |> filter(SampleID %in% sp$samples), sp$formula, sp$contrasts, sp$name, mgn)
}) |> bind_rows()
if (!nrow(res)) stop("No RELAD / RELAD2 model could be fitted.")
# approximate 95% CI from the moderated standard error (logFC / t)
res <- annotate_results(res, assay_map |> distinct(OlinkID, Assay), fdr) |>
  mutate(se = abs(logFC / t), ci_low = logFC - 1.96 * se, ci_high = logFC + 1.96 * se) |> select(-se)

groups_n <- map(specs, \(sp) {
  d <- info |> filter(SampleID %in% sp$samples)
  g <- if (is.null(sp$groups)) "all" else paste(names(table(d[[sp$groups]])), table(d[[sp$groups]]), sep = " = ", collapse = ", ")
  tibble(model = sp$name, samples = nrow(d), groups = g, formula = paste(deparse(sp$formula), collapse = ""))
}) |> bind_rows()
summ <- res |> group_by(model, contrast) |>
  summarise(proteins = n(), significant_FDR = sum(significant), up = sum(significant & logFC > 0),
            down = sum(significant & logFC < 0), nominal_p05 = sum(P.Value < 0.05),
            top_10 = paste(head(Assay[order(P.Value)], 10), collapse = ", "), .groups = "drop") |>
  left_join(groups_n, by = "model") |>
  mutate(model = factor(model, levels = map_chr(specs, "name"))) |> arrange(model) |> mutate(model = as.character(model))
print(summ |> select(model, contrast, samples, significant_FDR, nominal_p05) |> as.data.frame())
save_csv(res, cfg, "relad", "RELAD_RELAD2_serum_results.csv")
for (mdl in unique(res$model))
  save_plot(volcano(res |> filter(model == mdl), paste("Serum", mdl), fdr), cfg, "relad", "volcano",
            paste0(str_replace_all(mdl, "[^A-Za-z0-9]+", "_"), ".png"), width = 7, height = 6)

# ---- detection per cohort and group -------------------------------------------------------------------------
sd_ <- clean |> filter(cohort %in% cohorts, matrix == "Serum") |>
  mutate(grp = case_when(group == "HC" ~ "healthy", relapse %in% c("relapse", "non-relapse") ~ relapse,
                         relapse %in% "active" | clinical_state %in% "active AD" ~ "active AD", TRUE ~ "AD other"))
detection <- sd_ |> group_by(OlinkID, Assay, cohort, grp) |>
  summarise(n = sum(!is.na(value)), pct_above_LOD = round(100 * mean(!below_lod, na.rm = TRUE), 1), .groups = "drop") |>
  mutate(col = paste(cohort, grp, sep = ": ")) |> select(OlinkID, Assay, col, pct_above_LOD) |>
  pivot_wider(names_from = col, values_from = pct_above_LOD) |>
  left_join(sd_ |> group_by(OlinkID) |> summarise(pct_above_LOD_all = round(100 * mean(!below_lod, na.rm = TRUE), 1),
                                                  LOD_median = median(LOD, na.rm = TRUE), analysed = any(keep %in% TRUE)),
            by = "OlinkID") |>
  relocate(pct_above_LOD_all, LOD_median, analysed, .after = Assay) |> arrange(Assay)

# ---- key proteins: TNFRSF9, TNFSF9, mast cell markers, TNFRSF9 partners ---------------------------------------
tn <- cfg$tnfrsf9 %||% list()
key_names <- unique(c(unlist(cfg$key_questions$proteins %||% c("TNFRSF9", "TNFSF9")), tn$protein %||% "TNFRSF9",
                      unlist(tn$partners), unlist(cfg$key_questions$mast_cell_markers), unlist(cfg$focus_proteins)))
key_oid <- setNames(map_chr(key_names, \(n) find_assay(assay_map, n)), key_names)
key_res <- res |> filter(OlinkID %in% key_oid) |> arrange(match(OlinkID, key_oid), model)
key_missing <- names(key_oid)[is.na(key_oid)]

# TNFRSF9 vs its partners in RELAD / RELAD2 serum (Spearman, values as measured; one sample per person)
tp_oid <- find_assay(assay_map, tn$protein %||% "TNFRSF9")
tnf_corr <- NULL
if (!is.na(tp_oid)) {
  v <- clean |> filter(matrix == "Serum", cohort %in% cohorts) |> select(SampleID, OlinkID, value)
  x <- v |> filter(OlinkID == tp_oid) |> select(SampleID, x = value)
  strata <- list(`all RELAD + RELAD2` = \(d) d, `AD` = \(d) filter(d, group == "AD"), `healthy` = \(d) filter(d, group == "HC"),
                 `relapse` = \(d) filter(d, relapse %in% "relapse"), `non-relapse` = \(d) filter(d, relapse %in% "non-relapse"),
                 RELAD = \(d) filter(d, cohort == "RELAD"), RELAD2 = \(d) filter(d, cohort == "RELAD2"))
  tnf_corr <- map(setdiff(key_oid[!is.na(key_oid)], tp_oid), \(o) {
    d <- info |> select(SampleID, cohort, group, relapse) |> inner_join(x, by = "SampleID") |>
      inner_join(v |> filter(OlinkID == o) |> select(SampleID, y = value), by = "SampleID")
    imap(strata, \(f, nm) { dd <- f(d); spearman_row(dd$x, dd$y) |> mutate(stratum = nm, .before = 1) }) |> bind_rows() |>
      mutate(partner = assay_map$Assay[match(o, assay_map$OlinkID)], .before = 1)
  }) |> bind_rows() |> mutate(protein = tn$protein %||% "TNFRSF9", .before = 1)
}

# ---- samples and values ------------------------------------------------------------------------------------------
samples <- meta |> filter(matrix == "Serum", cohort %in% cohorts) |>
  select(SampleID, SubjectID, SampleName, cohort, group, relapse, relapse_raw, time_to_relapse, clinical_state, dropout,
         plate, well, volume_ul, any_of("sex"), flags) |>
  mutate(in_analysis = SampleID %in% info$SampleID) |> arrange(cohort, SubjectID)
vals <- raw |> filter(SampleID %in% samples$SampleID)
dup_assay <- vals |> distinct(OlinkID, Assay) |> count(Assay) |> filter(n > 1) |> pull(Assay)
vals <- vals |> mutate(col = if_else(Assay %in% dup_assay, paste(Assay, OlinkID, sep = "_"), Assay))
npx_wide <- vals |> select(SampleID, col, PCNormalizedNPX) |>
  pivot_wider(names_from = col, values_from = PCNormalizedNPX)
npx_wide <- samples |> select(SampleID, SubjectID, cohort, group, relapse, time_to_relapse) |>
  inner_join(npx_wide |> select(SampleID, all_of(sort(setdiff(names(npx_wide), "SampleID")))), by = "SampleID")

# ---- workbook ----------------------------------------------------------------------------------------------------------
readme <- tibble(sheet = c("README", "summary", "all_results", "<one sheet per model>", "key_proteins", "TNFRSF9_correlations",
                           "detection", "samples", "NPX_values"),
                 content = c(
                   "This workbook: all RELAD and RELAD2 serum results of the O-MicroAD Olink Explore HT run (step 18).",
                   "One row per model: samples and group sizes, formula, number of significant proteins (FDR < 0.05 and p < 0.05), top 10.",
                   "Every protein x every model: logFC (difference in NPX, log2; first group minus second), 95% CI, t, p, FDR (BH within model).",
                   "The same results per model, sorted by p.",
                   sprintf("The pre-specified proteins (TNFRSF9, TNFSF9, mast cell markers, TNFRSF9 partners, focus proteins).%s",
                           if (length(key_missing)) paste(" Not measured / not analysed in serum:", paste(key_missing, collapse = ", ")) else ""),
                   "Spearman correlation of TNFRSF9 with its partner proteins in RELAD / RELAD2 serum (one sample per person).",
                   "% of samples above LOD per protein, per cohort and group. 'analysed' = passed the serum detection filter (step 02).",
                   "Sample information (manifest Ver2 codes harmonised). in_analysis = passed QC and in the serum data.",
                   "PCNormalizedNPX (the values analysed) for all proteins and RELAD / RELAD2 samples, below-LOD values as measured."))
model_sheets <- map(unique(res$model), \(m) res |> filter(model == m) |> arrange(P.Value) |>
                      select(Assay, OlinkID, contrast, logFC, ci_low, ci_high, t, P.Value, adj.P.Val, significant, AveExpr, n_samples))
names(model_sheets) <- str_trunc(str_replace_all(unique(res$model), "[^A-Za-z0-9<>+ ]", ""), 31, ellipsis = "")
sheets <- c(list(README = readme, summary = summ, all_results = res |> arrange(model, P.Value)),
            model_sheets,
            list(key_proteins = key_res, TNFRSF9_correlations = tnf_corr, detection = detection, samples = samples,
                 NPX_values = npx_wide))
sheets <- Filter(\(x) !is.null(x) && nrow(x), sheets)
names(sheets) <- make.unique(str_replace_all(names(sheets), "[\\[\\]\\*\\?/\\\\:]", ""), sep = "_")
writexl::write_xlsx(sheets, out_path(cfg, "relad", "RELAD_RELAD2_serum_results.xlsx"))
msg("RELAD / RELAD2 serum results: %s", file.path(cfg$paths$output, "relad", "RELAD_RELAD2_serum_results.xlsx"))
