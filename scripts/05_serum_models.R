# 05 - Serum models
# Out: output/models/Serum_results.csv, Serum_summary.csv, volcano plots
#
# AD vs controls is fitted twice: against the in-study healthy controls (adjusted for cohort)
# and against the LEIP biobank controls. Biobank serum is collected and stored differently,
# so a protein is most credible when both comparisons agree (column 'agree_both_controls').

source("R/utils.R")
source("R/models.R")
source("R/design.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
expr <- wide$Serum
clear_outputs(cfg, "models", "^Serum_")

info <- serum_design(meta) |> filter(SampleID %in% colnames(expr))
specs <- serum_specs(info)

res <- run_model_specs(specs, expr, info, cfg, "Serum", assay_map)

# agreement of the two AD-vs-control comparisons
if (nrow(res) && all(c("AD_vs_HC_in_study", "AD_vs_Biobank") %in% res$model)) {
  agree <- res |>
    filter(model %in% c("AD_vs_HC_in_study", "AD_vs_Biobank")) |>
    select(OlinkID, Assay, contrast, logFC, adj.P.Val, significant) |>
    pivot_wider(names_from = contrast, values_from = c(logFC, adj.P.Val, significant)) |>
    mutate(agree_both_controls = significant_AD_vs_HC & significant_AD_vs_Biobank &
             sign(logFC_AD_vs_HC) == sign(logFC_AD_vs_Biobank))
  save_csv(agree, cfg, "models", "Serum_AD_vs_controls_agreement.csv")
  msg("AD vs controls: %d proteins significant against both control groups (same direction)",
      sum(agree$agree_both_controls, na.rm = TRUE))
} else msg("AD vs controls agreement skipped: AD_vs_HC and/or AD_vs_Biobank not estimated.")

# ---- RELAD / RELAD2: all serum results in one workbook ------------------------------------------------------
clear_outputs(cfg, "relad")
rl_models <- c("RELAD_relapse", "RELAD_only_relapse", "RELAD2_only_relapse", "RELAD_relapse_unflagged", "RELAD_AD_vs_HC")
rl_ids <- info |> filter(cohort %in% c("RELAD", "RELAD2")) |> pull(SampleID)
rl_long <- clean |> filter(SampleID %in% rl_ids, matrix == "Serum")
rl_samples <- meta |> filter(cohort %in% c("RELAD", "RELAD2"), matrix == "Serum") |>
  select(SampleID, SubjectID, cohort, any_of("relad2_no"), group, relapse, relapse_raw, time_to_relapse, clinical_state,
         any_of("sex"), date, plate, well, volume_ul, flags) |>
  mutate(in_analysis = SampleID %in% rl_ids)
rl_means <- rl_long |>
  mutate(set = case_when(group == "HC" ~ paste(cohort, "healthy"), !is.na(relapse) ~ paste(cohort, "AD", relapse), TRUE ~ paste(cohort, "AD other"))) |>
  group_by(OlinkID, Assay, set) |>
  summarise(n = sum(!is.na(value)), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
            pct_above_LOD = round(100 * mean(!below_lod, na.rm = TRUE), 1), .groups = "drop") |>
  pivot_wider(names_from = set, values_from = c(n, mean, sd, pct_above_LOD), names_glue = "{set}: {.value}")
rl_res <- map(setNames(rl_models, rl_models), \(m) res |> filter(model == m) |> arrange(P.Value) |>
                select(Assay, OlinkID, contrast, logFC, AveExpr, t, P.Value, FDR = adj.P.Val, significant, n_samples, n_subjects, formula)) |>
  keep(\(d) nrow(d) > 0)
rl_sum <- map(rl_res, \(d) tibble(n_samples = d$n_samples[1], proteins = nrow(d), significant = sum(d$significant),
                                  up = sum(d$significant & d$logFC > 0), down = sum(d$significant & d$logFC < 0))) |>
  bind_rows(.id = "model")
rl_wide <- rl_long |> filter(keep %in% TRUE) |> select(SampleID, Assay, value) |>
  pivot_wider(names_from = Assay, values_from = value) |>
  right_join(rl_samples |> filter(in_analysis) |> select(SampleID, cohort, group, relapse, time_to_relapse), by = "SampleID") |>
  relocate(SampleID, cohort, group, relapse, time_to_relapse)
readme <- tibble(sheet = c("summary", "samples", names(rl_res), "group_means", "NPX_values"),
                 content = c("number of significant proteins per model (FDR < stats$fdr)",
                             "all RELAD and RELAD2 serum samples with their labels; in_analysis = passed sample QC",
                             paste("all proteins, model", names(rl_res), "(relapse_vs_non: relapse minus non-relapse; AD_vs_HC: AD minus healthy; log2 NPX)"),
                             "per protein: n, mean, SD and % above LOD per cohort and group (PCNormalizedNPX)",
                             "PCNormalizedNPX of every analysed serum protein for RELAD and RELAD2 samples (values below LOD as measured)"))
writexl::write_xlsx(c(list(README = readme, summary = rl_sum, samples = rl_samples), rl_res,
                      list(group_means = rl_means, NPX_values = rl_wide)),
                    out_path(cfg, "relad", "RELAD_RELAD2_serum_results.xlsx"))
msg("RELAD/RELAD2 workbook: %s", out_path(cfg, "relad", "RELAD_RELAD2_serum_results.xlsx"))
