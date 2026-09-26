# 05 - Serum models
# Out: output/models/Serum_results.csv, Serum_summary.csv, volcano plots
#
# AD vs controls is fitted twice: against the in-study healthy controls (adjusted for cohort)
# and against the LEIP biobank controls. Biobank serum is collected and stored differently,
# so a protein is most credible when both comparisons agree (column 'agree_both_controls').

source("R/utils.R")
source("R/models.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
expr <- wide$Serum
clear_outputs(cfg, "models", "^Serum_")

info <- meta |>
  filter(matrix == "Serum", SampleID %in% colnames(expr)) |>
  mutate(
    # one sample per person for cross-sectional comparisons: MicroAD baseline + single-sample cohorts
    cross_sectional = cohort != "MicroAD" | visit == "V1",
    status = case_when(group == "AD" ~ "AD", group == "HC" ~ "HC", group == "Biobank" ~ "Biobank"),
    relapse2 = if_else(relapse %in% c("relapse", "non-relapse"), str_replace(relapse, "-", "_"), NA_character_),
    flagged = !is.na(flags) & str_detect(flags, "Relapse|Group")
  )
ids <- \(...) info |> filter(...) |> pull(SampleID)

specs <- list(
  list(name = "AD_vs_HC_in_study",
       samples = ids(cross_sectional, status %in% c("AD", "HC")),
       formula = ~ 0 + status + cohort + plate,
       contrasts = c(AD_vs_HC = "statusAD - statusHC")),
  list(name = "AD_vs_Biobank",
       samples = ids(cross_sectional, status %in% c("AD", "Biobank")),
       formula = ~ 0 + status + plate,
       contrasts = c(AD_vs_Biobank = "statusAD - statusBiobank")),
  list(name = "HC_vs_Biobank",                  # pre-analytical / source check
       samples = ids(cross_sectional, status %in% c("HC", "Biobank")),
       formula = ~ 0 + status + plate,
       contrasts = c(HC_vs_Biobank = "statusHC - statusBiobank")),
  list(name = "MicroAD_active_vs_cleared",      # serum when the tracked lesion is active vs cleared
       samples = ids(cohort == "MicroAD", group == "AD", lesion_state %in% c("active", "cleared")),
       formula = ~ 0 + lesion_state + plate + (1 | SubjectID),
       contrasts = c(active_vs_cleared = "lesion_stateactive - lesion_statecleared")),
  list(name = "MicroAD_relapse",                # exploratory: cleared visits, relapsers vs non-relapsers
       samples = ids(cohort == "MicroAD", lesion_state == "cleared", !is.na(relapse2)),
       formula = ~ 0 + relapse2 + plate + (1 | SubjectID),
       contrasts = c(relapse_vs_non = "relapse2relapse - relapse2non_relapse")),
  list(name = "RELAD_relapse",
       samples = ids(cohort %in% c("RELAD", "RELAD2"), group == "AD", !is.na(relapse2)),
       formula = ~ 0 + relapse2 + cohort + plate,
       contrasts = c(relapse_vs_non = "relapse2relapse - relapse2non_relapse")),
  list(name = "RELAD_relapse_unflagged",        # sensitivity: without samples with conflicting labels
       samples = ids(cohort %in% c("RELAD", "RELAD2"), group == "AD", !is.na(relapse2), !flagged),
       formula = ~ 0 + relapse2 + cohort + plate,
       contrasts = c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"))
)

res <- run_model_specs(specs, expr, info, cfg, "Serum", assay_map)

# agreement of the two AD-vs-control comparisons
if (nrow(res) && all(c("AD_vs_HC", "AD_vs_Biobank") %in% res$contrast)) {
  agree <- res |>
    filter(contrast %in% c("AD_vs_HC", "AD_vs_Biobank")) |>
    select(OlinkID, Assay, contrast, logFC, adj.P.Val, significant) |>
    pivot_wider(names_from = contrast, values_from = c(logFC, adj.P.Val, significant)) |>
    mutate(agree_both_controls = significant_AD_vs_HC & significant_AD_vs_Biobank &
             sign(logFC_AD_vs_HC) == sign(logFC_AD_vs_Biobank))
  save_csv(agree, cfg, "models", "Serum_AD_vs_controls_agreement.csv")
  msg("AD vs controls: %d proteins significant against both control groups (same direction)",
      sum(agree$agree_both_controls, na.rm = TRUE))
} else msg("AD vs controls agreement skipped: AD_vs_HC and/or AD_vs_Biobank not estimated.")
