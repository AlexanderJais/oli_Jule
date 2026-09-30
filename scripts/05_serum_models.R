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
if (nrow(res) && all(c("AD_vs_HC", "AD_vs_Biobank") %in% res$contrast)) {
  agree <- res |>
    filter(model %in% c("AD_vs_HC_in_study", "AD_vs_Biobank")) |>   # not AD_vs_HC_MicroAD (same contrast name)
    select(OlinkID, Assay, contrast, logFC, adj.P.Val, significant) |>
    pivot_wider(names_from = contrast, values_from = c(logFC, adj.P.Val, significant)) |>
    mutate(agree_both_controls = significant_AD_vs_HC & significant_AD_vs_Biobank &
             sign(logFC_AD_vs_HC) == sign(logFC_AD_vs_Biobank))
  save_csv(agree, cfg, "models", "Serum_AD_vs_controls_agreement.csv")
  msg("AD vs controls: %d proteins significant against both control groups (same direction)",
      sum(agree$agree_both_controls, na.rm = TRUE))
} else msg("AD vs controls agreement skipped: AD_vs_HC and/or AD_vs_Biobank not estimated.")
