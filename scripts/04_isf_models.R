# 04 - ISF (dermal interstitial fluid) models, MicroAD
# Per protein; subject as random effect (dream), plate as fixed effect in every model.
# Out: output/models/ISF_results.csv, ISF_summary.csv, volcano plots

source("R/utils.R")
source("R/models.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
expr <- wide$ISF

info <- meta |>
  filter(matrix == "ISF", SampleID %in% colnames(expr)) |>
  mutate(
    cond = case_when(
      group == "AD" & state == "lesional"     ~ "AD_L",
      group == "AD" & state == "ex-lesional"  ~ "AD_xL",
      group == "AD" & state == "non-lesional" ~ "AD_NL",
      group == "HC"                           ~ "HC",
      group == "CPUO" & state == "lesional"   ~ "CPUO_L",
      group == "CPUO" & state == "non-lesional" ~ "CPUO_NL"),
    weeks = days_since_v1 / 7,
    relapse2 = if_else(relapse %in% c("relapse", "non-relapse"), str_replace(relapse, "-", "_"), NA_character_)
  )
ids <- \(...) info |> filter(...) |> pull(SampleID)

specs <- list(
  list(name = "states_all_visits",
       samples = ids(!is.na(cond)),
       formula = ~ 0 + cond + plate + (1 | SubjectID),
       contrasts = c(AD_L_vs_NL = "condAD_L - condAD_NL", AD_xL_vs_NL = "condAD_xL - condAD_NL",
                     AD_L_vs_xL = "condAD_L - condAD_xL", AD_NL_vs_HC = "condAD_NL - condHC",
                     AD_L_vs_HC = "condAD_L - condHC", AD_L_vs_CPUO_L = "condAD_L - condCPUO_L",
                     CPUO_L_vs_NL = "condCPUO_L - condCPUO_NL")),
  list(name = "baseline_V1",
       samples = ids(visit == "V1", !is.na(cond)),
       formula = ~ 0 + cond + plate + (1 | SubjectID),
       contrasts = c(AD_L_vs_NL = "condAD_L - condAD_NL", AD_NL_vs_HC = "condAD_NL - condHC",
                     AD_L_vs_HC = "condAD_L - condHC", CPUO_L_vs_NL = "condCPUO_L - condCPUO_NL")),
  list(name = "time_ex_lesional",           # change per week in cleared lesional skin
       samples = ids(cond == "AD_xL"),
       formula = ~ weeks + plate + (1 | SubjectID),
       contrasts = c(per_week_xL = "weeks")),
  list(name = "time_non_lesional",
       samples = ids(cond == "AD_NL", !dropout),
       formula = ~ weeks + plate + (1 | SubjectID),
       contrasts = c(per_week_NL = "weeks")),
  list(name = "relapse_ex_lesional",        # exploratory: cleared skin, relapsers vs non-relapsers
       samples = ids(cond == "AD_xL", !is.na(relapse2)),
       formula = ~ 0 + relapse2 + weeks + plate + (1 | SubjectID),
       contrasts = c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"))
)

res <- run_model_specs(specs, expr, info, cfg, "ISF", assay_map)

# Relapse on the within-visit difference ex-lesional minus non-lesional (same subject, visit, plate:
# removes plate and systemic day-to-day variation).
pairs <- info |>
  filter(group == "AD", !is.na(relapse2), cond %in% c("AD_xL", "AD_NL")) |>
  select(SubjectID, visit, cond, SampleID, relapse2, plate, weeks) |>
  pivot_wider(names_from = cond, values_from = SampleID) |>
  filter(!is.na(AD_xL), !is.na(AD_NL)) |>
  mutate(SampleID = paste(SubjectID, visit, sep = "_"))
if (nrow(pairs) >= 4) {
  delta <- expr[, pairs$AD_xL, drop = FALSE] - expr[, pairs$AD_NL, drop = FALSE]
  colnames(delta) <- pairs$SampleID
  rd <- fit_contrasts(delta, pairs, ~ 0 + relapse2 + weeks + (1 | SubjectID),
                      c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"),
                      "relapse_delta_xL_minus_NL", cfg$stats$min_group_n)
  if (!is.null(rd)) {
    rd <- annotate_results(rd, assay_map, cfg$stats$fdr)
    save_csv(rd, cfg, "models", "ISF_relapse_delta_results.csv")
    msg("relapse_delta_xL_minus_NL: %d pairs, %d significant", nrow(pairs), sum(rd$significant))
  }
}
