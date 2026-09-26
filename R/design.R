# Study design shared by the proteome-wide models (04, 05) and the focus-protein analysis (12):
# group coding per sample and the model specifications. Change a model here and it changes
# everywhere.

#' ISF samples with model variables (cond = group x skin state, weeks since V1, relapse status).
isf_design <- function(meta) {
  meta |>
    filter(matrix == "ISF") |>
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
}

#' Serum samples with model variables.
serum_design <- function(meta) {
  meta |>
    filter(matrix == "Serum") |>
    mutate(
      # one sample per person for cross-sectional comparisons: MicroAD baseline + single-sample cohorts
      cross_sectional = cohort != "MicroAD" | visit == "V1",
      status = case_when(group == "AD" ~ "AD", group == "HC" ~ "HC", group == "Biobank" ~ "Biobank"),
      relapse2 = if_else(relapse %in% c("relapse", "non-relapse"), str_replace(relapse, "-", "_"), NA_character_),
      flagged = !is.na(flags) & str_detect(flags, "Relapse|Group")
    )
}

#' ISF model specifications: list(name, samples, formula, contrasts).
isf_specs <- function(info) {
  ids <- \(...) info |> filter(...) |> pull(SampleID)
  list(
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
}

#' Serum model specifications.
serum_specs <- function(info) {
  ids <- \(...) info |> filter(...) |> pull(SampleID)
  list(
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
}

#' Ex-lesional / non-lesional ISF pairs from the same AD visit (for the within-visit difference
#' xL - NL, which removes plate and systemic day-to-day variation). Pivot on subject + visit only,
#' so a pair is kept even if its two samples differ in plate or date; covariates come from the xL sample.
isf_delta_pairs <- function(info) {
  info |>
    filter(group == "AD", cond %in% c("AD_xL", "AD_NL")) |>
    select(SubjectID, visit, cond, SampleID) |>
    pivot_wider(id_cols = c(SubjectID, visit), names_from = cond, values_from = SampleID) |>
    filter(!is.na(AD_xL), !is.na(AD_NL)) |>
    left_join(info |> select(AD_xL = SampleID, relapse2, weeks), by = "AD_xL") |>
    filter(!is.na(relapse2)) |>
    mutate(SampleID = paste(SubjectID, visit, sep = "_"))
}
