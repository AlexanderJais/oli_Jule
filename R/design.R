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
    list(name = "MicroAD_AD_vs_HC",              # MicroAD only (no RELAD/RELAD2): the serum reference for dISF comparisons
         samples = ids(cohort == "MicroAD", visit == "V1", status %in% c("AD", "HC")),
         formula = ~ 0 + status + plate,
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
    list(name = "RELAD_only_relapse",             # each cohort on its own
         samples = ids(cohort == "RELAD", group == "AD", !is.na(relapse2)),
         formula = ~ 0 + relapse2 + plate,
         contrasts = c(relapse_vs_non = "relapse2relapse - relapse2non_relapse")),
    list(name = "RELAD2_only_relapse",
         samples = ids(cohort == "RELAD2", group == "AD", !is.na(relapse2)),
         formula = ~ 0 + relapse2 + plate,
         contrasts = c(relapse_vs_non = "relapse2relapse - relapse2non_relapse")),
    list(name = "RELAD_AD_vs_HC",                 # RELAD + RELAD2 only
         samples = ids(cohort %in% c("RELAD", "RELAD2"), status %in% c("AD", "HC")),
         formula = ~ 0 + status + cohort + plate,
         contrasts = c(AD_vs_HC = "statusAD - statusHC")),
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

#' Per-visit ISF comparisons (step 13): at each visit, the tracked lesion site (lesional or, after
#' clearing, ex-lesional) vs healthy skin, vs non-lesional skin of the same patients, and non-lesional
#' vs healthy. Healthy controls (single visit) are the reference at every visit.
visit_specs <- function(info, min_subjects = 5) {
  info <- info |> mutate(site_cond = case_when(group == "AD" & site == "L" ~ "Lsite",
                                               group == "AD" & site == "NL" ~ "NL",
                                               group == "HC" ~ "HC"))
  visits <- info |> filter(group == "AD", !is.na(visit_num)) |> count(visit_num, SubjectID) |>
    count(visit_num, name = "n_subjects") |> filter(n_subjects >= min_subjects) |> pull(visit_num)
  hc <- info$SampleID[info$site_cond %in% "HC"]
  specs <- map(visits, \(v) list(
    name = paste0("V", v),
    samples = c(info$SampleID[info$group == "AD" & info$visit_num %in% v & !is.na(info$site_cond)], hc),
    formula = ~ 0 + site_cond + plate + (1 | SubjectID),
    contrasts = c(Lsite_vs_HC = "site_condLsite - site_condHC", Lsite_vs_NL = "site_condLsite - site_condNL",
                  NL_vs_HC = "site_condNL - site_condHC")))
  list(info = info, specs = specs)
}

#' All pre-specified single-protein tests (steps 12 and 15): the models of steps 04/05 plus the
#' xL - NL relapse model, for one value per sample (a protein or a score).
#' @param isf_info,serum_info isf_design()/serum_design() output joined with a `value` column
prespecified_tests <- function(isf_info, serum_info, min_group_n = 3) {
  run_specs <- \(specs, info) map(specs, \(sp) {
    test_single(info |> filter(SampleID %in% sp$samples), sp$formula, sp$contrasts, sp$name, min_group_n)
  }) |> bind_rows()
  tests <- bind_rows(
    run_specs(isf_specs(isf_info), isf_info) |> mutate(matrix = "dISF"),
    run_specs(serum_specs(serum_info), serum_info) |> mutate(matrix = "serum"))
  pairs <- isf_delta_pairs(isf_info)
  if (nrow(pairs) >= 4) {
    dl <- pairs |>
      mutate(value = isf_info$value[match(AD_xL, isf_info$SampleID)] - isf_info$value[match(AD_NL, isf_info$SampleID)])
    tests <- bind_rows(tests, test_single(dl, ~ 0 + relapse2 + weeks + (1 | SubjectID),
                                          c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"),
                                          "relapse_delta_xL_minus_NL", min_group_n) |> mutate(matrix = "dISF"))
  }
  tests
}
