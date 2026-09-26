# End-to-end test on simulated data: simulate -> run all steps -> check that the built-in
# effects are found and that false positives stay rare.
#   Rscript tests/test_pipeline.R
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). ",
       "In RStudio open oli_Jule.Rproj, or run setwd(\"path/to/oli_Jule\") first. Current: ", getwd())
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
# fully synthetic: no study data needed
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
stopifnot(system2(rscript, "run_all.R") == 0)

out   <- "output_sim"
truth <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
rd <- \(f) read_csv(file.path(out, f), show_col_types = FALSE) |> left_join(truth |> select(OlinkID, role), by = "OlinkID")
hits <- \(r, mdl, ct, rl) with(r |> filter(model == mdl, contrast == ct), c(found = sum(significant & role == rl), of = sum(role == rl)))
fp   <- \(r, mdl, ct) with(r |> filter(model == mdl, contrast == ct), mean(significant[role == "null"]))

isf <- rd("models/ISF_results.csv"); ser <- rd("models/Serum_results.csv")
delta <- rd("models/ISF_relapse_delta_results.csv") |> mutate(model = "delta")
cor <- rd("isf_serum/isf_serum_correlation.csv")
checks <- list(
  "ISF lesional vs non-lesional"   = hits(isf, "states_all_visits", "AD_L_vs_NL", "ISF_lesional"),
  "ISF AD non-lesional vs healthy" = hits(isf, "states_all_visits", "AD_NL_vs_HC", "ISF_AD_vs_HC"),
  "ISF relapse (xL - NL delta)"    = hits(delta, "delta", "relapse_vs_non", "ISF_relapse"),
  "Serum AD vs in-study HC"        = hits(ser, "AD_vs_HC_in_study", "AD_vs_HC", "Serum_AD_vs_HC"),
  "Serum RELAD relapse"            = hits(ser, "RELAD_relapse", "relapse_vs_non", "Serum_relapse"),
  "ISF-serum coupling (L site)"    = with(cor |> filter(site == "L"), c(found = sum(significant & role == "ISF_serum_coupled"), of = sum(role == "ISF_serum_coupled")))
)
res <- tibble(check = names(checks), found = map_dbl(checks, "found"), of = map_dbl(checks, "of")) |>
  mutate(power = found / of)
print(res)
fpr <- c(ISF = fp(isf, "states_all_visits", "AD_L_vs_NL"), Serum = fp(ser, "AD_vs_HC_in_study", "AD_vs_HC"))
print(round(fpr, 3))

biobank <- ser |> filter(model == "AD_vs_HC_in_study", role == "Biobank_shift")
leip <- read_csv(file.path(out, "leip_reference/leip_clinical_associations.csv"), show_col_types = FALSE)
th2 <- read_csv(file.path(out, "enrichment/gsea_results.csv"), show_col_types = FALSE) |>
  filter(pathway == "CUSTOM_AD_TH2_AXIS", model %in% c("states_all_visits", "baseline_V1"), contrast == "AD_L_vs_NL")

prof <- read_csv(file.path(out, "isf_profile/isf_detection_profile.csv"), show_col_types = FALSE) |>
  left_join(truth |> select(OlinkID, role), by = "OlinkID")
enr <- read_csv(file.path(out, "matrix_comparison/relative_enrichment.csv"), show_col_types = FALSE) |>
  left_join(truth |> select(OlinkID, isf_offset), by = "OlinkID")
enr_r <- enr |> group_by(model) |> summarise(r = cor(rel_log2_isf_vs_serum, isf_offset))
traj <- read_csv(file.path(out, "trajectories/trajectory_results.csv"), show_col_types = FALSE) |>
  left_join(truth |> select(OlinkID, role), by = "OlinkID")
easi <- traj |> filter(model == "ISF lesional site: EASI")
print(enr_r)
cd137 <- readxl::read_excel(file.path(out, "focus/TNFRSF9/TNFRSF9_report.xlsx"), sheet = "prespecified_tests")
cd137_LvNL <- cd137 |> filter(model == "states_all_visits", contrast == "AD_L_vs_NL")

vc <- read_csv(file.path(out, "visit_course/consistency_across_visits.csv"), show_col_types = FALSE) |>
  left_join(truth |> select(OlinkID, role), by = "OlinkID") |> filter(contrast == "Lsite_vs_HC")
stopifnot(
  "no protein regulated at all visits found (step 13)" = sum(vc$all_visits_nominal & vc$role != "null") >= 3,
  "executive summary PDF missing (step 14)" = file.exists(file.path(out, "Executive_summary.pdf")) &&
    file.size(file.path(out, "Executive_summary.pdf")) > 20000,
  "lesion-restricted proteins not found (step 09)" =
    sum(prof$lesion_restricted & prof$role == "ISF_lesion_restricted", na.rm = TRUE) >= 4,
  "relative dISF/serum enrichment not recovered (step 10)" = all(enr_r$r > 0.9),
  "severity-linked proteins not found (step 11)" = mean(easi$significant[easi$role == "ISF_lesional"]) >= 0.7,
  "too many false trajectory hits (step 11)" =
    sum(traj$significant & traj$role == "null") <= max(3, 0.15 * sum(traj$significant)),
  "CD137 (TNFRSF9) lesional effect not found in the dedicated analysis (step 12)" =
    nrow(cd137_LvNL) == 1 && cd137_LvNL$p < 0.05 && cd137_LvNL$estimate > 0,
  "an effect was not recovered (power < 0.7)" = all(res$power >= 0.7),
  "too many false positives" = all(fpr <= 0.05),
  "biobank-only shift leaked into the in-study comparison" = sum(biobank$significant) <= 1,
  "BMI association in LEIP not found" = any(leip$significant & leip$parameter == "BMI"),
  "Th2 set not enriched in lesional skin" = any(th2$padj < 0.05 & th2$NES > 0)
)
msg("All pipeline checks passed.")
