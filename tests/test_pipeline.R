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
ov15 <- read_csv(file.path(out, "serum_vs_disf/overlap_summary.csv"), show_col_types = FALSE)
fov  <- read_csv(file.path(out, "focus/focus_overview.csv"), show_col_types = FALSE)
ex_rl <- read_csv(file.path(out, "export/RELAD2/RELAD2_Serum_NPX_wide.csv"), show_col_types = FALSE)
kqa <- read_csv(file.path(out, "key_questions/answers.csv"), show_col_types = FALSE)
lb    <- \(f) read_csv(file.path(out, "leip_biobank", f), show_col_types = FALSE)
lb_as <- lb("associations_all.csv.gz") |> left_join(truth |> select(OlinkID, role, bmi_slope), by = "OlinkID")
lb_ag <- lb("galanin/olink_vs_elisa.csv") |> filter(analysis == "all samples")
lb_gp <- lb("galanin/GAL_vs_proteins.csv") |> left_join(truth |> select(OlinkID, role), by = "OlinkID")
lb_ev <- lb("galanin/ELISA_vs_proteins.csv"); lb_bm <- lb("galanin/lab_vs_olink_benchmark.csv")
lb_pa <- lb("parameters.csv"); lb_sa <- lb("sanity_checks.csv"); lb_an <- lb("answers.csv")
why   <- \(p) coalesce(lb_pa$reason[lb_pa$parameter == p], "analysed")
stopifnot(
  "LEIP step 18: Olink GAL does not follow the simulated galanin ELISA" =
    lb_ag$rho > 0.4 && lb_ag$p < 0.01 && lb_ev$rank_by_rho[lb_ev$Assay == "GAL"] <= 3,
  "LEIP step 18: GAL partner proteins not among the top GAL correlates" = sum(head(lb_gp$role, 10) == "LEIP_GAL_partner") >= 2,
  "LEIP step 18: lab vs Olink benchmark wrong" =
    lb_bm$rho[lb_bm$olink_assay == "FABP4"] > 0.4 && lb_bm$status[lb_bm$olink_assay == "CRP"] == "not measured by Olink",
  "LEIP step 18: BMI-linked proteins not found" = mean(lb_as$significant[lb_as$parameter == "BMI" & lb_as$bmi_slope > 0]) >= 0.8,
  "LEIP step 18: false hits for unlinked parameters" = sum(lb_as$significant[lb_as$parameter %in% c("RESTRAINT", "c_tsh", "IL10")]) <= 3,
  "LEIP step 18: clinical parameter selection wrong" =
    str_detect(why("ln_BMI"), "log copy") && str_detect(why("gluk_0"), "same ranks as Gluc0_mg_dl") && why("t2d") == "no variation" &&
    str_detect(why("RE_BIN"), "config") && why("galanin_elisa") == "analysed" && why("sex_male") == "analysed",
  "LEIP step 18: sanity checks wrong" = all(lb_sa$status[lb_sa$protein != "NOT_ON_PANEL"] == "recovered") &&
    lb_sa$status[lb_sa$protein == "NOT_ON_PANEL"] == "protein not measured",
  "LEIP step 18: answers or summary PDF missing" = all(c("G1", "G2", "G3", "G4", "P1", "P2") %in% str_extract(lb_an$question, "^[A-Z][0-9]")) &&
    file.size(file.path(out, "leip_biobank", "LEIP_biobank_summary.pdf")) > 20000
)
stopifnot(
  "key questions incomplete (step 15)" = all(paste0("Q", 1:5) %in% str_extract(kqa$question, "^Q[0-9]")) &&
    any(str_detect(kqa$verdict[str_detect(kqa$question, "^Q1") & kqa$item == "Mast cell score"], "elevated in AD")),
  "export missing or incomplete (step 17)" = nrow(ex_rl) == 76 && all(ex_rl$cohort == "RELAD2") &&
    file.exists(file.path(out, "export/ISF_NPX_wide.csv")) && file.exists(file.path(out, "export/proteins.csv")),
  "serum vs dISF overlap missing (step 14)" = nrow(ov15) > 0 && any(ov15$visit == "all visits") &&
    file.exists(file.path(out, "serum_vs_disf/venn_AD_vs_healthy_nominal.png")),
  "focus overview missing a simulated focus protein (step 12)" =
    all(c("TNFRSF9", "TNFSF9", "KITLG", "CPA4", "FCER1A", "TPSAB1", "TPSD1", "POSTN") %in% fov$protein),
  "no protein regulated at all visits found (step 13)" = sum(vc$all_visits_nominal & vc$role != "null") >= 3,
  "executive summary PDF missing (step 16)" = file.exists(file.path(out, "Executive_summary.pdf")) &&
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
