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
md  <- read_csv(file.path(out, "metadata/sample_metadata.csv"), show_col_types = FALSE)
relad_x <- file.path(out, "relad/RELAD_RELAD2_serum_results.xlsx")
relad_rl <- rd("relad/RELAD_RELAD2_serum_results.csv") |> filter(model == "relapse RELAD+RELAD2")
tnf_pw <- read_csv(file.path(out, "tnfrsf9/proteome_wide_lesional_site.csv"), show_col_types = FALSE) |>
  left_join(truth |> select(OlinkID, role, tnf_coupled), by = "OlinkID")
tnf_ans <- read_csv(file.path(out, "tnfrsf9/answers.csv"), show_col_types = FALSE)
sig_ans <- read_csv(file.path(out, "signatures/answers.csv"), show_col_types = FALSE)
print(tnf_pw |> filter(tnf_coupled | rank <= 8) |> select(rank, Assay, partial_rho, mixed_fdr, role, tnf_coupled) |> as.data.frame(), digits = 2)
stopifnot(
  "manifest Ver2 codes not harmonised (step 01)" = any(md$cohort == "RELAD2" & md$relapse %in% "active") &&
    all(!is.na(md$sex[md$cohort == "MicroAD"])) && any(md$clinical_state %in% "remission") && !any(md$clinical_state %in% "helthy"),
  "serum AD vs healthy for dISF comparisons must be MicroAD only" =
    all(ser$n_samples[ser$model == "AD_vs_HC_MicroAD"] < sum(md$cohort == "MicroAD" & md$matrix == "Serum" & md$visit %in% "V1")),
  "RELAD / RELAD2 workbook incomplete (step 18)" = file.exists(relad_x) &&
    all(c("summary", "all_results", "NPX_values", "detection", "key_proteins") %in% readxl::excel_sheets(relad_x)) &&
    mean(relad_rl$significant[relad_rl$role == "Serum_relapse"]) >= 0.7,
  "TNFRSF9-coupled proteins not found proteome-wide (step 19)" = sum(tnf_pw$tnf_coupled & coalesce(tnf_pw$mixed_fdr, 1) < 0.05) >= 4,
  "group-driven proteins leak into the adjusted TNFRSF9 ranking (step 19)" = sum(!tnf_pw$tnf_coupled & coalesce(tnf_pw$mixed_fdr, 1) < 0.05) <= 3,
  "TNFRSF9 - IL33 within-group correlation not found (step 19)" =
    any(str_detect(tnf_ans$item, "IL33") & str_detect(tnf_ans$verdict, "within groups")),
  "serum vs dISF signature answers incomplete (step 20)" = all(as.character(1:6) %in% str_extract(sig_ans$question, "^[0-9]")) &&
    file.exists(file.path(out, "signatures/signatures.xlsx")),
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
