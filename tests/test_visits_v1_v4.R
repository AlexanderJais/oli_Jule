# Test of the per-visit analysis (run_visits_v1_v4.R) on simulated data, with the manifest saved as a single table like the
# 2026-10-15 manifest (one sheet "Sheet1", extra column ClinicalPhase, a stray column without header): the 5 proteins simulated
# higher in relapsers must be found at V1 and V2 and in the signature over the visits; AKR1B10 and the mTORC1 gene set must be
# higher at the lesion site at every visit; V4 (2 relapsers before their relapse) is not tested; figures are PNG and PDF;
# results made from another manifest are never reused.
#   Rscript tests/test_visits_v1_v4.R
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
m <- readxl::read_excel("data_sim/manifest.xlsx", sheet = "manifest", col_types = "text") |>
  mutate(ClinicalPhase = if_else(Relapse %in% "relapse" & ClinicalStateSkin %in% "lesional" & Visit != "V1", "Relapse", "other"),
         stray = if_else(row_number() == 3, "copy of another row", NA_character_))
names(m)[names(m) == "stray"] <- ""
writexl::write_xlsx(list(Sheet1 = m), "data_sim/manifest_single_sheet.xlsx")
cfg_sim <- yaml::read_yaml("data_sim/config_sim.yml"); cfg_sim$paths$manifest <- "data_sim/manifest_single_sheet.xlsx"
yaml::write_yaml(cfg_sim, "data_sim/config_sim_single_sheet.yml")
old_cfg <- Sys.getenv("OLINK_CONFIG", NA)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim_single_sheet.yml")
run_ok <- system2(rscript, "run_visits_v1_v4.R") == 0
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot("per-visit analysis failed on simulated data" = run_ok)

out <- "output_sim/visits_V1_V4"
rd <- \(f, sh) readxl::read_excel(file.path(out, f), sh)
truth <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
base <- truth$Assay[truth$role == "ISF_relapse_baseline"]
up_at <- \(sh) { r <- rd("1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx", sh); sum(r$Assay %in% base & r$p < 0.05 & r$logFC > 0) }
sg <- rd("1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx", "signature_over_visits")
akr <- rd("3_lesional_vs_nonlesional_DEG_V1-V4.xlsx", "AKR1B10")
mt <- rd("5_lesional_vs_nonlesional_GSEA_V1-V4.xlsx", "signature")
flags <- read_csv("output_sim/metadata/data_flags.csv", show_col_types = FALSE)
stopifnot(
  "single-table manifest: ClinicalPhase check flags consistent relapse visits" = !any(grepl("ClinicalPhase", flags$issue)),
  "simulated relapse markers not found at V1 / V2" = up_at("V1") >= 4 && up_at("V2") >= 3,
  "V4 should not be tested (2 relapsers before their relapse)" = grepl("not tested", rd("1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx", "V4")$note[1]),
  "relapse signature over the visits not found" = sum(sg$Assay %in% base & sg$direction %in% "higher in relapsers at >= 2 visits") >= 4,
  "AKR1B10 not higher at the lesion site at every visit" = nrow(akr) == 4 && all(akr$log2FC > 0 & akr$p < 0.05),
  "mTORC1 gene set not higher at the lesion site at every visit" = nrow(mt) == 4 && all(mt$NES > 0 & mt$p < 0.05),
  "files missing" = all(file.exists(file.path(out, c("1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx", "2_relapse_vs_nonrelapse_GSEA_V1-V4.xlsx",
    "3_lesional_vs_nonlesional_DEG_V1-V4.xlsx", "5_lesional_vs_nonlesional_GSEA_V1-V4.xlsx", "1_relapse_volcano_V1.png", "1_relapse_volcano_V1.pdf",
    "2_relapse_GSEA_V1.png", "3_lesional_volcano_V1.png", "4_AKR1B10_lesional_vs_nonlesional_over_visits.png", "5_mTORC1_lesional_vs_nonlesional_by_visit.png",
    "Visits_V1-V4_summary.pdf")))),
  "no PowerPoint file is written (figures are PNG / PDF)" = !any(grepl("\\.pptx$", list.files(out)))
)

# one manifest for all results: a resume must not reuse results made from another manifest (here: the other steps ran with
# data_sim/manifest.xlsx or not at all), and later steps stop if the manifest changed after step 01
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim_single_sheet.yml")
resume <- suppressWarnings(system2(rscript, c("-e", shQuote('start_at <- 12; source("run_all.R")')), stdout = TRUE, stderr = TRUE))
m_changed <- m; m_changed$Note <- "changed"
writexl::write_xlsx(list(Sheet1 = m_changed), "data_sim/manifest_single_sheet.xlsx")
changed <- suppressWarnings(system2(rscript, "visits_v1_v4/visits_v1_v4_analysis.R", stdout = TRUE, stderr = TRUE))
writexl::write_xlsx(list(Sheet1 = m), "data_sim/manifest_single_sheet.xlsx")
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot(
  "steps of the shortcut not logged with their manifest" = all(c("scripts/01_metadata.R", "scripts/02_import_qc.R", "visits_v1_v4/visits_v1_v4_analysis.R") %in%
                                                               read.csv("output_sim/metadata/steps_run.csv")$step),
  "resume with results of another manifest not stopped" = any(grepl("Cannot resume", resume)),
  "changed manifest not detected by a later step" = any(grepl("were made from the manifest", changed))
)
msg("All per-visit checks passed.")
