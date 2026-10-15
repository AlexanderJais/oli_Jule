# Test of the relapse follow-up (run_relapse_followup.R) on simulated data: the T-cell markers simulated lower in the dISF of
# relapsers must each come out lower; the 2 baseline markers also raised in all AD skin must be the druggable-target
# candidates meeting all criteria; the label and direction checks pass; the volcano tables label exactly the p < 0.05 proteins.
#   Rscript tests/test_relapse_followup.R
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
old_cfg <- Sys.getenv("OLINK_CONFIG", NA)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
run_ok <- system2(rscript, "run_relapse_followup.R") == 0
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot("relapse follow-up failed on simulated data" = run_ok)

out <- "output_sim/relapse_followup"
rd <- \(f, sh) readxl::read_excel(file.path(out, f), sh)
truth <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
both <- truth |> filter(role == "ISF_relapse_baseline") |> slice(4:5)          # also raised in all AD skin
lab <- rd("1_check_relapse_vs_nonrelapse.xlsx", "labels"); dirn <- rd("1_check_relapse_vs_nonrelapse.xlsx", "direction")
mkr <- rd("3_4_leukocytes_T_cells.xlsx", "markers") |> filter(comparison == "dISF lesion site V1")
cyt <- rd("5_cytokines_over_visits.xlsx", "group_means")
v6 <- rd("6_volcano_dISF_lesion_V1_data.xlsx", "volcano_data"); v7 <- rd("7_volcano_serum_V1_data.xlsx", "volcano_data")
tg <- rd("8_druggable_targets.xlsx", "targets")
t_low <- truth$Assay[truth$role == "ISF_T_cell_relapse_low"]
stopifnot(
  "label check flags consistent simulated labels" = all(lab$check[lab$in_comparison] == "ok"),
  "direction check failed" = all(dirn$same_sign),
  "simulated T-cell markers not lower in relapsers" = with(mkr |> filter(Assay %in% t_low), all(log2FC < 0) && mean(p_exact < 0.05) >= 0.75),
  "CD45 missing" = any(mkr$Assay == "PTPRC" & !is.na(mkr$log2FC)),
  "plate must not appear in the outputs" = !any(grepl("plate", unlist(map(c("1_check_relapse_vs_nonrelapse.xlsx", "3_4_leukocytes_T_cells.xlsx"),
                                                                         \(f) map(readxl::excel_sheets(file.path(out, f)), \(sh) names(rd(f, sh))))), ignore.case = TRUE)),
  "cytokines missing from the heatmap data" = all(c("IL4", "CCL27", "CCL5", "CCL8") %in% cyt$Assay),
  "volcano labels must be exactly the p < 0.05 proteins" = all(v6$labelled == (v6$p < 0.05)) && all(v7$labelled == (v7$p < 0.05)),
  "simulated AD-and-relapse targets not found" = all(both$OlinkID %in% tg$OlinkID[substr(tg$tier, 1, 1) %in% c("1", "2")]),
  "files missing" = all(file.exists(file.path(out, c("Relapse_followup_executive_summary.pdf", "1_sample_level_check.pdf",
                                                     "3_leukocyte_markers.pdf", "4_T_cell_markers.pdf", "4_T_cell_activation_markers.pdf", "5_cytokines_heatmap_group_means.pdf",
                                                     "5_cytokines_heatmap_per_patient.pdf", "6_volcano_dISF_lesion_V1.pdf", "7_volcano_serum_V1.pdf",
                                                     "8_druggable_targets.pdf"))))
)
msg("All relapse follow-up checks passed.")
