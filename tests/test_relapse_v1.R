# Test of the V1 relapse screen (run_relapse_v1.R) on simulated data: 5 proteins are simulated higher in the dISF of
# relapsers at every visit; they must top the target list, persist at V2, and the chance check must run.
#   Rscript tests/test_relapse_v1.R
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
old_cfg <- Sys.getenv("OLINK_CONFIG", NA)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
run_ok <- system2(rscript, "run_relapse_v1.R") == 0
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot("V1 relapse screen failed on simulated data" = run_ok)

out <- "output_sim/relapse_V1"; f <- file.path(out, "V1_relapse_targets.xlsx")
truth <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
real <- truth$Assay[truth$role == "ISF_relapse_baseline"]
up <- readxl::read_excel(f, "targets_up"); ch <- readxl::read_excel(f, "chance_check")
stopifnot(
  "simulated baseline markers not among the top 8 targets" = sum(head(up$Assay, 8) %in% real) >= 4,
  "simulated markers should persist at V2" = all(up$persists_V2_up[up$Assay %in% real]),
  "chance check missing" = any(ch$value == "210"),
  "prediction test missing or wrong" = grepl("better than chance", readxl::read_excel(f, "README")$note[2]) &&
    all(c("loo_correct", "cutoff_NPX", "rho_time_to_relapse") %in% names(up)),
  "files missing" = all(file.exists(file.path(out, c("V1_relapse_volcano.pdf", "V1_relapse_targets_course.pdf", "V1_relapse_chance_check.pdf", "V1_relapse_prediction_test.pdf", "V1_relapse_GSEA.xlsx"))))
)
msg("All V1 relapse checks passed.")
