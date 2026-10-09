# Test of the Leipzig analysis (run_leipzig.R) on simulated data: LEIP_35 has no clinical data in the files, but the
# simulation makes it male, 72 years, BMI 31, HDL 1.1, with 4 proteins raised and 3 lowered only in LEIP_35.
#   Rscript tests/test_leipzig.R
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
old_cfg <- Sys.getenv("OLINK_CONFIG", NA)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
run_ok <- system2(rscript, "run_leipzig.R") == 0
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot("Leipzig analysis failed on simulated data" = run_ok)

out <- "output_sim/leipzig"
prof <- readxl::read_excel(file.path(out, "LEIP_35_profile.xlsx"), "profile")
mk <- file.path(out, "LEIP_35_markers.xlsx")
hits <- tibble(Assay = c(readxl::read_excel(mk, "elevated")$Assay, readxl::read_excel(mk, "decreased")$Assay))
truth <- read.csv("data_sim/leip_case_truth.csv")
tr <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
up <- tr$Assay[tr$role == "LEIP_case_shift"]; down <- tr$Assay[tr$role == "LEIP_case_down"]; shift <- c(up, down)
age <- prof |> filter(parameter == "age")
npx <- readRDS(file.path(out, "data/leipzig_npx.rds"))
stopifnot(
  "only Leipzig samples read" = all(npx$cohort == "LEIP") && n_distinct(npx$SampleID) == 35,
  "sex not recovered" = prof$verdict[prof$parameter == "sex"] == "male",
  "age not recovered" = age$verdict == "estimate" && age$lo80 <= truth$age && age$hi80 >= truth$age,
  "glucose should be not predictable" = prof$verdict[prof$parameter == "Gluc0_mg_dl"] == "not predictable",
  "changed proteins not found, or false hits" = sum(hits$Assay %in% shift) >= 6 && sum(!hits$Assay %in% shift) <= 1,
  "raised proteins not in the elevated sheet" = sum(readxl::read_excel(mk, "elevated")$Assay %in% up) >= 3,
  "lowered proteins not in the decreased sheet" = sum(readxl::read_excel(mk, "decreased")$Assay %in% down) >= 2,
  "PDFs missing" = all(file.exists(file.path(out, c("LEIP_35_profile.pdf", "LEIP_35_markers.pdf", "LEIP_35_markers.png"))))
)
msg("All Leipzig checks passed.")
