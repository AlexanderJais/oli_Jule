# Test of the visit 3 analysis (run_visit3.R) on simulated data: the simulated AD / lesional effects are found,
# false positives stay rare, and all files are written.
#   Rscript tests/test_visit3.R
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
old_cfg <- Sys.getenv("OLINK_CONFIG", NA)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
run_ok <- system2(rscript, "run_visit3.R") == 0
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot("visit 3 analysis failed on simulated data" = run_ok)

out <- "output_sim/visit3"
truth <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
rd <- \(f, sh) readxl::read_excel(file.path(out, f), sh) |> left_join(truth |> select(OlinkID, role), by = "OlinkID")
hl <- rd("V3_DEGs_dISF.xlsx", "healthy_vs_lesional"); lnl <- rd("V3_DEGs_dISF.xlsx", "lesional_vs_nonlesional")
ser <- rd("V3_DEGs_serum.xlsx", "healthy_vs_AD")
# false positives = simulated null proteins with FDR < 0.05; Benjamini-Hochberg admits more of them the more true effects there
# are, so allow 2 or 10 % of the significant proteins (twice the nominal 5 %), whichever is larger
fp_ok <- \(r) sum(r$significant & r$role == "null", na.rm = TRUE) <= max(2, 0.1 * sum(r$significant, na.rm = TRUE))
stopifnot(
  "dISF lesion vs healthy: AD effects not found" = sum(hl$significant & hl$role %in% c("ISF_AD_vs_HC", "ISF_lesional")) >= 10,
  "dISF lesion vs non-lesional: lesional effects not found" = sum(lnl$significant & lnl$role == "ISF_lesional") >= 10,
  "serum AD vs healthy: AD effects not found" = sum(ser$significant & ser$role == "Serum_AD_vs_HC") >= 3,
  "too many false positives" = fp_ok(hl) && fp_ok(lnl) && fp_ok(ser),
  "serum has no lesional vs non-lesional sheet" = !"lesional_vs_nonlesional" %in% readxl::excel_sheets(file.path(out, "V3_DEGs_serum.xlsx")),
  "files missing" = all(file.exists(file.path(out, c("V3_GSEA.xlsx", "V3_dISF_vs_serum.xlsx", "V3_dISF_vs_serum.pdf",
                                                     "volcano/V3_volcano_dISF_healthy_vs_lesional.pdf", "volcano/V3_volcano_dISF_lesional_vs_nonlesional.pdf",
                                                     "volcano/V3_volcano_serum_healthy_vs_AD.pdf", "gsea/V3_GSEA_dISF_healthy_vs_lesional.pdf")))),
  "dISF vs serum workbook incomplete" = all(c("key_numbers", "proteins", "gene_sets") %in% readxl::excel_sheets(file.path(out, "V3_dISF_vs_serum.xlsx")))
)
msg("All visit 3 checks passed.")
