# End-to-end test of the LEIP galanin study with invented data (no study data needed):
# make the test data -> run all steps -> check that the built-in effects are found (see make_test_data.R).
#   In RStudio (project folder): source("tests/test_leip.R")   - takes about 2-4 minutes
# Must end with "All LEIP galanin checks passed".
if (!file.exists("run_all.R") || !file.exists("R/leip.R"))
  stop("Working directory must be the project folder (the one containing run_all.R and R/). ",
       "In RStudio open oli_Jule.Rproj, or run setwd(\"path/to/oli_Jule\") first. Current: ", getwd())
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_test_data.R") == 0)
unlink("tests/output", recursive = TRUE)
old <- Sys.getenv("LEIP_CONFIG")
Sys.setenv(LEIP_CONFIG = "tests/data/config_test.yml")
status <- system2(rscript, "run_all.R")
if (nzchar(old)) Sys.setenv(LEIP_CONFIG = old) else Sys.unsetenv("LEIP_CONFIG")
stopifnot("the LEIP study did not run through" = status == 0)

out   <- "tests/output"
rd    <- \(...) read_csv(file.path(out, ...), show_col_types = FALSE)
truth <- read_csv("tests/data/truth.csv", show_col_types = FALSE)
smp <- rd("0_data", "samples.csv"); chk <- rd("0_data", "sample_checks.csv")
det <- rd("0_data", "detection.csv") |> left_join(truth, by = c("OlinkID", "Assay"))
agr <- rd("1_elisa_validation", "agreement.csv") |> filter(analysis == "all samples")
ev  <- rd("1_elisa_validation", "ELISA_vs_all_proteins.csv"); bm <- rd("1_elisa_validation", "lab_vs_Olink_benchmark.csv")
ne  <- rd("2_olink_galanin", "neuroendocrine_proteins.csv"); gcl <- rd("2_olink_galanin", "olink_galanin_vs_clinical.csv")
pcs <- rd("2_olink_galanin", "proteome_axes.csv"); gpr <- rd("2_olink_galanin", "olink_galanin_vs_proteins.csv")
lip <- rd("3_galanin_hdl", "galanin_vs_lipids.csv"); hp <- rd("3_galanin_hdl", "lipoprotein_proteins.csv")
dsc <- rd("3_galanin_hdl", "discordance_vs_HDL.csv") |> filter(with == "HDL cholesterol", discordance == "signed", analysis == "all persons")
ahd <- rd("3_galanin_hdl", "agreement_by_HDL.csv"); agr_all <- rd("1_elisa_validation", "agreement.csv")
idt <- rd("1_elisa_validation", "sample_identity.csv"); swapped <- idt |> filter(SubjectID %in% c("LEIP_05", "LEIP_06"))
scr <- rd("4_clinical_screen", "associations_all.csv.gz"); san <- rd("4_clinical_screen", "sanity_checks.csv")
ans <- rd("answers.csv")
lp  <- \(m, l, a) lip |> filter(measure == m, lipid == l, analysis == a)
partners <- c("CHGA", "NPY", "SCG2")
# the interaction helper on simulated data: agreement that fades with h, and none at all
source("R/leip.R"); set.seed(1)
sim <- tibble(h = rnorm(300), x = rnorm(300), y = x * (0.7 - 0.5 * h) + rnorm(300, 0, 0.5))
sim_int <- agreement_interaction(sim$x, sim$y, sim$h); sim_null <- agreement_interaction(sim$x, rnorm(300), sim$h)
unlinked <- c("WHR", "c_fett", "Gluc0_mg_dl", "HOMA_IR", "c_CRP", "MDRD_kurz", "C_TRIGLY", "IL10", "RESTRAINT", "c_tsh")
stopifnot(
  "01: the failed sample was not left out, or the checks miss it" =
    !"S216" %in% smp$SampleID && grepl("failed", chk$problem[chk$SampleID == "S216"]) && grepl("warning", chk$problem[chk$SampleID == "S190"]),
  "01: proteins below LOD not recognised" = !any(det$measurable[det$role == "below_LOD"]) && all(det$measurable[det$role != "below_LOD"]),
  "01: clinical columns with micro sign / umlaut not read" = all(c("Ins0_uU_ml", "lipamisch", "sex_male", "galanin_elisa") %in% names(smp)),
  "02: Olink GAL does not follow the ELISA" = agr$rho > 0.4 && agr$p < 0.01,
  "02: Olink GAL is not among the top proteins the ELISA follows" = ev$rank_by_rho[ev$Assay == "GAL"] <= 10,
  "02: lab vs Olink benchmark wrong" = bm$rho[bm$olink_assay == "FABP4"] > 0.5 && bm$status[bm$olink_assay == "CRP"] == "not measured by Olink",
  "02: the swapped Olink samples of LEIP_05 and LEIP_06 not found, or too many false alarms" =
    nrow(swapped) == 2 && all(grepl("possible swap", swapped$fit)) && sum(grepl("possible swap", idt$fit)) <= 3,
  "02: agreement within sex or plate table missing" = all(c("women only", "men only", "adjusted for sex and Olink plate") %in% agr_all$analysis) &&
    file.exists(file.path(out, "1_elisa_validation", "plate_effects.csv")),
  "2a: the ELISA must not be used in the Olink-only analysis" = !"galanin_elisa" %in% gcl$parameter,
  "2a: built-in link of Olink GAL with age not found" = gcl$rho[gcl$parameter == "age"] < -0.3 && gcl$p[gcl$parameter == "age"] < 0.05,
  "2a: too many false hits among clinical variables without a link to Olink GAL" = sum(gcl$fdr[gcl$parameter %in% unlinked] < 0.05) <= 1,
  "2b: proteins released together with galanin not found (adjusted for age, sex, plate: GAL also falls with age)" =
    all(ne$p_adj[ne$gene %in% partners] < 0.05),
  "2b: false hits among the other neuroendocrine proteins" = sum(ne$p_adj[!ne$gene %in% partners] < 0.05) <= 1,
  "2b: protein table or proteome axes missing" = nrow(pcs) >= 3 && !any(gpr$Assay == "GAL") && nrow(gpr) > 200,
  "04: galanin ELISA - HDL link (built in within sex) not found" =
    lp("galanin ELISA", "C_HDL", "all persons")$rho > 0.4 && lp("galanin ELISA", "C_HDL", "adjusted for sex")$p < 0.05,
  "04: galanin ELISA wrongly linked to LDL" = abs(lp("galanin ELISA", "C_LDL", "adjusted for sex")$rho) < 0.4,
  "04: the Olink HDL protein score does not track the lab HDL" = hp$`rho with lab HDL-C`[hp$protein == "HDL protein score"][1] > 0.8,
  "04: ELISA-Olink discordance does not rise with HDL" = dsc$rho > 0.2,
  "04: agreement-by-HDL table incomplete" = all(c("HDL cholesterol", "ApoA-I (lab)") %in% ahd$with) && all(is.finite(ahd$interaction)),
  "04: interaction helper misses a built-in interaction, or finds one where there is none" = sim_int$p < 1e-4 && sim_int$interaction < 0 && sim_null$p > 0.01,
  "05: sanity checks wrong" = all(san$status[san$protein %in% c("LEP", "APOA1")] == "recovered") &&
    all(san$unadjusted[san$protein %in% c("LEP", "APOA1")] == "recovered") && all(!is.na(san$rho_adj[san$protein %in% c("LEP", "APOA1")])) &&
    san$status[san$protein == "NOT_ON_PANEL"] == "protein not measured",
  "05: too many false hits for clinical variables without a built-in link" = sum(scr$significant[scr$parameter %in% unlinked]) <= 2,
  "06: report or answers incomplete" = file.size(file.path(out, "LEIP_galanin_report.pdf")) > 50000 &&
    all(c("1", "2", "3", "S") %in% substr(ans$question, 1, 1)) &&
    any(grepl("agreement weaken", ans$item)) && any(grepl("correctly matched", ans$item)) &&
    any(grepl("^2a", ans$item)) && any(grepl("^2b", ans$item)) && any(grepl("follow instead", ans$item))
)
msg("All LEIP galanin checks passed.")
