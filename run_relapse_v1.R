# V1 relapse screen: which dISF proteins are already higher at baseline (V1) in the AD patients who later relapse?
# Ranked "actionable targets" with an exact relabelling check, leave-one-out, AUC, persistence at V2, serum, drug targets.
# Runs only the import / QC (steps 01, 02) and this analysis, not the rest of the pipeline (run_all.R). At the end it rebuilds
# the main Executive_summary.pdf (step 18) with the predictor plots, if run_all.R has been run before with the same manifest.
#   In RStudio: source("run_relapse_v1.R")      Terminal: Rscript run_relapse_v1.R
# Settings: config.yml -> relapse_v1. Results: output/relapse_V1/
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("scripts/01_metadata.R", "scripts/02_import_qc.R", "relapse_v1/relapse_v1_analysis.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
  log_step(cfg_run, s)
}
# main executive summary (step 18): only if run_all.R has run with this manifest - otherwise most of its pages would be empty,
# or made from another manifest
steps_done <- steps_log(cfg_run); same <- steps_done$manifest_md5 == manifest_md5(cfg_run)
if ("scripts/17_disf_serum_signatures.R" %in% steps_done$step[same] &&
    all(same | steps_done$step %in% c("scripts/18_summary_report.R", "scripts/19_export_data.R"))) {
  message("\n==== scripts/18_summary_report.R (main executive summary, with the predictor plots) ====")
  if (system2(rscript, "scripts/18_summary_report.R") != 0) stop("scripts/18_summary_report.R failed")
  log_step(cfg_run, "scripts/18_summary_report.R")
  message("Executive_summary.pdf updated: ", file.path(cfg_run$paths$output, "Executive_summary.pdf"))
} else message("\nExecutive_summary.pdf not updated: run source(\"run_all.R\") once with this manifest first (the predictor plots are added on every later run of this script).")
message("\nV1 relapse screen finished: ", file.path(cfg_run$paths$output, "relapse_V1"), "/")
