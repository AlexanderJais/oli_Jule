# V1 relapse screen: which dISF proteins are already higher at baseline (V1) in the AD patients who later relapse?
# Ranked "actionable targets" with an exact relabelling check, leave-one-out, AUC, persistence at V2, serum, drug targets.
# Runs only the import / QC (steps 01, 02) and this analysis, not the rest of the pipeline (run_all.R).
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
}
message("\nV1 relapse screen finished: ", file.path(cfg_run$paths$output, "relapse_V1"), "/")
