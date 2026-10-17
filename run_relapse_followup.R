# Relapse follow-up: eight questions on relapse vs non-relapse (MicroAD) - is the comparison right (labels, direction,
# the lower immune signalling in relapsers), how it is calculated, leukocytes (CD45), T cells (LCK, CD3, activation),
# IL4 / CCL27 / CCL5 / CCL8 over the visits, volcano plots with every p < 0.05 labelled (dISF, serum), druggable targets.
# Runs only the import / QC (steps 01, 02) and this analysis, not the rest of the pipeline (run_all.R).
#   In RStudio: source("run_relapse_followup.R")      Terminal: Rscript run_relapse_followup.R
# Settings: config.yml -> relapse_followup. Results: output/relapse_followup/
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("scripts/01_metadata.R", "scripts/02_import_qc.R", "relapse_followup/relapse_followup_analysis.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
  log_step(cfg_run, s)
}
message("\nRelapse follow-up finished: ", file.path(cfg_run$paths$output, "relapse_followup"), "/  (start with Relapse_followup_executive_summary.pdf)")
