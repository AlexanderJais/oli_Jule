# Thresholds: which thresholds remove samples or proteins or decide results, and are they too strict for these data? The detection
# (LOD) filter, FDR, minimum group sizes, sample QC and the relapse criteria, each measured on the data at other values.
# Part of run_all.R - this is only a shortcut to rerun this one analysis (import / QC, steps 01-02, plus the analysis).
#   In RStudio: source("run_thresholds.R")      Terminal: Rscript run_thresholds.R
# Settings: config.yml -> thresholds (only the values compared). Results: output/thresholds/
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("scripts/01_metadata.R", "scripts/02_import_qc.R", "thresholds/thresholds_analysis.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
  log_step(cfg_run, s)
}
message("\nThreshold analysis finished: ", file.path(cfg_run$paths$output, "thresholds"), "/  (start with Thresholds_summary.pdf)")
