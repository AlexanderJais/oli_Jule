# Relapse and lesional signatures per visit (V1-V4): relapse vs non-relapse at the lesion site and lesion site vs non-lesional
# skin at each visit (limma, volcano, GSEA), AKR1B10 and the mTORC1 gene set over the visits; figures as PNG for PowerPoint.
# Part of run_all.R - this is only a shortcut to rerun this one analysis (import / QC, steps 01-02, plus the analysis).
#   In RStudio: source("run_visits_v1_v4.R")      Terminal: Rscript run_visits_v1_v4.R
# Settings: config.yml -> visits_v1_v4. Results: output/visits_V1_V4/
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("scripts/01_metadata.R", "scripts/02_import_qc.R", "visits_v1_v4/visits_v1_v4_analysis.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
  log_step(cfg_run, s)
}
message("\nPer-visit analysis finished: ", file.path(cfg_run$paths$output, "visits_V1_V4"), "/  (start with Visits_V1-V4_summary.pdf)")
