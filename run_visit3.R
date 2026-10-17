# Dedicated analysis of visit 3 (MicroAD): differential proteins in dISF (lesion site vs healthy skin, lesion site vs
# non-lesional skin) and serum (AD vs healthy), volcano plots, GSEA, and what dISF shows that serum does not.
# Runs only the import / QC (steps 01, 02) and this analysis, not the rest of the pipeline (run_all.R).
#   In RStudio: source("run_visit3.R")      Terminal: Rscript run_visit3.R
# Settings: config.yml -> visit3 (visit, number of labelled proteins). Results: output/visit3/
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("scripts/01_metadata.R", "scripts/02_import_qc.R", "visit3/visit3_analysis.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
  log_step(cfg_run, s)
}
message("\nVisit 3 analysis finished: ", file.path(cfg_run$paths$output, "visit3"), "/")
