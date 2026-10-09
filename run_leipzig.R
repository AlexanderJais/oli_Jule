# Dedicated analysis of the Leipzig biobank samples (separate from the MicroAD pipeline run_all.R).
# One Leipzig sample (default: LEIP_35, the one without clinical data) against the other Leipzig samples:
# which proteins differ, and what the proteome says about the person. Only Leipzig serum samples are read.
#   In RStudio: source("run_leipzig.R")      Terminal: Rscript run_leipzig.R
# Settings: config.yml -> paths (manifest, NPX files, LEIP clinical file) and leip_case. Results: output/leipzig/
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and leipzig/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("leipzig/01_import.R", "leipzig/02_case.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
}
message("\nLeipzig analysis finished. Results are in: ", file.path(cfg_run$paths$output, "leipzig"), "/")
