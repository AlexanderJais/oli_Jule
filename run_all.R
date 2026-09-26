# Run the whole pipeline in order.
#   Rscript run_all.R                 # real data (config.yml)
#   OLINK_CONFIG=data_sim/config_sim.yml Rscript run_all.R   # simulated data
# Resume from a later step (earlier results are reused): in R/RStudio
#   start_at <- 8; source("run_all.R")
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). ",
       "In RStudio open oli_Jule.Rproj, or run setwd(\"path/to/oli_Jule\") first. Current: ", getwd())
source("R/utils.R")
check_inputs(load_config())        # stops with a clear list if a required input is missing
rscript <- file.path(R.home("bin"), "Rscript")   # works on Windows/RStudio without Rscript on PATH
steps <- c("scripts/01_metadata.R", "scripts/02_import_qc.R", "scripts/03_explore.R",
           "scripts/04_isf_models.R", "scripts/05_serum_models.R", "scripts/06_isf_vs_serum.R",
           "scripts/07_enrichment.R", "scripts/08_leip_reference.R", "scripts/09_isf_profile.R",
           "scripts/10_matrix_comparison.R", "scripts/11_trajectories.R",
           "scripts/12_focus_proteins.R")
if (!exists("start_at")) start_at <- 1
if (start_at > 1) message("Starting at step ", start_at, " (reusing earlier results)")
for (s in steps[start_at:length(steps)]) {
  message("\n==== ", s, " ====")
  status <- system2(rscript, s)
  if (status != 0) stop(s, " failed (exit ", status, ")")
}
message("\nAll steps finished. Results are in the output folder set in the config.")
