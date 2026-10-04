# Run the whole pipeline in order.
#   Rscript run_all.R                 # real data (config.yml)
#   OLINK_CONFIG=data_sim/config_sim.yml Rscript run_all.R   # simulated data
# Resume from a later step (earlier results are reused): in R/RStudio
#   start_at <- 8; source("run_all.R")
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). ",
       "In RStudio open oli_Jule.Rproj, or run setwd(\"path/to/oli_Jule\") first. Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)              # stops with a clear list if a required input is missing
message("Settings: ", cfg_run$config_file, "  ->  results in: ", cfg_run$paths$output, "/")
if (cfg_run$config_file != "config.yml")
  message("NOTE: not the standard config.yml (environment variable OLINK_CONFIG is set). ",
          "For the real analysis restart R or run Sys.unsetenv(\"OLINK_CONFIG\").")
rscript <- file.path(R.home("bin"), "Rscript")   # works on Windows/RStudio without Rscript on PATH
steps <- c("scripts/01_metadata.R", "scripts/02_import_qc.R", "scripts/03_explore.R",
           "scripts/04_isf_models.R", "scripts/05_serum_models.R", "scripts/06_isf_vs_serum.R",
           "scripts/07_enrichment.R", "scripts/08_leip_reference.R", "scripts/09_isf_profile.R",
           "scripts/10_matrix_comparison.R", "scripts/11_trajectories.R",
           "scripts/12_focus_proteins.R", "scripts/13_visit_course.R", "scripts/14_serum_vs_disf.R", "scripts/15_key_questions.R",
           "scripts/16_tnfrsf9_correlation.R", "scripts/17_disf_serum_signatures.R",
           "scripts/18_summary_report.R", "scripts/19_export_data.R")
if (!exists("start_at")) start_at <- 1
if (!is.numeric(start_at) || length(start_at) != 1 || is.na(start_at) || start_at < 1 || start_at > length(steps))
  stop("start_at must be a step number from 1 to ", length(steps), " (it is ", format(start_at), "). Type rm(start_at) for a full run.")
if (start_at > 1) message("Starting at step ", start_at, " (reusing earlier results)")
for (s in steps[start_at:length(steps)]) {
  message("\n==== ", s, " ====")
  status <- system2(rscript, s)
  if (status != 0) stop(s, " failed (exit ", status, ")")
}
# forget start_at after a successful run, so the next source("run_all.R") is a full run again
if (exists("start_at", envir = globalenv(), inherits = FALSE)) rm("start_at", envir = globalenv())
message("\nAll steps finished. Results are in: ", cfg_run$paths$output, "/")
