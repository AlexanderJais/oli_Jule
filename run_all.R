# Run the whole pipeline in order - the ONE script to run. It includes the dedicated analyses (QC overview = step 02b,
# then visit 3, V1 relapse screen, relapse follow-up, visits V1-V4 right after it: they only need the imported data, so their folders
# appear in the first minutes). The separate run_*.R scripts are only shortcuts to rerun one of them quickly.
#   Rscript run_all.R                 # real data (config.yml)
#   OLINK_CONFIG=data_sim/config_sim.yml Rscript run_all.R   # simulated data
# Resume from a later step (earlier results are reused): in R/RStudio
#   start_at <- 8; source("run_all.R")        # step numbers = script numbers; start_at <- "02b" starts at 02b
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
steps <- c("scripts/01_metadata.R", "scripts/02_import_qc.R", "scripts/02b_qc_overview.R", "scripts/03_explore.R",
           "scripts/04_isf_models.R", "scripts/05_serum_models.R", "scripts/06_isf_vs_serum.R",
           "scripts/07_enrichment.R", "scripts/08_leip_reference.R", "scripts/09_isf_profile.R",
           "scripts/10_matrix_comparison.R", "scripts/11_trajectories.R",
           "scripts/12_focus_proteins.R", "scripts/13_visit_course.R", "scripts/14_serum_vs_disf.R", "scripts/15_key_questions.R",
           "scripts/16_tnfrsf9_correlation.R", "scripts/17_disf_serum_signatures.R",
           "scripts/18_summary_report.R", "scripts/19_export_data.R")
if (!exists("start_at")) start_at <- 1
run_steps <- select_steps(steps, start_at)      # by script number: start_at <- 3 runs 03, 04, ...; "02b" starts at 02b
if (length(run_steps) < length(steps)) message("Starting at ", basename(run_steps[1]), " (reusing earlier results)")
# dedicated analyses (not numbered): right after step 02b, as they only need steps 01-02 (also before the summary, step 18)
dedicated <- c("visit3/visit3_analysis.R", "relapse_v1/relapse_v1_analysis.R", "relapse_followup/relapse_followup_analysis.R",
               "visits_v1_v4/visits_v1_v4_analysis.R")
i02b <- match("scripts/02b_qc_overview.R", run_steps)
if (!is.na(i02b)) run_steps <- append(run_steps, dedicated, after = i02b)
# results reused when resuming must come from the manifest config.yml names now (log: output/metadata/steps_run.csv)
steps_done <- steps_log(cfg_run)
stale <- setdiff(setdiff(c(steps, dedicated), run_steps), steps_done$step[steps_done$manifest_md5 == manifest_md5(cfg_run)])
if (length(stale))
  stop("Cannot resume: these results were made from another sample manifest (or before the manifest was logged): ",
       paste(basename(stale), collapse = ", "), ". Run everything from step 01: rm(start_at); source(\"run_all.R\")", call. = FALSE)
for (s in run_steps) {
  message("\n==== ", s, " ====")
  status <- system2(rscript, s)
  if (status != 0) stop(s, " failed (exit ", status, ")")
  log_step(cfg_run, s)
}
# forget start_at after a successful run, so the next source("run_all.R") is a full run again
if (exists("start_at", envir = globalenv(), inherits = FALSE)) rm("start_at", envir = globalenv())
message("\nAll steps finished. Results are in: ", cfg_run$paths$output, "/  (start with Executive_summary.pdf)")
