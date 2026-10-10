# QC overview only: proteins above / below LOD in serum and dISF (ring charts) and their distribution across
# protein classes. Runs just the three steps this figure needs, not the rest of the pipeline (run_all.R).
#   In RStudio: source("run_qc_overview.R")      Terminal: Rscript run_qc_overview.R
# Results: output/qc_overview/ (qc_overview.pdf / .png, qc_overview.xlsx). Protein classes need the Human Protein
# Atlas table once: source("tools/download_hpa.R").
if (!file.exists("R/utils.R"))
  stop("Working directory must be the project folder (the one containing R/ and scripts/). Current: ", getwd())
source("R/utils.R")
cfg_run <- load_config()
check_inputs(cfg_run)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in c("scripts/01_metadata.R", "scripts/02_import_qc.R", "scripts/02b_qc_overview.R")) {
  message("\n==== ", s, " ====")
  if (system2(rscript, s) != 0) stop(s, " failed")
}
message("\nQC overview finished: ", file.path(cfg_run$paths$output, "qc_overview"), "/")
