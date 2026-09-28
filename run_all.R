# LEIP galanin study - runs all steps in order.
#   In RStudio: open oli_Jule.Rproj, then            source("run_all.R")
#   Resume at a step (earlier results are reused):   start_at <- 3; source("run_all.R")
#   Check that everything works, with invented data: source("tests/test_leip.R")
if (!file.exists("R/leip.R") || !file.exists("scripts/01_data.R"))
  stop("Working directory must be the project folder (the one containing run_all.R, R/ and scripts/). ",
       "In RStudio open oli_Jule.Rproj, or run setwd(\"path/to/oli_Jule\") first. Current: ", getwd())
source("R/utils.R")
cfg <- load_config()

p <- cfg$paths
npx <- if (dir.exists(p$npx_dir)) list.files(p$npx_dir, "\\.parquet$") else character()
inputs <- tibble(input = c("Olink NPX parquet file(s)", "LEIP clinical file", "Olink fixed LOD file", "manifest"),
                 required = c("yes", "yes", "recommended", "no"),
                 path = c(file.path(p$npx_dir, "*.parquet"), p$leip_clinical %||% "", p$fixed_lod %||% "", p$manifest %||% ""),
                 found = c(length(npx) > 0, file.exists(p$leip_clinical %||% ""), file.exists(p$fixed_lod %||% ""),
                           file.exists(p$manifest %||% "")))
message("LEIP galanin study - input files (working directory: ", getwd(), "):")
for (i in seq_len(nrow(inputs)))
  message(sprintf("  [%s] %-26s %s", if (inputs$found[i]) "ok" else if (inputs$required[i] == "yes") "MISSING" else "--",
                  inputs$input[i], inputs$path[i]))
if (any(!inputs$found & inputs$required == "yes"))
  stop("Missing input file(s) - see data/README.md.", call. = FALSE)

rscript <- file.path(R.home("bin"), "Rscript")   # works on Windows/RStudio without Rscript on PATH
steps <- file.path("scripts", c("01_data.R", "02_elisa_validation.R", "03_olink_galanin.R",
                                "04_galanin_hdl.R", "05_clinical_screen.R", "06_report.R", "07_platelet_figures.R"))
if (!exists("start_at")) start_at <- 1
if (start_at > 1) message("Starting at step ", start_at, " (reusing earlier results)")
for (s in steps[start_at:length(steps)]) {
  message("\n==== ", s, " ====")
  status <- system2(rscript, s)
  if (status != 0) stop(s, " failed (exit ", status, ")")
}
message("\nDone. Start with ", file.path(p$output, "LEIP_galanin_report.pdf"), "; paper figures: ", file.path(p$output, "paper_platelets"))
