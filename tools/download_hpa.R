# Download the Human Protein Atlas gene table (public reference data, ~7 MB) for the protein classes of the
# QC overview (step 02b) and the tissue-origin annotation of step 17. Run once from the project folder:  source("tools/download_hpa.R")
# The file goes to data/reference/ (git-ignored); config.yml -> paths$hpa points to it.
dest <- "data/reference/proteinatlas.tsv.zip"
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
options(timeout = max(600, getOption("timeout")))
download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", dest, mode = "wb")
message("Saved ", dest, " (", round(file.size(dest) / 1e6, 1), " MB). Re-run from step 02b: start_at <- \"02b\"; source(\"run_all.R\")")
