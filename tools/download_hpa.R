# Download the Human Protein Atlas gene table (public reference data, ~7 MB) for the tissue-origin
# annotation of step 17. Run once from the project folder:  source("tools/download_hpa.R")
# The file goes to data/reference/ (git-ignored); config.yml -> paths$hpa points to it.
dest <- "data/reference/proteinatlas.tsv.zip"
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
options(timeout = max(600, getOption("timeout")))
download.file("https://www.proteinatlas.org/download/proteinatlas.tsv.zip", dest, mode = "wb")
message("Saved ", dest, " (", round(file.size(dest) / 1e6, 1), " MB). Re-run step 17: start_at <- 17; source(\"run_all.R\")")
