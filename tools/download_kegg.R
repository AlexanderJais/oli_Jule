# Download the current KEGG human pathways (pathway -> genes, with gene symbols) from the KEGG REST server
# (rest.kegg.jp) for the KEGG analysis of the QC overview (step 02b). Run once from the project folder:
#   source("tools/download_kegg.R")
# Only KEGG's public lists are downloaded; nothing is sent. KEGG is free for academic use
# (https://www.kegg.jp/kegg/legal.html); the file goes to data/reference/ (git-ignored) and is not shared.
# config.yml -> paths$kegg points to it. Re-run it to update KEGG; the release date is stored in the file.
dest <- "data/reference/kegg_hsa_pathways.tsv"
dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
options(timeout = max(600, getOption("timeout")))
get <- function(what) {
  f <- tempfile(); on.exit(unlink(f))
  download.file(paste0("https://rest.kegg.jp/", what), f, quiet = TRUE, mode = "wb")
  readLines(f, warn = FALSE, encoding = "UTF-8")
}
split_tab <- function(x, n) { m <- do.call(rbind, lapply(strsplit(x, "\t", fixed = TRUE), `length<-`, n)); m[is.na(m)] <- ""; m }
message("Downloading KEGG human pathways from rest.kegg.jp ...")
pw  <- split_tab(get("list/pathway/hsa"), 2)          # hsa04060 | Cytokine-cytokine receptor interaction - Homo sapiens (human)
lnk <- split_tab(get("link/pathway/hsa"), 2)          # hsa:3558 | path:hsa04060
gn  <- split_tab(get("list/hsa"), 4)                  # hsa:3558 | CDS | position | IL2, IL-2, TCGF; interleukin 2
info <- get("info/kegg")
release <- sub("^\\s*pathway\\s+[0-9,]+\\s+", "", grep("^\\s*pathway\\s", info, value = TRUE)[1])
if (is.na(release) || !nzchar(release)) release <- paste("downloaded", Sys.Date())
symbol <- ifelse(grepl(";", gn[, 4]), trimws(sub(",.*$", "", sub(";.*$", "", gn[, 4]))), NA_character_)
out <- data.frame(pathway_id = sub("^path:", "", lnk[, 2]), gene_id = lnk[, 1])
out$pathway_name <- sub(" - Homo sapiens \\(human\\)$", "", pw[match(out$pathway_id, pw[, 1]), 2])
out$gene_symbol <- symbol[match(out$gene_id, gn[, 1])]
out$kegg_release <- paste("KEGG pathway database", release)
out <- out[!is.na(out$gene_symbol) & !is.na(out$pathway_name), c("pathway_id", "pathway_name", "gene_id", "gene_symbol", "kegg_release")]
write.table(out, dest, sep = "\t", row.names = FALSE, quote = FALSE)
message(sprintf("Saved %s: %d pathways, %d genes (%s). Re-run the QC overview: source(\"run_qc_overview.R\")",
                dest, length(unique(out$pathway_id)), length(unique(out$gene_symbol)), out$kegg_release[1]))
