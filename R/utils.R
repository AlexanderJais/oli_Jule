# Shared helpers: configuration, output paths, logging.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

#' Load config.yml (or the file named in env var OLINK_CONFIG).
load_config <- function(path = Sys.getenv("OLINK_CONFIG", "config.yml")) {
  if (!file.exists(path)) stop("Config file not found: ", path)
  cfg <- yaml::read_yaml(path)
  cfg$config_file <- path
  options(olink.cores = cfg$stats$cores %||% 1)
  cfg$paths <- resolve_data_paths(cfg$paths)
  cfg
}

# If a configured file is missing, look for it in the same folder under its usual name,
# so the original file names from Olink / the lab can be used without renaming.
data_file_patterns <- list(
  manifest      = "(manifest|Sample[ _-]?Submission[ _-]?Sheet).*\\.xlsx$",
  leip_clinical = "LEIP.*clinical.*\\.xlsx$",
  fixed_lod     = "Fixed[ _-]?LOD.*\\.csv$",        # e.g. "Explore HT_Fixed LOD.csv", "Explore_HT_Fixed_LOD.csv"
  severity      = "(severity|SCORAD|EASI).*\\.(xlsx|csv)$"
)

resolve_data_paths <- function(paths) {
  for (key in names(data_file_patterns)) {
    p <- paths[[key]]
    if (is.null(p) || file.exists(p)) next
    dir <- dirname(p)
    hits <- if (dir.exists(dir)) list.files(dir, pattern = data_file_patterns[[key]], ignore.case = TRUE, full.names = TRUE) else character()
    hits <- hits[!startsWith(basename(hits), "~$")]          # skip Excel lock files
    if (length(hits) == 1) {
      paths[[key]] <- hits
    } else if (length(hits) > 1) {
      stop(sprintf("Several candidate files for '%s' in %s: %s\nKeep one, or set paths$%s in config.yml.",
                   key, dir, paste(basename(hits), collapse = ", "), key))
    }
  }
  paths
}

#' Check that the input files exist; print what is found and what is missing.
check_inputs <- function(cfg) {
  p <- cfg$paths
  npx <- if (dir.exists(p$npx_dir)) list.files(p$npx_dir, "\\.parquet$") else character()
  status <- tibble(
    input = c("manifest", "NPX parquet files", "LEIP clinical data", "Olink fixed LOD file", "severity scores", "Human Protein Atlas"),
    required = c("yes", "yes", "no", "recommended", "no", "no"),
    path = c(p$manifest, file.path(p$npx_dir, "*.parquet"), p$leip_clinical %||% "", p$fixed_lod %||% "", p$severity %||% "", p$hpa %||% ""),
    found = c(file.exists(p$manifest %||% ""), length(npx) > 0, file.exists(p$leip_clinical %||% ""),
              file.exists(p$fixed_lod %||% ""), file.exists(p$severity %||% ""), file.exists(p$hpa %||% ""))
  )
  message("Input files (working directory: ", getwd(), "):")
  for (i in seq_len(nrow(status)))
    message(sprintf("  [%s] %-22s %s%s", if (status$found[i]) "ok" else if (status$required[i] == "yes") "MISSING" else "--",
                    status$input[i], status$path[i],
                    if (i == 2 && length(npx)) paste0("  (", paste(npx, collapse = ", "), ")") else ""))
  missing <- status$input[!status$found & status$required == "yes"]
  if (length(missing)) {
    dd <- dirname(p$manifest)
    message("\nFiles currently in ", dd, "/: ",
            if (dir.exists(dd)) paste(list.files(dd), collapse = ", ") else "(folder does not exist)")
    stop("Missing: ", paste(missing, collapse = ", "), ". See data/README.md for where each file goes.", call. = FALSE)
  }
  invisible(status)
}

#' Step number of a script from its file name: "02_import_qc.R" -> 2, "02b_qc_overview.R" -> 2.2.
step_number <- function(f) {
  m <- str_match(basename(f), "^(\\d+)([a-z]?)_")
  as.numeric(m[, 2]) + if_else(is.na(m[, 3]) | m[, 3] == "", 0, match(m[, 3], letters) / 10)
}

#' Steps to run when resuming: all scripts from step `start_at` on (a number such as 3, or a name such as "02b").
select_steps <- function(steps, start_at = 1) {
  key <- step_number(steps)
  s <- if (is.character(start_at)) step_number(paste0(start_at, "_")) else if (is.numeric(start_at)) start_at else NA
  if (length(s) != 1 || is.na(s) || !s %in% key)       # must name an existing script
    stop("start_at must be a step number from 1 to ", floor(max(key)), " (or e.g. \"02b\"); it is ",
         paste(format(start_at), collapse = ", "), ". Type rm(start_at) for a full run.", call. = FALSE)
  steps[key >= s]
}

#' Human Protein Atlas gene table (proteinatlas.tsv or the downloaded .zip); NULL if the file is missing.
read_hpa <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  h <- if (str_detect(path, "\\.zip$")) readr::read_tsv(unz(path, "proteinatlas.tsv"), show_col_types = FALSE, guess_max = 1e5)
       else readr::read_tsv(path, show_col_types = FALSE, guess_max = 1e5)
  attr(h, "downloaded") <- format(as.Date(file.mtime(path)))
  h
}

#' Path inside the output folder; creates the sub-folder if needed.
out_path <- function(cfg, ...) {
  p <- file.path(cfg$paths$output, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

#' Parallel back-end for dream/variancePartition (stats: cores in config.yml).
bpparam_cores <- function() {
  n <- getOption("olink.cores", 1)
  if (n > 1 && .Platform$OS.type != "windows") BiocParallel::MulticoreParam(n)
  else if (n > 1) BiocParallel::SnowParam(n)
  else BiocParallel::SerialParam()
}

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...))

#' Remove results of an earlier run so a skipped model cannot leave stale files behind.
clear_outputs <- function(cfg, subdir, pattern = ".*") {
  d <- file.path(cfg$paths$output, subdir)
  if (!dir.exists(d)) return(invisible(0))
  f <- list.files(d, pattern = pattern, full.names = TRUE, recursive = TRUE)
  unlink(f)
  invisible(length(f))
}

#' Write a data frame as CSV and return it invisibly.
save_csv <- function(df, cfg, ...) {
  p <- out_path(cfg, ...)
  readr::write_csv(df, p, na = "")
  invisible(df)
}

save_plot <- function(p, cfg, ..., width = 8, height = 6) {
  f <- out_path(cfg, ...)
  # never exceed ggsave's 50-inch limit, however many proteins/contrasts a plot has
  ggsave(f, p, width = min(width, 45), height = min(height, 45), dpi = 150)
  invisible(f)
}

#' Read an intermediate result written by an earlier script; fail with a clear message.
read_step <- function(cfg, ..., step) {
  p <- file.path(cfg$paths$output, ...)
  if (!file.exists(p)) stop(sprintf("%s not found - run %s first.", p, step))
  readRDS(p)
}

#' Optional severity scores (paths$severity): one row per SubjectID and visit, numeric score columns.
#' Visit may be given as "V1" or as 1; SubjectID as text or number. Returns NULL if no file is set.
read_severity <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  sev <- if (str_detect(path, "\\.xlsx?$")) readxl::read_excel(path) else read_csv(path, show_col_types = FALSE)
  names(sev)[tolower(names(sev)) == "subjectid"] <- "SubjectID"
  names(sev)[tolower(names(sev)) == "visit"] <- "visit"
  if (!all(c("SubjectID", "visit") %in% names(sev))) stop("Severity file needs columns SubjectID and Visit: ", path)
  names(sev) <- make.names(names(sev))            # e.g. "itch NRS" -> "itch.NRS"
  sev |> mutate(SubjectID = as.character(SubjectID),
                visit = if_else(str_detect(as.character(visit), "^\\d+$"), paste0("V", visit), as.character(visit)))
}
severity_scores <- function(sev) setdiff(names(sev)[vapply(sev, is.numeric, logical(1))], c("SubjectID", "visit"))

theme_set(theme_bw(base_size = 11))
