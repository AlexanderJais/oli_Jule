# Shared helpers: packages, configuration, output files, logging.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

#' Load config.yml (or the file named in env var LEIP_CONFIG, e.g. the test configuration).
load_config <- function(path = Sys.getenv("LEIP_CONFIG", "config.yml")) {
  if (!file.exists(path)) stop("Config file not found: ", path)
  cfg <- yaml::read_yaml(path)
  cfg$config_file <- path
  cfg$paths <- resolve_data_paths(cfg$paths)
  cfg
}

# If a configured file is missing, look for it in the same folder under its usual name,
# so the original file names from Olink / the lab can be used without renaming.
data_file_patterns <- list(
  manifest      = "(manifest|Sample[ _-]?Submission[ _-]?Sheet).*\\.xlsx$",
  leip_clinical = "LEIP.*clinical.*\\.xlsx$",
  fixed_lod     = "Fixed[ _-]?LOD.*\\.csv$"        # e.g. "Explore HT_Fixed LOD.csv", "Explore_HT_Fixed_LOD.csv"
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

#' Path inside the output folder; creates the sub-folder if needed.
out_path <- function(cfg, ...) {
  p <- file.path(cfg$paths$output, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

msg <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", sprintf(...))

#' Remove results of an earlier run so a skipped analysis cannot leave stale files behind.
clear_outputs <- function(cfg, subdir, pattern = ".*") {
  d <- file.path(cfg$paths$output, subdir)
  if (!dir.exists(d)) return(invisible(0))
  f <- list.files(d, pattern = pattern, full.names = TRUE, recursive = TRUE)
  unlink(f)
  invisible(length(f))
}

#' Write a data frame as CSV (gzip-compressed if the name ends in .gz) and return it invisibly.
save_csv <- function(df, cfg, ...) {
  p <- out_path(cfg, ...)
  readr::write_csv(df, p, na = "")
  invisible(df)
}

save_plot <- function(p, cfg, ..., width = 8, height = 6) {
  f <- out_path(cfg, ...)
  # never exceed ggsave's 50-inch limit
  ggsave(f, p, width = min(width, 45), height = min(height, 45), dpi = 150)
  invisible(f)
}

#' Read an intermediate result written by an earlier script; fail with a clear message.
read_step <- function(cfg, ..., step) {
  p <- file.path(cfg$paths$output, ...)
  if (!file.exists(p)) stop(sprintf("%s not found - run %s first.", p, step))
  readRDS(p)
}

theme_set(theme_bw(base_size = 11))
