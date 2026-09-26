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
  cfg
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
  ggsave(f, p, width = width, height = height, dpi = 150)
  invisible(f)
}

#' Read an intermediate result written by an earlier script; fail with a clear message.
read_step <- function(cfg, ..., step) {
  p <- file.path(cfg$paths$output, ...)
  if (!file.exists(p)) stop(sprintf("%s not found - run %s first.", p, step))
  readRDS(p)
}

theme_set(theme_bw(base_size = 11))
