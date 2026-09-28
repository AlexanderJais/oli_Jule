# Reading the Olink data: the NPX parquet files and the limit of detection (LOD).
# add_lod() is the LOD routine of the O-MicroAD analysis, unchanged (tested on the real Olink files).

#' Olink rows of the LEIP samples: reads every NPX parquet file that contains one of `ids`, adds the
#' LOD (add_lod() below: Olink fixed LOD file, per sample for count-based assays), and returns the
#' sample x assay rows of these samples. Control IDs (PC1, NC1 ...) repeat on every plate and are made
#' unique per plate first.
read_leip_npx <- function(npx_dir, ids, fixed_lod = NULL, value_col = "PCNormalizedNPX") {
  files <- list.files(npx_dir, pattern = "\\.parquet$", full.names = TRUE)
  if (!length(files)) stop("No .parquet files in ", npx_dir)
  d <- map(files, \(f) {
    x <- OlinkAnalyze::read_npx(f) |> as_tibble()
    if (!any(x$SampleID %in% ids)) { msg("%s: no LEIP sample - skipped", basename(f)); return(NULL) }
    msg("%s: %d LEIP samples", basename(f), n_distinct(x$SampleID[x$SampleID %in% ids]))
    req <- c("SampleID", "SampleType", "PlateID", "OlinkID", "Assay", "AssayType", value_col, "SampleQC", "DataAnalysisRefID")
    miss <- setdiff(req, names(x))
    if (length(miss)) stop(basename(f), " is missing columns: ", paste(miss, collapse = ", "))
    x <- x |> mutate(SampleID = if_else(SampleType == "SAMPLE", SampleID, paste(SampleID, PlateID, sep = "@")))
    add_lod(x, fixed_lod, "auto", value_col = value_col) |> mutate(source_file = basename(f))
  }) |> compact()
  if (!length(d)) stop("None of the LEIP samples of the clinical file is in the NPX files in ", npx_dir)
  bind_rows(d) |>
    filter(SampleType == "SAMPLE", AssayType == "assay", SampleID %in% ids) |>
    mutate(value = .data[[value_col]], below_lod = if_else(is.na(LOD), NA, value < LOD))
}

#' Add a per-row LOD (PC-normalised NPX scale) to the full NPX data.
#' Uses OlinkAnalyze's own calculation, which is sample-specific for count-based assays
#' (LODMethod = lod_count: log2(LODCount / extension-control count) - plate-control median).
#'   1. Olink fixed LOD file (paths$fixed_lod), matched on DataAnalysisRefID
#'   2. rows without a fixed LOD: the Olink negative-control method on the pooled negative
#'      controls. OlinkAnalyze requires >= 10 NCs for this; this study has 8, so the internal
#'      routine is called with the lower minimum (less precise - reported as such).
#' The LOD is returned on the scale of `value_col`: PC-normalised, or - for intensity-normalised
#' NPX - shifted by the plate median of the samples, as OlinkAnalyze does.
add_lod <- function(d, fixed_lod_path = NULL, method = "auto", value_col = "PCNormalizedNPX") {
  d <- d |> mutate(.row = row_number(), LOD = NA_real_, LOD_source = NA_character_)
  use_fixed <- method %in% c("auto", "FixedLOD") && !is.null(fixed_lod_path) && file.exists(fixed_lod_path)
  if (method == "FixedLOD" && !use_fixed) stop("lod_method = FixedLOD but no fixed LOD file at ", fixed_lod_path)
  if (!use_fixed) msg("WARNING: no Olink fixed LOD file at %s - LOD from negative controls (less precise)", fixed_lod_path %||% "(not set)")

  if (use_fixed) {
    lf <- utils::read.table(fixed_lod_path, sep = ";", header = TRUE)
    msg("Fixed LOD file: version %s, %d DataAnalysisRefIDs", paste(unique(lf$Version), collapse = ", "),
        n_distinct(lf$DataAnalysisRefID))
    dar <- unique(d$DataAnalysisRefID)
    miss <- setdiff(dar, lf$DataAnalysisRefID)
    if (length(miss)) msg("WARNING: DataAnalysisRefID(s) %s not in the fixed LOD file - NC-based LOD used for them",
                          paste(miss, collapse = ", "))
    if (length(setdiff(dar, miss))) {
      fx <- suppressMessages(OlinkAnalyze::olink_lod(d |> select(-LOD, -LOD_source), lod_file_path = fixed_lod_path,
                                                     lod_method = "FixedLOD"))
      fx <- fx |> as_tibble() |> select(.row, LOD_fixed = PCNormalizedLOD)
      d <- d |> left_join(fx, by = ".row") |>
        mutate(LOD_source = if_else(!is.na(LOD_fixed), "Olink fixed LOD", LOD_source), LOD = LOD_fixed) |>
        select(-LOD_fixed)
    }
  }

  need_nc <- is.na(d$LOD) & d$AssayType == "assay"
  if (any(need_nc)) {
    ns <- asNamespace("OlinkAnalyze")
    n_nc <- n_distinct(d$SampleID[d$SampleType == "NEGATIVE_CONTROL" & d$SampleQC != "FAIL"])
    nc <- tryCatch({
      lod_data <- get("olink_nc_lod", ns)(d, min_num_nc = 2L)
      get("pc_norm_count", ns)(d |> select(-LOD, -LOD_source), lod_data) |>
        as_tibble() |> select(.row, LOD_nc = PCNormalizedLOD)
    }, error = \(e) {
      msg("OlinkAnalyze NC routine failed (%s); using median + max(0.2, 3 SD) of NC NPX", conditionMessage(e))
      d |> filter(SampleType == "NEGATIVE_CONTROL", AssayType == "assay", SampleQC != "FAIL") |>
        group_by(OlinkID, DataAnalysisRefID) |>
        summarise(LOD_nc = median(PCNormalizedNPX, na.rm = TRUE) + max(0.2, 3 * sd(PCNormalizedNPX, na.rm = TRUE)),
                  .groups = "drop") |>
        right_join(d |> select(.row, OlinkID, DataAnalysisRefID), by = c("OlinkID", "DataAnalysisRefID")) |>
        select(.row, LOD_nc)
    })
    d <- d |> left_join(nc, by = ".row") |>
      mutate(fill = is.na(LOD) & !is.na(LOD_nc),
             LOD = if_else(fill, LOD_nc, LOD),
             LOD_source = if_else(fill, sprintf("negative controls (n=%d)", n_nc), LOD_source)) |>
      select(-LOD_nc, -fill)
  }
  if (value_col == "NPX" && any(d$Normalization == "Intensity", na.rm = TRUE)) {
    pm <- d |> filter(SampleType == "SAMPLE", AssayType == "assay") |>
      group_by(OlinkID, PlateID) |>
      summarise(plate_median = median(PCNormalizedNPX, na.rm = TRUE), .groups = "drop")
    d <- d |> left_join(pm, by = c("OlinkID", "PlateID")) |>
      mutate(LOD = if_else(Normalization == "Intensity", LOD - plate_median, LOD)) |>
      select(-plate_median)
    msg("LOD converted to the intensity-normalised NPX scale (npx_column = NPX)")
  }
  d |> select(-.row)
}
