# Import and quality control of Olink Explore HT NPX data.

#' Matrix named in an NPX file name ("ISF" / "Serum"), NA if neither.
file_matrix <- function(f) {
  f <- toupper(basename(f))
  case_when(str_detect(f, "ISF") ~ "ISF", str_detect(f, "SERUM") ~ "Serum", TRUE ~ NA_character_)
}

#' Read all parquet files in a folder with OlinkAnalyze and bind them.
#' Olink may deliver one file per matrix (e.g. O-MicroAD_ISF_NPX_*.parquet, O-MicroAD_Serum_NPX_*.parquet).
#' Control wells of a mixed plate then appear in both files; they are kept per file (LOD is
#' computed within each file) and de-duplicated only where controls are pooled.
import_npx <- function(npx_dir) {
  files <- list.files(npx_dir, pattern = "\\.parquet$", full.names = TRUE)
  if (!length(files)) stop("No .parquet files found in ", npx_dir)
  msg("Reading %d NPX file(s): %s", length(files), paste(basename(files), collapse = ", "))
  d <- map(files, \(f) OlinkAnalyze::read_npx(f) |> as_tibble() |>
             mutate(source_file = basename(f), file_matrix = file_matrix(f))) |>
    bind_rows()
  dup <- d |> filter(SampleType == "SAMPLE") |> distinct(SampleID, source_file) |> count(SampleID) |> filter(n > 1)
  if (nrow(dup)) stop("Sample(s) present in more than one NPX file: ", paste(head(dup$SampleID, 10), collapse = ", "))
  req <- c("SampleID", "SampleType", "PlateID", "OlinkID", "Assay", "AssayType", "NPX",
           "PCNormalizedNPX", "Normalization", "SampleQC", "AssayQC", "DataAnalysisRefID")
  miss <- setdiff(req, names(d))
  if (length(miss)) stop("NPX data is missing columns: ", paste(miss, collapse = ", "))
  d
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

#' One row per assay: LOD summary across samples (count-based LODs vary by sample).
lod_summary <- function(d) {
  d |> filter(AssayType == "assay", SampleType == "SAMPLE") |>
    group_by(OlinkID, Assay, DataAnalysisRefID) |>
    summarise(LOD_source = paste(unique(na.omit(LOD_source)), collapse = "; "),
              LOD_median = median(LOD, na.rm = TRUE), LOD_min = min(LOD, na.rm = TRUE),
              LOD_max = max(LOD, na.rm = TRUE), .groups = "drop") |>
    mutate(across(starts_with("LOD_m"), \(x) if_else(is.finite(x), x, NA_real_)))
}

#' Per-sample QC summary with the Olink outlier rule (median / IQR beyond mean +/- k SD, per matrix).
sample_qc <- function(s, k = 3) {
  s |>
    group_by(SampleID, matrix) |>
    summarise(SampleQC = case_when(any(SampleQC == "FAIL") ~ "FAIL", any(SampleQC == "WARN") ~ "WARN", TRUE ~ "PASS"),
              median_npx = median(value, na.rm = TRUE),
              iqr_npx = IQR(value, na.rm = TRUE),
              frac_below_lod = mean(below_lod, na.rm = TRUE),
              .groups = "drop") |>
    group_by(matrix) |>
    mutate(outlier = abs(median_npx - mean(median_npx[SampleQC != "FAIL"])) > k * sd(median_npx[SampleQC != "FAIL"]) |
                     abs(iqr_npx - mean(iqr_npx[SampleQC != "FAIL"])) > k * sd(iqr_npx[SampleQC != "FAIL"])) |>
    ungroup()
}

#' Detection group used for the assay filter: ISF by group and skin state, serum by group.
detect_group <- function(matrix, group, state) {
  if_else(matrix == "ISF", paste(group, coalesce(state, ""), sep = ":"), group)
}

#' Fraction of samples above LOD per assay, matrix and detection group; keep flag.
assay_detection <- function(s, min_frac = 0.5) {
  by_group <- s |>
    group_by(matrix, OlinkID, Assay, det_group) |>
    summarise(n = sum(!is.na(value)), frac_detected = mean(!below_lod, na.rm = TRUE), .groups = "drop")
  overall <- s |>
    group_by(matrix, OlinkID, Assay) |>
    summarise(frac_detected_all = mean(!below_lod, na.rm = TRUE), .groups = "drop")
  by_group |>
    group_by(matrix, OlinkID, Assay) |>
    summarise(max_group_frac = suppressWarnings(max(frac_detected, na.rm = TRUE)),
              best_group = if (all(is.nan(frac_detected))) NA_character_ else det_group[which.max(frac_detected)],
              .groups = "drop") |>
    left_join(overall, by = c("matrix", "OlinkID", "Assay")) |>
    # assays without any LOD cannot be judged: kept, and marked lod_available = FALSE
    mutate(lod_available = is.finite(max_group_frac),
           keep = !lod_available | max_group_frac >= min_frac)
}

#' CV of the sample controls (linear scale), across all plates and within plates.
control_cv <- function(d, value_col) {
  sc <- d |> filter(SampleType == "SAMPLE_CONTROL", AssayType == "assay") |>
    mutate(lin = 2^.data[[value_col]])
  inter <- sc |> group_by(OlinkID, Assay) |>
    summarise(inter_cv = sd(lin, na.rm = TRUE) / mean(lin, na.rm = TRUE), .groups = "drop")
  intra <- sc |> group_by(OlinkID, Assay, PlateID) |>
    summarise(cv = sd(lin, na.rm = TRUE) / mean(lin, na.rm = TRUE), .groups = "drop") |>
    group_by(OlinkID, Assay) |> summarise(intra_cv = mean(cv, na.rm = TRUE), .groups = "drop")
  left_join(inter, intra, by = c("OlinkID", "Assay"))
}

#' Wide matrix (assays x samples) for one matrix type from the cleaned long data.
wide_matrix <- function(s, which_matrix, assays_keep) {
  w <- s |>
    filter(matrix == which_matrix, OlinkID %in% assays_keep) |>
    select(OlinkID, SampleID, value) |>
    pivot_wider(names_from = SampleID, values_from = value)
  m <- as.matrix(w[, -1]); rownames(m) <- w$OlinkID
  m
}

#' PCA scores (samples) from an assays x samples matrix; missing values set to the assay median.
pca_scores <- function(m, n = 4) {
  m <- t(apply(m, 1, \(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))
  m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
  p <- prcomp(t(m), center = TRUE, scale. = TRUE)
  ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  sc <- as_tibble(p$x[, seq_len(min(n, ncol(p$x)))], rownames = "SampleID")
  attr(sc, "var_explained") <- ve[seq_len(min(n, length(ve)))]
  sc
}
