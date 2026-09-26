# Import and quality control of Olink Explore HT NPX data.

#' Read all parquet files in a folder with OlinkAnalyze and bind them.
import_npx <- function(npx_dir) {
  files <- list.files(npx_dir, pattern = "\\.parquet$", full.names = TRUE)
  if (!length(files)) stop("No .parquet files found in ", npx_dir)
  msg("Reading %d NPX file(s): %s", length(files), paste(basename(files), collapse = ", "))
  d <- map(files, \(f) OlinkAnalyze::read_npx(f) |> as_tibble() |> mutate(source_file = basename(f))) |>
    bind_rows()
  req <- c("SampleID", "SampleType", "PlateID", "OlinkID", "Assay", "AssayType", "NPX",
           "PCNormalizedNPX", "Normalization", "SampleQC", "AssayQC", "DataAnalysisRefID")
  miss <- setdiff(req, names(d))
  if (length(miss)) stop("NPX data is missing columns: ", paste(miss, collapse = ", "))
  d
}

#' LOD per assay. Uses Olink's fixed LOD file when available; otherwise the Olink
#' negative-control formula (median + max(0.2, 3 SD) of PC-normalised NPX) on all
#' negative controls pooled - OlinkAnalyze itself requires >= 10 NCs, this study has 8.
compute_lod <- function(d, fixed_lod_path = NULL, method = "auto") {
  nc <- d |>
    filter(SampleType == "NEGATIVE_CONTROL", AssayType == "assay", SampleQC != "FAIL") |>
    group_by(OlinkID, DataAnalysisRefID) |>
    summarise(n_nc = n_distinct(SampleID),
              LOD_NC = median(PCNormalizedNPX, na.rm = TRUE) + max(0.2, 3 * sd(PCNormalizedNPX, na.rm = TRUE)),
              .groups = "drop")

  fixed <- NULL
  if (method %in% c("auto", "FixedLOD") && !is.null(fixed_lod_path) && file.exists(fixed_lod_path)) {
    lf <- utils::read.table(fixed_lod_path, sep = ";", header = TRUE)
    fixed <- lf |> transmute(OlinkID, DataAnalysisRefID, LOD_fixed = LODNPX)
    if ("Version" %in% names(lf)) msg("Fixed LOD file version: %s", paste(unique(lf$Version), collapse = ", "))
  } else if (method == "FixedLOD") {
    stop("lod_method = FixedLOD but no fixed LOD file at ", fixed_lod_path)
  }

  lod <- nc
  if (!is.null(fixed)) lod <- full_join(lod, fixed, by = c("OlinkID", "DataAnalysisRefID"))
  else lod$LOD_fixed <- NA_real_
  lod |>
    mutate(LOD = coalesce(LOD_fixed, LOD_NC),
           LOD_source = case_when(!is.na(LOD_fixed) ~ "Olink fixed LOD",
                                  !is.na(LOD_NC) ~ sprintf("negative controls (n=%d)", n_nc),
                                  TRUE ~ "none"))
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
    summarise(max_group_frac = max(frac_detected), best_group = det_group[which.max(frac_detected)], .groups = "drop") |>
    left_join(overall, by = c("matrix", "OlinkID", "Assay")) |>
    mutate(keep = max_group_frac >= min_frac)
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
