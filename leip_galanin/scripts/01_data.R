# 01 - LEIP data: the Olink serum values of the LEIP samples (with LOD and QC) and their clinical data.
# The LEIP samples are the rows of the clinical file (column Olink_SampleID).
# In:  paths$npx_dir (parquet), paths$fixed_lod, paths$leip_clinical, optional paths$manifest
# Out: output/0_data/leip_data.rds and samples.csv, sample_checks.csv, detection.csv, parameters.csv

source("R/utils.R")
source("R/qc.R")                      # add_lod(): the Olink LOD, as in the main analysis
source("leip_galanin/R/leip.R")
cfg <- leip_config()
clear_outputs(cfg, "0_data")

# ---- clinical file (defines the LEIP samples) -----------------------------------------------------------------------
cl   <- read_leip_clinical(cfg$paths$leip_clinical, cfg$galanin$elisa_column %||% "Galanin [pg/mL]")
clin <- cl$data
ids  <- clin$SampleID
msg("Clinical file: %d LEIP samples, %d variables (%d in Key_parameters)", length(ids), ncol(clin) - 1, length(cl$key))

# ---- Olink values ----------------------------------------------------------------------------------------------------------
lp <- read_leip_npx(cfg$paths$npx_dir, ids, cfg$paths$fixed_lod, cfg$npx_column %||% "PCNormalizedNPX")
if (!"AssayQC" %in% names(lp)) lp$AssayQC <- NA_character_
sq <- lp |>
  group_by(SampleID, plate = PlateID) |>
  summarise(SampleQC = case_when(any(SampleQC == "FAIL") ~ "FAIL", any(SampleQC == "WARN") ~ "WARN", TRUE ~ "PASS"),
            median_npx = median(value, na.rm = TRUE), iqr_npx = IQR(value, na.rm = TRUE),
            frac_below_lod = mean(below_lod, na.rm = TRUE), .groups = "drop") |>
  mutate(ok = SampleQC != "FAIL",
         qc_outlier = abs(median_npx - mean(median_npx[ok])) > 3 * sd(median_npx[ok]) |
           abs(iqr_npx - mean(iqr_npx[ok])) > 3 * sd(iqr_npx[ok])) |>
  select(-ok)
drop <- if (isTRUE(cfg$drop_failed_samples %||% TRUE)) sq$SampleID[sq$SampleQC == "FAIL"] else character()
lp <- lp |> filter(!SampleID %in% drop)
empty <- lp |> group_by(OlinkID) |> summarise(any_value = any(!is.na(value))) |> filter(!any_value) |> pull(OlinkID)
lp <- lp |> filter(!OlinkID %in% empty)                      # e.g. assays Olink excluded (no values)
msg("Olink: %d LEIP samples (%d failed QC and are left out), %d assays (%d without values removed)",
    n_distinct(lp$SampleID), length(drop), n_distinct(lp$OlinkID), length(empty))

det <- lp |>
  group_by(OlinkID, Assay, across(any_of("UniProt"))) |>
  summarise(n_samples = sum(!is.na(value)),
            frac_above_lod = if (all(is.na(below_lod))) NA_real_ else mean(!below_lod, na.rm = TRUE),
            median_npx = median(value, na.rm = TRUE), median_lod = median(LOD, na.rm = TRUE),
            assay_qc_warn = any(AssayQC %in% "WARN"), .groups = "drop") |>
  mutate(measurable = n_samples >= 10 & coalesce(frac_above_lod >= (cfg$min_detect_frac %||% 0.5), TRUE),
         protein = if_else(duplicated(Assay) | duplicated(Assay, fromLast = TRUE), paste0(Assay, " (", OlinkID, ")"), Assay)) |>
  arrange(Assay)
W <- lp |> select(SampleID, OlinkID, value) |> pivot_wider(names_from = OlinkID, values_from = value) |> arrange(SampleID)
Y <- as.matrix(W[, -1]); rownames(Y) <- W$SampleID

# ---- one table per sample: Olink QC + clinical data ---------------------------------------------------------------------------
S <- tibble(SampleID = rownames(Y)) |>
  left_join(sq |> select(SampleID, plate, SampleQC, qc_outlier, frac_below_lod), by = "SampleID") |>
  left_join(clin, by = "SampleID")

# ---- checks: every clinical row, and the Olink samples ---------------------------------------------------------------------
key_vals <- setdiff(intersect(cl$key, names(clin)), leip_technical)
checks <- tibble(SampleID = ids, SubjectID = as.character(col_or(clin, "SubjectID")),
                 clinical_plate = as.character(col_or(clin, "Olink_plate")),
                 n_clinical_values = if (length(key_vals)) rowSums(!is.na(clin[key_vals])) else 0L,
                 low_serum = col_or(clin, "low_serum"), lipaemic = col_or(clin, "lipamisch")) |>
  left_join(sq |> select(SampleID, olink_plate = plate, SampleQC, qc_outlier), by = "SampleID")
if (!is.null(cfg$paths$manifest) && file.exists(cfg$paths$manifest)) {
  man <- readxl::read_excel(cfg$paths$manifest, sheet = "manifest", col_types = "text", na = c("", "EMPTY", "NA"))
  if (all(c("SampleID", "SubjectID") %in% names(man)))
    checks <- checks |> left_join(man |> select(SampleID, manifest_SubjectID = SubjectID), by = "SampleID")
}
# the plate grouping of the clinical file must match the Olink plates (the labels may differ, e.g. "Plate 3" / "Plate3")
pc <- checks |> filter(!is.na(clinical_plate), !is.na(olink_plate)) |> distinct(olink_plate, clinical_plate)
bad_plate <- (checks$olink_plate %in% pc$olink_plate[duplicated(pc$olink_plate)]) |
  (checks$clinical_plate %in% pc$clinical_plate[duplicated(pc$clinical_plate)])
issues <- list(
  if_else(is.na(checks$SampleQC), "not in the Olink data", NA_character_),
  if_else(checks$SampleQC %in% "FAIL", if (length(drop)) "Olink QC failed (left out)" else "Olink QC failed", NA_character_),
  if_else(checks$SampleQC %in% "WARN", "Olink QC warning (kept)", NA_character_),
  if_else(coalesce(checks$qc_outlier, FALSE), "Olink outlier: median/IQR beyond 3 SD (kept)", NA_character_),
  if_else(checks$n_clinical_values == 0, "no clinical values", NA_character_),
  if_else(bad_plate, "Olink plate differs from the clinical file", NA_character_),
  if_else(coalesce(checks$low_serum == 1, FALSE), "low serum volume (kept)", NA_character_),
  if_else(coalesce(checks$lipaemic == 1, FALSE), "lipaemic (kept)", NA_character_),
  if ("manifest_SubjectID" %in% names(checks))
    if_else(!is.na(checks$manifest_SubjectID) & checks$manifest_SubjectID != checks$SubjectID, "SubjectID differs from the manifest", NA_character_))
checks$problem <- pmap_chr(compact(issues), \(...) { v <- na.omit(c(...)); if (length(v)) paste(v, collapse = "; ") else NA_character_ })
if (any(!is.na(checks$problem))) {
  msg("Sample checks (0_data/sample_checks.csv):")
  print(checks |> filter(!is.na(problem)) |> select(SampleID, SubjectID, problem), n = Inf)
}

# ---- clinical variables to analyse ----------------------------------------------------------------------------------------------
pinfo <- select_parameters(clin |> filter(SampleID %in% S$SampleID), cl$key, unlist(cfg$clinical$exclude),
                           cfg$clinical$min_n %||% 20, cfg$clinical$min_group %||% 5)

gal <- find_assay(det, cfg$galanin$olink_assay %||% "GAL")
msg("%d of %d proteins measurable in LEIP serum (>= %.0f%% above LOD); %d of %d clinical variables analysed",
    sum(det$measurable), nrow(det), 100 * (cfg$min_detect_frac %||% 0.5), sum(pinfo$analysed), nrow(pinfo))
if (length(gal)) msg("Galanin: %s (%s), %.0f%% of LEIP samples above LOD", paste(det$Assay[det$OlinkID %in% gal], collapse = ", "),
                     paste(gal, collapse = ", "), 100 * det$frac_above_lod[det$OlinkID == gal[1]]) else
  msg("WARNING: galanin (%s) is not in the Olink data", cfg$galanin$olink_assay %||% "GAL")
if (!"galanin_elisa" %in% names(S)) msg("WARNING: no galanin ELISA column '%s' in the clinical file", cfg$galanin$elisa_column)

saveRDS(list(S = S, Y = Y, det = det, lod = lp |> select(SampleID, OlinkID, LOD, below_lod), pinfo = pinfo, key = cl$key,
             checks = checks), out_path(cfg, "0_data", "leip_data.rds"))
save_csv(S, cfg, "0_data", "samples.csv")
save_csv(checks, cfg, "0_data", "sample_checks.csv")
save_csv(det, cfg, "0_data", "detection.csv")
save_csv(pinfo, cfg, "0_data", "parameters.csv")
