# 02 - Import NPX data and quality control
# In:  parquet file(s) in paths$npx_dir, output of 01
# Out: output/data/npx_clean.rds (long), output/data/npx_wide.rds (per matrix),
#      output/qc/*.csv and plots

source("R/utils.R")
source("R/qc.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
vcol <- cfg$npx_column

d <- import_npx(cfg$paths$npx_dir)

# ---- normalisation check -------------------------------------------------------------
norm <- count(d |> filter(SampleType == "SAMPLE"), Normalization)
save_csv(norm, cfg, "qc", "normalization.csv")
msg("Normalization in file: %s", paste(norm$Normalization, norm$n, collapse = ", "))
if (any(norm$Normalization == "Intensity") && vcol == "NPX")
  warning("Data are intensity normalised and npx_column = NPX. With matrix-separated plates ",
          "PCNormalizedNPX is recommended (see config.yml).")
msg("Analysing column: %s", vcol)

# ---- sample IDs vs manifest ------------------------------------------------------------
ids_npx <- unique(d$SampleID[d$SampleType == "SAMPLE"])
id_check <- bind_rows(
  tibble(SampleID = setdiff(ids_npx, meta$SampleID), problem = "in NPX data, not in manifest"),
  tibble(SampleID = setdiff(meta$SampleID, ids_npx), problem = "in manifest, not in NPX data")
)
save_csv(id_check, cfg, "qc", "sample_id_mismatches.csv")
if (nrow(id_check)) msg("WARNING: %d sample ID mismatches - see qc/sample_id_mismatches.csv", nrow(id_check))

# ---- LOD ---------------------------------------------------------------------------------
d <- add_lod(d, cfg$paths$fixed_lod, cfg$qc$lod_method)
lod <- lod_summary(d)
save_csv(lod, cfg, "qc", "lod.csv")
msg("LOD source (assays): %s", paste(names(table(lod$LOD_source)), table(lod$LOD_source), collapse = ", "))

# ---- control samples ----------------------------------------------------------------------
cv <- control_cv(d, vcol)
save_csv(cv, cfg, "qc", "sample_control_cv.csv")
msg("Sample controls: median inter-plate CV %.1f%%, intra-plate CV %.1f%%",
    100 * median(cv$inter_cv, na.rm = TRUE), 100 * median(cv$intra_cv, na.rm = TRUE))

# ---- study samples --------------------------------------------------------------------------
s <- d |>
  filter(SampleType == "SAMPLE", AssayType == "assay") |>
  inner_join(meta, by = "SampleID") |>
  mutate(value = .data[[vcol]],
         below_lod = !is.na(LOD) & value < LOD,
         det_group = detect_group(matrix, group, state))

if (cfg$qc$drop_assay_qc_warn) s <- s |> mutate(value = if_else(AssayQC == "WARN", NA_real_, value))

sq <- sample_qc(s, cfg$qc$outlier_sd) |>
  left_join(meta |> select(SampleID, SubjectID, cohort, group, plate, volume_ul, flags), by = "SampleID")
save_csv(sq, cfg, "qc", "sample_qc.csv")

drop <- sq |>
  filter(SampleQC == "FAIL" |
           (cfg$qc$drop_sample_qc_warn & SampleQC == "WARN") |
           (cfg$qc$drop_outliers & outlier)) |>
  pull(SampleID)
msg("Samples: %d FAIL, %d WARN, %d outliers -> %d excluded",
    sum(sq$SampleQC == "FAIL"), sum(sq$SampleQC == "WARN"), sum(sq$outlier, na.rm = TRUE), length(drop))
s <- s |> filter(!SampleID %in% drop)

# ---- assay detection filter (per matrix) -----------------------------------------------------
det <- assay_detection(s, cfg$qc$min_detect_frac)
save_csv(det, cfg, "qc", "assay_detection.csv")
print(det |> group_by(matrix) |> summarise(assays = n(), kept = sum(keep)))

aqc <- s |> filter(AssayQC != "PASS") |> distinct(OlinkID, Assay, PlateID, AssayQC)
save_csv(aqc, cfg, "qc", "assay_qc_flags.csv")

keep <- split(det$OlinkID[det$keep], det$matrix[det$keep])
s <- s |> left_join(det |> select(matrix, OlinkID, keep), by = c("matrix", "OlinkID"))

# ---- plots ---------------------------------------------------------------------------------------
p <- ggplot(sq, aes(median_npx, iqr_npx, colour = SampleQC, shape = outlier)) +
  geom_point() + facet_wrap(~matrix, scales = "free") +
  labs(title = "Sample QC: median vs IQR of NPX", x = "median NPX", y = "IQR NPX")
save_plot(p, cfg, "qc", "sample_qc_median_iqr.png", width = 10, height = 5)

p <- det |> ggplot(aes(frac_detected_all, fill = matrix)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = cfg$qc$min_detect_frac, linetype = 2) +
  labs(title = "Share of samples above LOD per assay", x = "fraction above LOD", y = "assays")
save_plot(p, cfg, "qc", "assay_detection.png")

pc <- map(c("ISF", "Serum"), \(mx) {
  if (!length(keep[[mx]])) return(NULL)
  sc <- pca_scores(wide_matrix(s, mx, keep[[mx]]))
  ve <- attr(sc, "var_explained")
  sc |> left_join(meta, by = "SampleID") |>
    mutate(label = sprintf("%s  (PC1 %.0f%%, PC2 %.0f%%)", mx, ve[1], ve[2]))
}) |> bind_rows()
p <- ggplot(pc, aes(PC1, PC2, colour = plate, shape = group)) + geom_point() +
  facet_wrap(~label, scales = "free") + labs(title = "PCA per matrix, coloured by plate")
save_plot(p, cfg, "qc", "pca_by_plate.png", width = 11, height = 5)

lv <- s |> filter(!is.na(volume_ul)) |> distinct(SampleID, matrix, volume_ul) |>
  left_join(sq |> select(SampleID, median_npx, frac_below_lod), by = "SampleID") |>
  mutate(low_volume = volume_ul < unlist(cfg$qc$low_volume_ul)[matrix])
save_csv(lv |> group_by(matrix, low_volume) |>
           summarise(n = n(), median_npx = median(median_npx), frac_below_lod = median(frac_below_lod), .groups = "drop"),
         cfg, "qc", "low_volume_summary.csv")

# ---- save --------------------------------------------------------------------------------------------
saveRDS(s, out_path(cfg, "data", "npx_clean.rds"))
saveRDS(map(c(ISF = "ISF", Serum = "Serum"), \(mx) wide_matrix(s, mx, keep[[mx]])),
        out_path(cfg, "data", "npx_wide.rds"))
msg("Done: %d samples, ISF %d assays, serum %d assays kept.",
    n_distinct(s$SampleID), length(keep$ISF), length(keep$Serum))
