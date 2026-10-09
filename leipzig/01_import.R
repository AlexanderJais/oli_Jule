# Leipzig analysis, part 1 (run with run_leipzig.R): read the Olink data of the Leipzig (LEIP) serum samples only,
# with LOD and sample QC as in the main pipeline. No other cohort is analysed.
# In:  manifest, NPX parquet file(s), LEIP clinical file (config.yml -> paths)
# Out: output/leipzig/data/leipzig_npx.rds, output/leipzig/qc/sample_qc.csv

source("R/utils.R")
source("R/metadata.R")
source("R/qc.R")
cfg  <- load_config()
vcol <- cfg$npx_column
clear_outputs(cfg, "leipzig")

meta <- build_metadata(cfg$paths$manifest, cfg$paths$leip_clinical) |> filter(cohort == "LEIP", matrix == "Serum")
if (!nrow(meta)) stop("No Leipzig serum samples (Study = LEIP) in the manifest.", call. = FALSE)
msg("%d Leipzig serum samples in the manifest; clinical data: %s", nrow(meta),
    if (file.exists(cfg$paths$leip_clinical %||% "")) cfg$paths$leip_clinical else "none (set leip_case$case)")

# only the NPX file(s) holding Leipzig samples; LOD per file, with that file's controls, as in the main pipeline
d <- import_npx(cfg$paths$npx_dir)
files <- unique(d$source_file[d$SampleID %in% meta$SampleID])
d <- d |> filter(source_file %in% files)
miss <- setdiff(meta$SampleID, d$SampleID)
if (length(miss)) msg("WARNING: %d Leipzig samples of the manifest are not in the NPX data: %s", length(miss), paste(miss, collapse = ", "))
d <- d |> group_split(source_file) |>
  map(\(x) add_lod(x, cfg$paths$fixed_lod, cfg$qc$lod_method, value_col = vcol)) |> bind_rows()

s <- d |> filter(SampleType == "SAMPLE", AssayType == "assay", SampleID %in% meta$SampleID) |>
  inner_join(meta, by = "SampleID") |>
  mutate(value = .data[[vcol]], below_lod = if_else(is.na(LOD), NA, value < LOD))
if (cfg$qc$drop_assay_qc_warn) s <- s |> mutate(value = if_else(AssayQC == "WARN", NA_real_, value))

sq <- sample_qc(s, cfg$qc$outlier_sd) |> left_join(meta |> select(SampleID, SubjectID, plate), by = "SampleID")
save_csv(sq, cfg, "leipzig", "qc", "sample_qc.csv")
drop <- sq |> filter(SampleQC == "FAIL" | (cfg$qc$drop_sample_qc_warn & SampleQC == "WARN") | (cfg$qc$drop_outliers & outlier)) |> pull(SampleID)
msg("Leipzig samples: %d FAIL, %d WARN, %d outliers -> %d excluded%s", sum(sq$SampleQC == "FAIL"), sum(sq$SampleQC == "WARN"),
    sum(sq$outlier, na.rm = TRUE), length(drop), if (length(drop)) paste0(" (", paste(sq$SubjectID[sq$SampleID %in% drop], collapse = ", "), ")") else "")
s <- s |> filter(!SampleID %in% drop)

# assays without any value (Olink EXCLUDED) are removed
empty <- s |> group_by(OlinkID) |> filter(all(is.na(value))) |> distinct(OlinkID) |> pull(OlinkID)
s <- s |> filter(!OlinkID %in% empty)
saveRDS(s, out_path(cfg, "leipzig", "data", "leipzig_npx.rds"))
msg("Leipzig data: %d samples, %d proteins (%d without values removed)", n_distinct(s$SampleID), n_distinct(s$OlinkID), length(empty))
