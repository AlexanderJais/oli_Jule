# 17 - Export of the Olink results (all proteins, all samples) as CSV
# Out (output/export/):
#   proteins.csv                      one row per protein: OlinkID, Assay, UniProt, block, LOD, detection per matrix
#   samples.csv                       one row per sample: all metadata + Olink sample QC + whether excluded in QC
#   ISF_NPX_wide.csv, Serum_NPX_wide.csv                     NPX as delivered by Olink (samples x proteins)
#   ISF_PCNormalizedNPX_wide.csv, Serum_PCNormalizedNPX_wide.csv  the values used in the analysis
#   NPX_long.csv.gz                   everything in long format incl. LOD, below-LOD and QC flags (optional)
#   RELAD2/RELAD2_samples.csv, RELAD2/RELAD2_Serum_NPX_wide.csv, RELAD2/RELAD2_Serum_PCNormalizedNPX_wide.csv
# Values below LOD are exported as measured (Olink's recommendation); the long table flags them.
# Proteins Olink excluded (no values) are listed in proteins.csv but have no data columns.

source("R/utils.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
raw  <- read_step(cfg, "data", "npx_all_samples.rds", step = "scripts/02_import_qc.R")
sq   <- read_csv(file.path(cfg$paths$output, "qc", "sample_qc.csv"), show_col_types = FALSE)
det  <- read_csv(file.path(cfg$paths$output, "qc", "assay_detection.csv"), show_col_types = FALSE)
clear_outputs(cfg, "export")
sep  <- cfg$export$sep %||% ","

write_table <- function(df, ...) {
  p <- out_path(cfg, "export", ...)
  if (sep == ";") readr::write_csv2(df, p, na = "") else readr::write_csv(df, p, na = "")
  msg("  %s  (%d rows x %d columns)", file.path("export", ...), nrow(df), ncol(df))
}

# ---- samples ------------------------------------------------------------------------------------------------
clean_ids <- unique(readRDS(file.path(cfg$paths$output, "data", "npx_clean.rds"))$SampleID)
samples <- meta |>
  select(SampleID, SubjectID, SampleName, matrix, cohort, group, visit, visit_num, date, days_since_v1, site, state,
         lesion_state, clinical_state, relapse, relapse_raw, time_to_relapse, relapse_visit, dropout, plate, well, volume_ul,
         any_of(c("sex", "age", "BMI", "relad2_no")), flags) |>
  left_join(sq |> select(SampleID, SampleQC, qc_outlier = outlier), by = "SampleID") |>
  mutate(in_npx_data = SampleID %in% raw$SampleID, excluded_in_qc = in_npx_data & !SampleID %in% clean_ids)
write_table(samples, "samples.csv")

# ---- proteins -------------------------------------------------------------------------------------------------
proteins <- raw |>
  group_by(OlinkID, Assay, UniProt, across(any_of(c("Panel", "Block")))) |>
  summarise(n_values = sum(!is.na(NPX)), LOD_median = suppressWarnings(median(LOD, na.rm = TRUE)),
            normalization = paste(unique(Normalization), collapse = ";"), .groups = "drop") |>
  left_join(det |> select(matrix, OlinkID, frac_detected_all, keep) |>
              pivot_wider(names_from = matrix, values_from = c(frac_detected_all, keep)), by = "OlinkID") |>
  mutate(column_name = if_else(duplicated(Assay) | duplicated(Assay, fromLast = TRUE), paste(Assay, OlinkID, sep = "_"), Assay)) |>
  arrange(Assay)
write_table(proteins, "proteins.csv")

# ---- wide tables ------------------------------------------------------------------------------------------------
id_cols <- c("SampleID", "SubjectID", "matrix", "cohort", "group", "visit", "site", "state", "relapse", "plate",
             "SampleQC", "excluded_in_qc")
wide_of <- function(d, value) {
  w <- d |> filter(OlinkID %in% proteins$OlinkID[proteins$n_values > 0]) |>
    left_join(proteins |> select(OlinkID, column_name), by = "OlinkID") |>
    select(SampleID, column_name, v = all_of(value)) |>
    pivot_wider(names_from = column_name, values_from = v)
  w <- w |> select(SampleID, all_of(sort(setdiff(names(w), "SampleID"))))      # proteins alphabetically
  samples |> select(all_of(id_cols)) |> inner_join(w, by = "SampleID") |> arrange(SampleID)
}
raw_m <- raw |> left_join(meta |> select(SampleID, matrix, cohort), by = "SampleID")
for (mx in c("ISF", "Serum")) for (v in c("NPX", "PCNormalizedNPX")) {
  d <- raw_m |> filter(matrix == mx)
  if (!nrow(d)) next
  write_table(wide_of(d, v), sprintf("%s_%s_wide.csv", mx, v))
}

# ---- long table ----------------------------------------------------------------------------------------------------
if (isTRUE(cfg$export$long_format %||% TRUE)) {
  long <- raw |> left_join(samples |> select(all_of(id_cols[1:10])), by = "SampleID") |>
    relocate(all_of(id_cols[1:10])) |> arrange(matrix, SampleID, Assay)
  p <- out_path(cfg, "export", "NPX_long.csv.gz")
  if (sep == ";") readr::write_csv2(long, p, na = "") else readr::write_csv(long, p, na = "")
  msg("  export/NPX_long.csv.gz  (%d rows)", nrow(long))
}

# ---- RELAD2 ---------------------------------------------------------------------------------------------------------
rl <- samples |> filter(cohort == "RELAD2")
write_table(rl, "RELAD2", "RELAD2_samples.csv")
d <- raw_m |> filter(cohort == "RELAD2")
if (nrow(d)) for (v in c("NPX", "PCNormalizedNPX")) write_table(wide_of(d, v), "RELAD2", sprintf("RELAD2_Serum_%s_wide.csv", v))
msg("Export done: %s (%d RELAD2 samples)", file.path(cfg$paths$output, "export"), nrow(rl))
