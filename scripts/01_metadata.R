# 01 - Sample metadata
# In:  manifest workbook (+ optional LEIP clinical file)
# Out: output/metadata/sample_metadata.{rds,csv}, output/metadata/data_flags.csv

source("R/utils.R")
source("R/metadata.R")
cfg <- load_config()

msg("Reading manifest: %s", cfg$paths$manifest)
meta <- build_metadata(cfg$paths$manifest, cfg$paths$leip_clinical)
layout <- read_plate_layout(cfg$paths$manifest)
validate_metadata(meta, layout)

flags <- flag_metadata(meta, cfg$qc$low_volume_ul)
meta <- meta |>
  left_join(flags |> group_by(SampleID) |> summarise(flags = paste(issue, collapse = "; ")),
            by = "SampleID")

saveRDS(meta, out_path(cfg, "metadata", "sample_metadata.rds"))
save_csv(meta, cfg, "metadata", "sample_metadata.csv")
save_csv(flags |> left_join(meta |> select(SampleID, SubjectID, cohort, matrix, visit), by = "SampleID") |>
           relocate(issue, .after = last_col()),
         cfg, "metadata", "data_flags.csv")

msg("%d samples: %s", nrow(meta), paste(names(table(meta$matrix)), table(meta$matrix), collapse = ", "))
print(count(meta, matrix, cohort, group))
msg("%d flagged sample/issue pairs:", nrow(flags))
print(count(flags, issue))
