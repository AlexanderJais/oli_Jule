# 02b - QC overview: how many Olink proteins are above / below LOD in serum and dISF, overall and per protein class
# In:  output of steps 01-02; optional Human Protein Atlas table (paths$hpa, tools/download_hpa.R) for the protein classes
# Out: output/qc_overview/qc_overview.pdf + .png (ring charts + class bars), qc_overview.xlsx, csv tables, plots.rds
# Settings: config.yml -> qc_overview (rule, serum cohorts, font, colours); the threshold is qc$min_detect_frac.

source("R/utils.R")
source("R/report.R")
source("R/qc_overview.R")
cfg   <- load_config()
clear_outputs(cfg, "qc_overview")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
det   <- read_csv(file.path(cfg$paths$output, "qc", "assay_detection.csv"), show_col_types = FALSE)
qo    <- cfg$qc_overview %||% list()
rule  <- qo$rule %||% "all_samples"
min_f <- cfg$qc$min_detect_frac %||% 0.5
serum_cohorts <- unlist(qo$serum_cohorts %||% "MicroAD")
font  <- qo$font %||% "Nimbus Sans"
colours <- qo$colours %||% list(above = "#69005F", below = "#FF506E", no_lod = "#BDBDBD")
if (!rule %in% c("all_samples", "analysis_filter")) stop("qc_overview$rule must be all_samples or analysis_filter, not ", rule)
avail <- sort(unique(clean$cohort[clean$matrix == "Serum"]))
bad <- setdiff(serum_cohorts, avail)
if (length(bad)) stop("qc_overview$serum_cohorts: no serum samples for ", paste(bad, collapse = ", "),
                      ". Available serum cohorts: ", paste(avail, collapse = ", "), " ([] = all).", call. = FALSE)
# the step 02 filter (rule analysis_filter) is computed over all serum cohorts, so the serum set follows the rule
fig_cohorts <- if (rule == "analysis_filter") character() else serum_cohorts
serum_lab <- if (length(fig_cohorts)) paste(fig_cohorts, collapse = " + ") else "all cohorts"
def_all <- sprintf(">= %.0f%% of all samples", 100 * min_f)
def_filter <- sprintf("analysis filter of step 02 (>= %.0f%% in at least one group)", 100 * min_f)

st <- lod_status(clean, det, rule, min_f, fig_cohorts)
msg("Proteins above LOD (%s, >= %.0f%%; serum: %s): %s", rule, 100 * min_f, serum_lab,
    paste(st |> count(matrix, status) |> mutate(t = paste(matrix, status, n)) |> pull(t), collapse = "; "))

# ---- protein classes (Human Protein Atlas) ------------------------------------------------------------------------
hpa <- read_hpa(cfg$paths$hpa)
cls <- NULL; not_in_hpa <- tibble()
if (!is.null(hpa)) {
  ac <- assay_protein_classes(st |> distinct(OlinkID, Assay, UniProt), hpa)
  cls <- if (nrow(ac$classes)) ac$classes else NULL
  not_in_hpa <- st |> distinct(OlinkID, Assay, UniProt) |> filter(!OlinkID %in% ac$matched)
  msg("Human Protein Atlas (downloaded %s): %d of %d assays matched; %d not in the HPA", attr(hpa, "downloaded"),
      length(ac$matched), n_distinct(st$OlinkID), nrow(not_in_hpa))
  if (is.null(cls)) msg("WARNING: the Human Protein Atlas file %s gives no protein classes (column 'Protein class' missing or no assay matched)", cfg$paths$hpa)
} else msg("No Human Protein Atlas file at %s - protein classes skipped (run source(\"tools/download_hpa.R\") once).", cfg$paths$hpa %||% "(not set)")

# ---- tables ---------------------------------------------------------------------------------------------------------
overall <- st |> count(matrix, status) |> group_by(matrix) |> mutate(total = sum(n), pct = round(100 * n / total, 1)) |> ungroup()
by_class <- if (is.null(cls)) tibble() else
  st |> filter(status != "no LOD") |> inner_join(cls, by = "OlinkID", relationship = "many-to-many") |>
  count(matrix, class, status) |> complete(nesting(matrix, class), status = c("above LOD", "below LOD"), fill = list(n = 0)) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0) |>
  mutate(total = `above LOD` + `below LOD`, pct_above_LOD = round(100 * `above LOD` / total, 1),
         class = factor(class, names(qc_protein_classes))) |> arrange(matrix, class) |> mutate(class = as.character(class))
# every definition next to each other, so each number elsewhere in the outputs can be traced
cfg_lab <- if (length(serum_cohorts)) paste(serum_cohorts, collapse = " + ") else "all cohorts"
variants <- bind_rows(
  lod_status(clean, det, "all_samples", min_f, serum_cohorts) |> mutate(definition = def_all, serum = cfg_lab),
  if (length(serum_cohorts)) lod_status(clean, det, "all_samples", min_f, character()) |> filter(matrix == "Serum") |>
    mutate(definition = def_all, serum = "all cohorts"),
  lod_status(clean, det, "analysis_filter", min_f, character()) |> mutate(definition = def_filter, serum = "all cohorts")) |>
  count(definition, serum, matrix, status) |> pivot_wider(names_from = status, values_from = n, values_fill = 0) |>
  mutate(used_in_figure = definition == (if (rule == "all_samples") def_all else def_filter) & (matrix == "dISF" | serum == serum_lab))
proteins <- st |> select(matrix, OlinkID, Assay, UniProt, status, frac_above_LOD, n_samples, analysis_filter) |>
  pivot_wider(names_from = matrix, values_from = c(status, frac_above_LOD, n_samples, analysis_filter), names_glue = "{matrix}_{.value}")
if (!is.null(cls)) proteins <- proteins |>
  left_join(cls |> mutate(v = TRUE) |> pivot_wider(names_from = class, values_from = v, values_fill = FALSE), by = "OlinkID") |>
  mutate(across(any_of(names(qc_protein_classes)), \(x) coalesce(x, FALSE)), in_HPA = !OlinkID %in% not_in_hpa$OlinkID)

save_csv(overall, cfg, "qc_overview", "detection_overall.csv")
if (nrow(by_class)) save_csv(by_class, cfg, "qc_overview", "detection_by_class.csv")
save_csv(proteins, cfg, "qc_overview", "protein_detection.csv")
readme <- tibble(item = c("figure", "above LOD", "samples", "proteins", "protein classes", "class definitions", "not in HPA", "sheets"),
                 note = c("qc_overview.pdf / .png: ring charts (share of Olink proteins above / below LOD, total in the middle) and the distribution across protein classes (serum left, dISF right)",
                          if (rule == "all_samples") sprintf("a protein is above LOD in a matrix if >= %.0f%% of the samples of that matrix are above its LOD (qc$min_detect_frac); 'no LOD' = no LOD available", 100 * min_f)
                          else sprintf("the detection filter of step 02: >= %.0f%% above LOD in at least one group", 100 * min_f),
                          sprintf("dISF: all dISF samples passing QC. Serum: %s, samples passing QC%s", serum_lab,
                                  if (rule == "analysis_filter") " (the step 02 filter covers all serum cohorts)" else ""),
                          "one Olink assay = one protein; assays without any value (Olink EXCLUDED) are not counted",
                          if (is.null(hpa)) "not available: no Human Protein Atlas file (run source(\"tools/download_hpa.R\"))" else
                          if (is.null(cls)) "not available: the Human Protein Atlas table gives no 'Protein class' values for these assays" else
                            sprintf("Human Protein Atlas, column 'Protein class' (proteinatlas.org, file downloaded %s; CC BY 4.0). A protein can belong to several classes, so the bars do not add up to the total.", attr(hpa, "downloaded")),
                          paste(sprintf("%s = %s", names(qc_protein_classes), map_chr(qc_protein_classes, paste, collapse = " + ")), collapse = "; "),
                          if (is.null(hpa)) "not applicable: no Human Protein Atlas file" else if (!nrow(not_in_hpa)) "all assays matched to the HPA (by UniProt, gene symbol or synonym)" else
                            sprintf("%d assays could not be matched to the HPA (by UniProt, gene symbol or synonym); listed in sheet not_in_HPA", nrow(not_in_hpa)),
                          paste0("overall = ring charts", if (nrow(by_class)) "; by_class = bars" else "",
                                 "; definitions = counts under each above-LOD definition (used_in_figure marks the figure's); proteins = one row per assay",
                                 if (nrow(not_in_hpa)) "; not_in_HPA = assays not found in the HPA" else "")))
writexl::write_xlsx(Filter(\(x) nrow(x) > 0, list(README = readme, overall = overall, by_class = by_class, definitions = variants,
                                                 proteins = proteins, not_in_HPA = not_in_hpa)),
                    out_path(cfg, "qc_overview", "qc_overview.xlsx"))

# ---- figure -----------------------------------------------------------------------------------------------------------
plots <- qc_overview_plots(st, cls, colours)
note <- if (!is.null(cls)) NULL else if (is.null(hpa))
  "Protein classes: Human Protein Atlas table not found - run source(\"tools/download_hpa.R\") once, then rerun from step 02b (start_at <- \"02b\")." else
  "Protein classes: no Olink assay could be given a class from the Human Protein Atlas table - check the file and its 'Protein class' column."
methods <- sprintf("Detection overview (step 02b): above LOD = %s. Samples: %s. Protein classes: %s", readme$note[2], readme$note[3], readme$note[5])
saveRDS(list(plots = plots, note = note, methods = methods), out_path(cfg, "qc_overview", "plots.rds"))
fi <- qc_font_setup(font)
if (fi$mode != "showtext" && grepl("^nimbus ?sans", font, ignore.case = TRUE) && file.exists("fonts/NimbusSans-Regular.otf"))
  stop("Package 'showtext' is needed to draw the figure in Nimbus Sans (fonts/): run source(\"install_packages.R\"), ",
       "restart R, then start_at <- \"02b\"; source(\"run_all.R\").", call. = FALSE)
for (type in c("pdf", "png")) {
  open_device(out_path(cfg, "qc_overview", paste0("qc_overview.", type)), type, fi)
  draw_qc_overview(plots, no_class_note = note, family = fi$family)
  close_device(fi)
}
msg("Font: %s (%s)", font, switch(fi$mode, showtext = "from fonts/, drawn as outlines - no font listed in the PDF, text not selectable",
                                   cairo = "installed font, embedded", builtin = "Helvetica instead"))
if (!is.null(fi$message)) msg("WARNING: %s", fi$message)
msg("QC overview: %s", file.path(cfg$paths$output, "qc_overview"))
