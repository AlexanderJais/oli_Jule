# 11 - Aim 3: longitudinal trajectories across the disease course (MicroAD, AD patients)
# Disease course of the tracked lesion: lesional (V1) -> cleared (ex-lesional) -> relapse or not.
#   A  normalisation after clearance: does the residual lesional signal in cleared skin
#      (ex-lesional minus non-lesional, same visit) fade with weeks since clearance?
#   B  approach to relapse (relapsers only, exploratory): does it rise in the weeks before relapse?
#   C  the same for serum at cleared visits
#   D  severity (optional): proteins tracking clinical scores over visits, if paths$severity is set
# plus per-patient trajectory plots (dISF lesional site, dISF non-lesional site, serum) for top proteins.
# Out: output/trajectories/*

source("R/utils.R")
source("R/models.R")
cfg   <- load_config()
clear_outputs(cfg, "trajectories")
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)

ad <- meta |> filter(cohort == "MicroAD", group == "AD")

# disease-course anchors per patient: date the tracked lesion was first cleared, date of relapse
anchors <- ad |>
  filter(matrix == "ISF", site == "L") |>
  group_by(SubjectID) |>
  summarise(clear_date = suppressWarnings(min(date[state %in% "ex-lesional"], na.rm = TRUE)),
            relapse_date = if (any(!is.na(relapse_visit))) min(date[visit_num == relapse_visit[1]]) else as.Date(NA),
            .groups = "drop") |>
  mutate(clear_date = if_else(is.finite(clear_date), clear_date, as.Date(NA)))
ad <- ad |> left_join(anchors, by = "SubjectID") |>
  mutate(weeks_since_clear = as.numeric(date - clear_date) / 7,
         weeks_to_relapse = as.numeric(date - relapse_date) / 7)   # negative = before relapse

# ex-lesional minus non-lesional per visit (removes plate and day-to-day systemic variation)
isf_pairs <- ad |>
  filter(matrix == "ISF", SampleID %in% colnames(wide$ISF)) |>
  select(SubjectID, visit, site, state, SampleID) |>
  pivot_wider(id_cols = c(SubjectID, visit), names_from = site, values_from = c(SampleID, state)) |>
  filter(!is.na(SampleID_L), !is.na(SampleID_NL), state_L %in% "ex-lesional") |>
  left_join(ad |> select(SampleID_L = SampleID, weeks_since_clear, weeks_to_relapse, relapse), by = "SampleID_L") |>
  mutate(SampleID = paste(SubjectID, visit, sep = "_"))
delta <- wide$ISF[, isf_pairs$SampleID_L, drop = FALSE] - wide$ISF[, isf_pairs$SampleID_NL, drop = FALSE]
colnames(delta) <- isf_pairs$SampleID

results <- list()
results$A <- fit_contrasts(delta, isf_pairs |> filter(relapse != "dropout" | is.na(relapse)),
                           ~ weeks_since_clear + (1 | SubjectID), c(per_week = "weeks_since_clear"),
                           "ISF_xL_minus_NL_vs_weeks_since_clearance", cfg$stats$min_group_n)
results$B <- fit_contrasts(delta, isf_pairs |> filter(!is.na(weeks_to_relapse)),
                           ~ weeks_to_relapse + (1 | SubjectID), c(per_week = "weeks_to_relapse"),
                           "ISF_xL_minus_NL_vs_weeks_to_relapse", cfg$stats$min_group_n)
serum_info <- ad |> filter(matrix == "Serum", SampleID %in% colnames(wide$Serum), lesion_state == "cleared",
                           !is.na(weeks_to_relapse))
results$C <- fit_contrasts(wide$Serum, serum_info, ~ weeks_to_relapse + plate + (1 | SubjectID),
                           c(per_week = "weeks_to_relapse"), "Serum_vs_weeks_to_relapse", cfg$stats$min_group_n)

# ---- D: clinical severity (optional) -------------------------------------------------------------------
sev_path <- cfg$paths$severity
if (!is.null(sev_path) && file.exists(sev_path)) {
  sev <- if (str_detect(sev_path, "\\.xlsx?$")) readxl::read_excel(sev_path) else read_csv(sev_path, show_col_types = FALSE)
  names(sev)[tolower(names(sev)) == "subjectid"] <- "SubjectID"
  names(sev)[tolower(names(sev)) == "visit"] <- "visit"
  names(sev) <- make.names(names(sev))            # e.g. "itch NRS" -> "itch.NRS"
  scores <- setdiff(names(sev)[vapply(sev, is.numeric, logical(1))], c("SubjectID", "visit"))
  msg("Severity file: %d rows, scores: %s", nrow(sev), paste(scores, collapse = ", "))
  sev_info <- ad |> left_join(sev, by = c("SubjectID", "visit"))
  for (sc in scores) {
    for (what in c("ISF lesional site", "ISF non-lesional site", "Serum")) {
      ii <- sev_info |> filter(!is.na(.data[[sc]]),
                               case_when(what == "Serum" ~ matrix == "Serum",
                                         what == "ISF lesional site" ~ matrix == "ISF" & site == "L",
                                         TRUE ~ matrix == "ISF" & site == "NL"))
      ex <- if (what == "Serum") wide$Serum else wide$ISF
      results[[paste(sc, what)]] <- fit_contrasts(
        ex, ii, as.formula(sprintf("~ %s + plate + (1 | SubjectID)", sc)), c(per_unit = sc),
        sprintf("%s: %s", what, sc), cfg$stats$min_group_n)
    }
  }
} else msg("No severity file (paths$severity) - severity models skipped.")

res <- bind_rows(results)
if (nrow(res)) {
  res <- annotate_results(res, assay_map, cfg$stats$fdr)
  save_csv(res, cfg, "trajectories", "trajectory_results.csv")
  summ <- res |> group_by(model, n_samples, n_subjects) |>
    summarise(n_sig = sum(significant), n_up = sum(significant & logFC > 0), n_down = sum(significant & logFC < 0),
              .groups = "drop")
  save_csv(summ, cfg, "trajectories", "trajectory_summary.csv")
  print(as.data.frame(summ))
}

# ---- per-patient trajectory plots for top proteins -------------------------------------------------------
tops <- character()
take <- \(f, mdl, ct, n) {
  p <- file.path(cfg$paths$output, "models", f)
  if (!file.exists(p)) return(character())
  read_csv(p, show_col_types = FALSE) |> filter(model == mdl, contrast == ct, significant) |>
    slice_min(P.Value, n = n, with_ties = FALSE) |> pull(OlinkID)
}
tops <- c(take("ISF_results.csv", "states_all_visits", "AD_L_vs_NL", 6),
          take("ISF_relapse_delta_results.csv", "relapse_delta_xL_minus_NL", "relapse_vs_non", 3),
          take("Serum_results.csv", "MicroAD_active_vs_cleared", "active_vs_cleared", 3))
if (nrow(res)) tops <- c(tops, res |> filter(significant) |> slice_min(P.Value, n = 3, with_ties = FALSE) |> pull(OlinkID))
tops <- unique(tops)

if (length(tops)) {
  long <- clean |>
    filter(OlinkID %in% tops, cohort == "MicroAD", group == "AD") |>
    left_join(anchors, by = "SubjectID") |>
    mutate(series = case_when(matrix == "Serum" ~ "serum",
                              site == "L" ~ "dISF lesional site",
                              TRUE ~ "dISF non-lesional site"),
           point_state = if_else(matrix == "ISF" & site == "L", state, "(other series)"),
           panel = sprintf("%s (%s)", SubjectID, coalesce(relapse, "?")))
  for (a in tops) {
    d <- long |> filter(OlinkID == a)
    rl <- d |> distinct(panel, relapse_date) |> filter(!is.na(relapse_date))
    first_day <- d |> group_by(panel) |> summarise(d0 = min(date, na.rm = TRUE))
    d <- d |> left_join(first_day, by = "panel") |> mutate(day = as.numeric(date - d0))
    rl <- rl |> left_join(first_day, by = "panel") |> mutate(day = as.numeric(relapse_date - d0))
    p <- ggplot(d, aes(day, value, colour = series)) +
      geom_vline(data = rl, aes(xintercept = day), linetype = 2, colour = "grey40") +
      geom_line() + geom_point(aes(shape = point_state), size = 2) +
      scale_colour_manual(values = c(`dISF lesional site` = "firebrick", `dISF non-lesional site` = "steelblue",
                                     serum = "darkgoldenrod")) +
      scale_shape_manual(values = c(lesional = 17, `ex-lesional` = 2, `(other series)` = 16)) +
      facet_wrap(~panel) +
      labs(title = sprintf("%s - per-patient course (dashed line: relapse)", d$Assay[1]),
           x = "days since V1", y = "NPX", colour = NULL, shape = "lesional-site state") +
      theme(legend.position = "bottom")
    save_plot(p, cfg, "trajectories", "plots", sprintf("%s_%s.png", d$Assay[1], a), width = 11, height = 8)
  }
  msg("Trajectory plots for %d proteins in trajectories/plots/", length(tops))
}
