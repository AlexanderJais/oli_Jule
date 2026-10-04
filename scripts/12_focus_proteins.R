# 12 - Dedicated analysis of pre-specified proteins (config: focus_proteins; default CD137 = TNFRSF9)
# For each focus protein, independent of the detection filter of step 02:
#   1. detection per matrix and group (is it measurable, and where?)
#   2. pre-specified tests: the same models as steps 04/05 (and the xL - NL relapse model), fitted for
#      this one protein. Because the protein was chosen in advance, the unadjusted p-value is the
#      primary test; the proteome-wide FDR from steps 04/05 is shown alongside.
#   3. ISF-serum correlation on matched visits; optional severity; LEIP clinical associations
#   4. every result for the protein from the proteome-wide steps, collected in one table
#   5. figures: skin states, serum groups, per-patient course, ISF vs serum
# Out: output/focus/<protein>/<protein>_report.xlsx + figures

source("R/utils.R")
source("R/models.R")
source("R/design.R")
cfg   <- load_config()
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
clear_outputs(cfg, "focus")
focus <- unlist(cfg$focus_proteins %||% "TNFRSF9")
aliases <- c(TNFRSF9 = "CD137 / 4-1BB", TNFSF9 = "CD137L / 4-1BBL", KITLG = "SCF / KIT ligand",
             CPA4 = "carboxypeptidase A4", FCER1A = "FceRI alpha", TPSAB1 = "tryptase alpha/beta-1",
             MS4A2 = "FceRI beta", TPSD1 = "tryptase delta-1", PNOC = "prepronociceptin", POSTN = "periostin")
overview <- list(); found <- character()

for (fp in focus) {
  # case-insensitive; also finds the protein inside combined assay names (e.g. "TPSAB1_TPSB2")
  FP <- toupper(fp)
  pd <- clean |> filter(toupper(Assay) == FP | OlinkID == fp | UniProt == fp |
                          str_detect(toupper(Assay), paste0("(^|_)", FP, "($|_)")))
  if (!nrow(pd)) { msg("%s: NOT in the NPX data (checked Assay, OlinkID, UniProt) - skipped", fp); next }
  if (n_distinct(pd$OlinkID) > 1) {
    msg("%s: several assays match (%s) - using %s", fp, paste(unique(pd$Assay), collapse = ", "), pd$Assay[1])
    pd <- pd |> filter(OlinkID == OlinkID[1])
  }
  nm <- pd$Assay[1]; oid <- pd$OlinkID[1]
  found <- c(found, nm)
  lab <- if (nm %in% names(aliases)) sprintf("%s (%s)", nm, aliases[[nm]]) else nm
  fpath <- function(...) file.path("focus", nm, ...)
  msg("==== Focus protein %s (%s) ====", nm, oid)
  val <- pd |> select(SampleID, value, below_lod, LOD)

  # ---- 1. detection -----------------------------------------------------------------------------
  det <- pd |>
    mutate(where = case_when(matrix == "ISF" ~ paste("dISF", det_group),
                             TRUE ~ paste("serum", cohort, group))) |>
    group_by(matrix, where) |>
    summarise(n = n(), pct_above_LOD = 100 * mean(!below_lod, na.rm = TRUE),
              median_NPX = median(value, na.rm = TRUE), median_LOD = median(LOD, na.rm = TRUE), .groups = "drop")
  keep <- read_csv(file.path(cfg$paths$output, "qc", "assay_detection.csv"), show_col_types = FALSE) |>
    filter(OlinkID == oid) |> select(matrix, kept_in_proteome_analysis = keep, frac_detected_all)
  print(det, n = Inf)
  msg("%s passes the detection filter: %s", nm,
      paste(keep$matrix, if_else(keep$kept_in_proteome_analysis, "yes", "NO"), collapse = ", "))

  # ---- 2. pre-specified tests -----------------------------------------------------------------------
  isf_info   <- isf_design(meta) |> inner_join(val, by = "SampleID")
  serum_info <- serum_design(meta) |> inner_join(val, by = "SampleID")
  tests <- prespecified_tests(isf_info, serum_info, cfg$stats$min_group_n)
  # proteome-wide FDR for the same contrasts, where the protein was in steps 04/05
  pw <- map(c("ISF_results.csv", "Serum_results.csv", "ISF_relapse_delta_results.csv"), \(f) {
    p <- file.path(cfg$paths$output, "models", f)
    if (file.exists(p)) read_csv(p, show_col_types = FALSE) |> filter(OlinkID == oid) |>
      select(model, contrast, proteome_wide_FDR = adj.P.Val) else NULL
  }) |> bind_rows()
  if (nrow(tests)) {
    tests <- tests |> left_join(pw, by = c("model", "contrast")) |>
      mutate(significant_p05 = p < 0.05) |> relocate(matrix, .before = model)
    print(tests |> select(matrix, model, contrast, estimate, p, proteome_wide_FDR) |> as.data.frame(), digits = 3)
  }

  # ---- 3. ISF-serum correlation, severity, LEIP ----------------------------------------------------
  micro <- meta |> filter(cohort == "MicroAD") |> inner_join(val, by = "SampleID")
  pr <- micro |> filter(matrix == "ISF") |>
    transmute(SubjectID, visit, group, site = if_else(site == "L", "lesional site", "non-lesional / healthy skin"),
              isf = value) |>
    inner_join(micro |> filter(matrix == "Serum") |> select(SubjectID, visit, serum = value), by = c("SubjectID", "visit"))
  corr <- pr |> group_by(site) |> group_modify(\(d, k) {
    multi <- d |> filter(SubjectID %in% names(which(table(SubjectID) >= 2)))
    w <- if (n_distinct(multi$SubjectID) >= 2) suppressWarnings(rmcorr::rmcorr(participant = SubjectID, measure1 = isf,
                                                                               measure2 = serum, dataset = as.data.frame(multi))) else NULL
    m <- d |> group_by(SubjectID) |> summarise(isf = mean(isf), serum = mean(serum))
    b <- suppressWarnings(cor.test(m$isf, m$serum, method = "spearman", exact = FALSE))
    tibble(n_pairs = nrow(d), n_subjects = n_distinct(d$SubjectID),
           r_within = w$r %||% NA_real_, p_within = w$p %||% NA_real_,
           rho_between = unname(b$estimate), p_between = b$p.value)
  }) |> ungroup()

  sev_tab <- NULL
  sev <- read_severity(cfg$paths$severity)
  if (!is.null(sev)) {
    scores <- severity_scores(sev)
    si <- micro |> filter(group == "AD") |> left_join(sev, by = c("SubjectID", "visit")) |>
      mutate(series = case_when(matrix == "Serum" ~ "serum", site == "L" ~ "dISF lesional site", TRUE ~ "dISF non-lesional site"))
    sev_tab <- map(scores, \(sc) map(unique(si$series), \(se) {
      test_single(si |> filter(series == se), as.formula(sprintf("~ %s + plate + (1 | SubjectID)", sc)),
                  c(per_unit = sc), paste(se, sc, sep = ": "), cfg$stats$min_group_n)
    }) |> bind_rows()) |> bind_rows()
  }

  leip_tab <- NULL
  leip <- meta |> filter(cohort == "LEIP") |> inner_join(val, by = "SampleID")
  params <- intersect(c("age", "BMI", "WHR", "c_fett", "HOMA_IR", "c_CRP", "C_CHOL", "C_HDL", "C_LDL",
                        "C_TRIGLY", "MDRD_kurz"), names(leip))
  if (nrow(leip) >= 10 && length(params)) {
    leip_tab <- map(params, \(pm) {
      ok <- !is.na(leip[[pm]]) & !is.na(leip$value)
      if (sum(ok) < 10) return(NULL)
      ct <- suppressWarnings(cor.test(leip$value[ok], leip[[pm]][ok], method = "spearman", exact = FALSE))
      tibble(parameter = pm, n = sum(ok), rho = unname(ct$estimate), p = ct$p.value)
    }) |> bind_rows()
    if ("sex" %in% names(leip) && all(c("M", "F") %in% leip$sex)) {
      w <- suppressWarnings(wilcox.test(value ~ sex, data = leip |> filter(sex %in% c("M", "F")), exact = FALSE))
      leip_tab <- bind_rows(leip_tab, tibble(parameter = "sex (M - F, median difference)", n = sum(leip$sex %in% c("M", "F")),
                                             rho = median(leip$value[leip$sex == "M"], na.rm = TRUE) -
                                               median(leip$value[leip$sex == "F"], na.rm = TRUE), p = w$p.value))
    }
  }

  # ---- 4. everything the proteome-wide steps say about this protein -------------------------------------
  collect <- \(path, what) {
    p <- file.path(cfg$paths$output, path)
    if (!file.exists(p)) return(NULL)
    read_csv(p, show_col_types = FALSE) |> filter(OlinkID == oid) |> mutate(across(everything(), as.character)) |>
      mutate(source = what, .before = 1)
  }
  pipeline <- list(
    isf_models = collect("models/ISF_results.csv", "04 ISF models"),
    serum_models = collect("models/Serum_results.csv", "05 serum models"),
    isf_serum = collect("isf_serum/isf_serum_correlation.csv", "06 ISF-serum correlation"),
    matrix_enrich = collect("matrix_comparison/relative_enrichment.csv", "10 relative dISF/serum level"),
    trajectories = collect("trajectories/trajectory_results.csv", "11 trajectories"),
    leip = collect("leip_reference/leip_reference_summary.csv", "08 LEIP reference")
  ) |> compact()

  # ---- 5. figures -----------------------------------------------------------------------------------------
  lod_isf <- median(pd$LOD[pd$matrix == "ISF"], na.rm = TRUE)
  lod_ser <- median(pd$LOD[pd$matrix == "Serum"], na.rm = TRUE)
  cond_lab <- c(HC = "healthy skin", AD_NL = "AD non-lesional", AD_xL = "AD ex-lesional", AD_L = "AD lesional",
                CPUO_NL = "CPUO non-lesional", CPUO_L = "CPUO lesional")
  pi <- isf_info |> filter(!is.na(cond)) |> mutate(cond = factor(cond_lab[cond], levels = cond_lab))
  p1 <- ggplot(pi, aes(cond, value)) +
    geom_hline(yintercept = lod_isf, linetype = 2, colour = "grey50") +
    geom_boxplot(outlier.shape = NA, fill = "grey92") + geom_jitter(aes(colour = below_lod), width = 0.15, size = 1.2) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "grey65"), labels = c("above LOD", "below LOD"), name = NULL) +
    labs(title = sprintf("%s in dISF by skin state (dashed: median LOD)", lab), x = NULL, y = "NPX") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p1, cfg, fpath(sprintf("%s_dISF_skin_states.png", nm)), width = 8, height = 5)

  v1 <- isf_info |> filter(visit == "V1", group %in% c("AD", "CPUO"), site %in% c("L", "NL"))
  if (nrow(v1)) {
    p2 <- ggplot(v1, aes(site, value, group = SubjectID, colour = group)) + geom_line(alpha = 0.6) + geom_point() +
      labs(title = sprintf("%s at V1: lesional vs non-lesional skin, per patient", lab), x = "skin site", y = "NPX (dISF)")
    save_plot(p2, cfg, fpath(sprintf("%s_dISF_V1_paired.png", nm)), width = 6, height = 5)
  }

  ps <- serum_info |> filter(cross_sectional | cohort == "MicroAD") |>
    mutate(grp = case_when(status == "HC" ~ "healthy (in-study)", status == "Biobank" ~ "LEIP biobank",
                           TRUE ~ paste("AD", cohort)),
           relapse_lab = coalesce(relapse, "")) |>
    filter(cross_sectional)
  p3 <- ggplot(ps, aes(grp, value)) +
    geom_hline(yintercept = lod_ser, linetype = 2, colour = "grey50") +
    geom_boxplot(outlier.shape = NA, fill = "grey92") + geom_jitter(aes(colour = relapse2), width = 0.15, size = 1.2) +
    labs(title = sprintf("%s in serum (one sample per person; MicroAD = V1)", lab), x = NULL, y = "NPX", colour = "relapse") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p3, cfg, fpath(sprintf("%s_serum_groups.png", nm)), width = 8, height = 5)

  tr <- micro |> filter(group == "AD") |>
    mutate(series = case_when(matrix == "Serum" ~ "serum", site == "L" ~ "dISF lesional site", TRUE ~ "dISF non-lesional site"),
           panel = sprintf("%s (%s)", SubjectID, coalesce(relapse, "?"))) |>
    group_by(SubjectID) |> mutate(day = as.numeric(date - min(date, na.rm = TRUE))) |> ungroup()
  rl <- tr |> filter(!is.na(relapse_visit), visit_num == relapse_visit) |> distinct(panel, day)
  p4 <- ggplot(tr, aes(day, value, colour = series)) +
    geom_vline(data = rl, aes(xintercept = day), linetype = 2, colour = "grey40") +
    geom_line() + geom_point(aes(shape = if_else(series == "dISF lesional site", state, "(other)")), size = 2) +
    scale_colour_manual(values = c(`dISF lesional site` = "firebrick", `dISF non-lesional site` = "steelblue", serum = "darkgoldenrod")) +
    scale_shape_manual(values = c(lesional = 17, `ex-lesional` = 2, `(other)` = 16)) +
    facet_wrap(~panel) + labs(title = sprintf("%s per patient over time (dashed: relapse)", lab), x = "days since V1",
                              y = "NPX", colour = NULL, shape = "lesional-site state") +
    theme(legend.position = "bottom")
  save_plot(p4, cfg, fpath(sprintf("%s_trajectories.png", nm)), width = 11, height = 8)

  if (nrow(pr)) {
    p5 <- ggplot(pr, aes(isf, serum, colour = SubjectID)) + geom_point() +
      geom_line(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, alpha = 0.6) +
      facet_wrap(~site) + guides(colour = "none") +
      labs(title = sprintf("%s: dISF vs serum at matched visits (lines: within-patient fits)", lab), x = "dISF NPX", y = "serum NPX")
    save_plot(p5, cfg, fpath(sprintf("%s_dISF_vs_serum.png", nm)), width = 10, height = 5)
  }

  sheets <- c(list(detection = det, detection_filter = keep, prespecified_tests = tests,
                   isf_serum_correlation = corr),
              if (!is.null(sev_tab) && nrow(sev_tab)) list(severity = sev_tab),
              if (!is.null(leip_tab) && nrow(leip_tab)) list(leip_clinical = leip_tab),
              setNames(pipeline, paste0("pipeline_", names(pipeline))),
              list(sample_values = pd |> select(SampleID, SubjectID, matrix, cohort, group, visit, site, state,
                                                plate, value, LOD, below_lod)))
  writexl::write_xlsx(sheets, out_path(cfg, fpath(sprintf("%s_report.xlsx", nm))))
  msg("%s: report and figures in %s", nm, file.path(cfg$paths$output, "focus", nm))
  overview[[nm]] <- list(
    tests = if (nrow(tests)) tests |> mutate(protein = nm, label = lab, .before = 1) else NULL,
    det = det |> mutate(protein = nm, .before = 1),
    values = isf_info |> filter(!is.na(cond)) |> transmute(protein = lab, cond, value, below_lod))
}

# ---- overview across all focus proteins ----------------------------------------------------------------------
if (length(overview)) {
  ov_tests <- map(overview, "tests") |> bind_rows()
  key <- c("states_all_visits AD_L_vs_NL", "states_all_visits AD_xL_vs_NL", "states_all_visits AD_L_vs_HC",
           "states_all_visits AD_NL_vs_HC", "AD_vs_HC_in_study AD_vs_HC", "MicroAD_AD_vs_HC AD_vs_HC", "MicroAD_active_vs_cleared active_vs_cleared",
           "relapse_delta_xL_minus_NL relapse_vs_non", "MicroAD_relapse relapse_vs_non")
  ov <- ov_tests |> filter(paste(model, contrast) %in% key) |>
    transmute(protein, label, matrix, comparison = paste0(matrix, ": ", model, " ", contrast), estimate, ci_low, ci_high, p,
              proteome_wide_FDR) |>
    mutate(comparison = factor(comparison, levels = unique(comparison[order(matrix != "dISF")])))
  save_csv(ov, cfg, "focus", "focus_overview.csv")
  ov_det <- map(overview, "det") |> bind_rows()
  writexl::write_xlsx(list(overview = ov, all_tests = ov_tests, detection = ov_det), out_path(cfg, "focus", "focus_overview.xlsx"))
  vals <- map(overview, "values") |> bind_rows() |>
    mutate(cond = factor(cond_lab[cond], levels = cond_lab))
  p <- ggplot(vals, aes(cond, value)) + geom_boxplot(outlier.shape = NA, fill = "grey92") +
    geom_jitter(aes(colour = below_lod), width = 0.15, size = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "grey65"), labels = c(`FALSE` = "above LOD", `TRUE` = "below LOD"), name = NULL) +
    facet_wrap(~protein, scales = "free_y") +
    labs(title = "Focus proteins in dISF by skin state", x = NULL, y = "NPX") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
  save_plot(p, cfg, "focus", "focus_skin_states.png", width = 13, height = 9)
  hm <- ov |> mutate(z = sign(estimate) * pmin(-log10(p), 10))
  p <- ggplot(hm, aes(comparison, label, fill = estimate)) + geom_tile() +
    geom_text(aes(label = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")), size = 4) +
    scale_fill_gradient2(low = "steelblue", high = "firebrick") +
    labs(title = "Focus proteins: effect (log2) per comparison (* p < 0.05, ** < 0.01, *** < 0.001; single-protein tests)",
         x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 35, hjust = 1))
  save_plot(p, cfg, "focus", "focus_overview_heatmap.png", width = 12, height = 2.5 + 0.4 * n_distinct(hm$label))
}
missing <- setdiff(toupper(focus), toupper(found))
if (length(missing)) msg("Focus proteins not in the data: %s", paste(missing, collapse = ", "))
