# 06 - Report: output/LEIP_galanin_report.pdf and output/answers.csv (all answers of steps 02-05).
# Every section is optional: if a step did not run, its page says so.

source("R/utils.R")
source("R/report.R")                  # PDF page helpers (page_text, page_table, page_plot, section)
source("R/leip.R")
cfg <- load_config()
d <- leip_load(cfg)
S <- d$S
dirs <- c("1_elisa_validation", "2_olink_galanin", "3_galanin_hdl", "4_clinical_screen")
answers <- map(dirs, \(x) out_csv(cfg, x, "answers.csv")) |> compact() |> bind_rows()
save_csv(answers, cfg, "answers.csv")
figs <- map(setNames(dirs, dirs), \(x) { f <- file.path(cfg$paths$output, x, "figures.rds"); if (file.exists(f)) readRDS(f) else list() })
show <- \(dir, name) { p <- figs[[dir]][[name]]; if (is.null(p)) stop("figure not available (step did not run?)"); page_plot(p) }
tab <- \(dir, file) { x <- out_csv(cfg, dir, file); if (is.null(x)) stop(file, " not found - did the step run?"); x }
n_clin <- sum(rowSums(!is.na(S[intersect(c("age", "BMI", "galanin_elisa", "C_HDL"), names(S))])) > 0)
q <- \(v, dg = 1) if (all(is.na(v))) "n/a" else
  sprintf(paste0("%.", dg, "f (IQR %.", dg, "f-%.", dg, "f)"), median(v, na.rm = TRUE), quantile(v, 0.25, na.rm = TRUE), quantile(v, 0.75, na.rm = TRUE))

out_file <- out_path(cfg, "LEIP_galanin_report.pdf")
pdf_open(out_file)
section("Overview", {
  sx <- col_or(S, "sex_male"); gdet <- d$det |> filter(OlinkID %in% d$gal)
  probs <- d$checks |> filter(!is.na(problem)) |> group_by(problem) |>
    summarise(which = if (n() <= 5) paste(SampleID, collapse = ", ") else sprintf("%d samples", n()), .groups = "drop")
  hdl_sex <- if (all(c("C_HDL", "sex_male") %in% names(S)))
    sprintf("HDL cholesterol: women %s, men %s mmol/l; galanin ELISA: women %s, men %s pg/mL.", q(S$C_HDL[S$sex_male %in% 0], 2),
            q(S$C_HDL[S$sex_male %in% 1], 2), q(d$elisa[S$sex_male %in% 0], 0), q(d$elisa[S$sex_male %in% 1], 0)) else NULL
  page_text("LEIP cohort: galanin - Olink validation of the ELISA, Olink galanin, galanin and HDL", c(
    "## Questions",
    "1. Does Olink confirm the galanin ELISA?  2. Olink galanin on its own, the ELISA ignored: 2a which clinical parameters and 2b which other Olink proteins go with it?  3. Is galanin related to HDL - could it bind to HDL particles?  Supplementary: all proteins vs all clinical parameters, with sanity checks.",
    "## Samples",
    sprintf("%d LEIP biobank sera in the Olink data after QC, %d with clinical data: %d women, %d men; age %s years; BMI %s.",
            nrow(S), n_clin, sum(sx == 0, na.rm = TRUE), sum(sx == 1, na.rm = TRUE), q(col_or(S, "age")), q(col_or(S, "BMI"))),
    hdl_sex,
    if (nrow(probs)) paste0("Sample checks: ", paste(sprintf("%s (%s)", probs$problem, probs$which), collapse = "; "), ". Details: 0_data/sample_checks.csv.")
    else "Sample checks: no problems.",
    "## Olink",
    sprintf("%d of %d proteins are measurable in LEIP serum (>= %.0f%% of the samples above LOD). Galanin: %s.", sum(d$det$measurable), nrow(d$det),
            100 * (cfg$min_detect_frac %||% 0.5), if (nrow(gdet)) sprintf("assay %s (%s), %.0f%% of the samples above LOD", gdet$Assay, gdet$OlinkID,
                                                                           100 * gdet$frac_above_lod) else "not measured"),
    "## Statistics",
    sprintf("Spearman rank correlation (robust, no transformation needed); adjusted = partial Spearman correlation (all variables ranked, covariates regressed out), usually for %s. Galanin questions are pre-specified: single-test p-values, with the FDR shown alongside. Proteome-wide lists: Benjamini-Hochberg FDR < %g.",
            covs_label(d$cov), cfg$fdr %||% 0.05),
    sprintf("Power: with n = %d a correlation needs about |rho| >= %.2f for p < 0.05, and about |rho| >= %.2f to pass the FDR over %d proteins. 'Not significant' does not mean 'no correlation'.",
            n_clin, r_crit(0.05, n_clin), r_crit(0.05 / max(ncol(d$M), 1), n_clin), ncol(d$M))),
    subtitle = sprintf("Olink Explore HT (%s); clinical data from the LEIP clinical file; generated %s", cfg$npx_column %||% "PCNormalizedNPX",
                       format(Sys.time(), "%Y-%m-%d %H:%M")))
})
section("Answers", {
  if (!nrow(answers)) stop("no answers - did steps 02-05 run?")
  items <- character()
  for (qq in unique(answers$question)) {
    a <- answers |> filter(question == qq)
    items <- c(items, paste("##", qq), sprintf("%s: %s.%s", a$item, a$verdict, if_else(is.na(a$evidence) | a$evidence == "", "", paste0(" ", a$evidence))))
  }
  page_text("Answers", items, subtitle = "derived automatically from the tests - the evidence is in the tables of each folder", size = 9)
})
section("1 Olink vs ELISA", {
  show("1_elisa_validation", "scatter")
  page_table("1  Agreement of Olink GAL with the ELISA", tab("1_elisa_validation", "agreement.csv") |>
               select(analysis, n, rho, ci_low, ci_high, p, npx_per_doubling),
             note = "rho = Spearman; npx_per_doubling = change of Olink NPX per doubling of the ELISA value (1 = same fold-change)")
  show("1_elisa_validation", "z_agreement")
  pf <- out_csv(cfg, "1_elisa_validation", "plate_effects.csv")
  if (!is.null(pf) && nrow(pf))
    page_table("1  Technical factors: do the values differ between plates?",
               pf |> mutate(test = str_remove(test, " \\(.*$")) |> separate_rows(medians, sep = "; ") |>
                 mutate(plate = str_remove(medians, ": [^:]*$"), median = as.numeric(str_extract(medians, "[^ ]+$"))) |>
                 group_by(test) |> mutate(across(c(n, groups, p), \(v) if_else(row_number() == 1, v, NA))) |> ungroup() |>
                 mutate(test = if_else(test == lag(test, default = ""), "", test)) |> select(test, n, groups, p, plate, median),
               note = "Kruskal-Wallis test; median per plate (ELISA in pg/mL, Olink GAL in NPX). ELISA values that differ between Olink plates: the samples were not placed on the plates at random")
  show("1_elisa_validation", "benchmark")
  idt <- out_csv(cfg, "1_elisa_validation", "sample_identity.csv")
  if (!is.null(idt) && nrow(idt))
    page_table("1  Per person: does the Olink sample fit the person's lab values?",
               idt |> mutate(top_galanin = rank(-abs(galanin_discordance), na.last = "keep") <= 3) |>
                 filter(row_number() <= 15 | coalesce(top_galanin, FALSE)) |>
                 select(SubjectID, plate, proteins, own_distance, own_rank, of, best_match, galanin_discordance, fit),
               note = "worst fits first, plus the 3 largest galanin disagreements; own_rank 1 = the own Olink sample fits best; possible swap: own sample not among the best 25%; all persons: 1_elisa_validation/sample_identity.csv")
  show("1_elisa_validation", "profiles")
  show("1_elisa_validation", "elisa_volcano")
  ec <- out_csv(cfg, "1_elisa_validation", "ELISA_vs_clinical.csv")
  if (!is.null(ec)) page_table("1  What does the ELISA follow instead? Clinical parameters with p < 0.05",
                               ec |> filter(p < 0.05) |> select(parameter = label, n, rho, p, fdr, rho_adj, p_adj),
                               note = sprintf("adjusted: partial Spearman for %s; about %.1f of %d parameters reach p < 0.05 by chance", covs_label(setdiff(d$cov, "plate")),
                                              0.05 * nrow(ec), nrow(ec)))
  eg <- out_csv(cfg, "1_elisa_validation", "ELISA_gene_sets.csv")
  if (!is.null(eg) && nrow(eg)) page_table("1  What does the ELISA follow instead? Gene sets (GSEA, top 15 by p)",
                                           eg |> slice_min(pval, n = 15, with_ties = FALSE) |> transmute(pathway = str_trunc(pathway, 60), size, NES, pval, padj),
                                           note = "GSEA on the ranking of all proteins by their correlation with the ELISA")
})
section("2 Olink galanin on its own (ELISA ignored)", {
  g2 <- "2_olink_galanin"
  show_if <- \(n) if (!is.null(figs[[g2]][[n]])) show(g2, n)
  show_if("clinical")
  cl <- tab(g2, "olink_galanin_vs_clinical.csv")
  page_table("2a  Olink galanin vs clinical parameters: p < 0.05 unadjusted or adjusted",
             cl |> filter(coalesce(p < 0.05 | p_adj < 0.05, FALSE)) |>
               select(parameter = label, n, rho, p, fdr, rho_adj, p_adj, rho_women, rho_men, loo_min, loo_max, robust),
             note = sprintf("adjusted: partial Spearman for %s; loo = range of rho leaving out one person; robust: see the answers. All: %s/olink_galanin_vs_clinical.csv",
                            covs_label(d$cov), g2))
  show_if("clinical_scatter")
  tf <- out_csv(cfg, g2, "technical_factors.csv")
  if (!is.null(tf)) page_table("2a  Technical factors: plates and flagged samples",
                               tf |> separate_rows(detail, sep = "; ") |>
                                 mutate(across(c(n_flagged, p), \(v) if_else(factor == lag(factor, default = ""), NA, v)),
                                        factor = if_else(factor == lag(factor, default = ""), "", factor)),
                               note = "Olink plate: median GAL NPX per plate, Kruskal-Wallis p; flagged samples: GAL z-score (p only with >= 3 flagged)")
  show_if("volcano")
  page_table("2b  Proteins most strongly correlated with Olink galanin", tab(g2, "olink_galanin_vs_proteins.csv") |> head(30) |>
               select(protein, n, rho, p, fdr, rho_adj, p_adj, frac_above_lod),
             note = sprintf("adjusted: partial Spearman for %s; all proteins: %s/olink_galanin_vs_proteins.csv", covs_label(d$cov), g2))
  show_if("top_proteins")
  show_if("heatmap")
  show_if("neuroendocrine")
  gs <- out_csv(cfg, g2, "gene_sets.csv")
  if (!is.null(gs) && nrow(gs)) page_table("2b  Gene sets (GSEA on the ranking by correlation with Olink galanin)",
                                           gs |> slice_min(pval, n = 20, with_ties = FALSE) |> transmute(pathway = str_trunc(pathway, 60), size, NES, pval, padj))
  show_if("proteome_axes")
  pa <- out_csv(cfg, g2, "proteome_axes.csv")
  if (!is.null(pa) && nrow(pa)) page_table("2b  Main axes of the serum proteome (principal components)",
                                           pa |> transmute(component, variance_pct, rho_GAL, p_GAL, plate_p, strongest_clinical, top_proteins = str_trunc(top_proteins, 38)),
                                           note = "top_proteins: the proteins most correlated with each component (all in proteome_axes.csv)")
})
section("3 Galanin and HDL", {
  show("3_galanin_hdl", "hdl_by_sex")
  show("3_galanin_hdl", "lipids")
  page_table("3  Galanin vs HDL cholesterol and ApoA-I", tab("3_galanin_hdl", "galanin_vs_lipids.csv") |>
               filter(lipid %in% c(cfg$hdl$hdl %||% "C_HDL", cfg$hdl$apoa1 %||% "c_apo")) |>
               select(measure, lipid = label, analysis, n, rho, ci_low, ci_high, p), rows_per_page = 24)
  sp <- out_csv(cfg, "3_galanin_hdl", "HDL_specificity.csv")
  if (!is.null(sp) && nrow(sp)) page_table("3  Is galanin closer to HDL than to the other lipids? (bootstrap)",
                                           sp |> select(measure, comparison, analysis, n, difference, ci_low, ci_high, p_boot),
                                           note = "difference > 0 with a CI above 0: galanin follows the first lipid more closely than the second")
  if (!is.null(figs[["3_galanin_hdl"]]$hdl_proteins)) show("3_galanin_hdl", "hdl_proteins")
  show("3_galanin_hdl", "hdl_volcano")
  if (!is.null(figs[["3_galanin_hdl"]]$discordance)) show("3_galanin_hdl", "discordance")
  if (!is.null(figs[["3_galanin_hdl"]]$agreement_by_hdl)) show("3_galanin_hdl", "agreement_by_hdl")
  ds <- out_csv(cfg, "3_galanin_hdl", "discordance_vs_HDL.csv")
  if (!is.null(ds) && nrow(ds)) page_table("3  ELISA-Olink disagreement vs HDL", ds |> select(with, discordance, analysis, n, rho, ci_low, ci_high, p),
                                           note = "signed = z(ELISA) - z(Olink GAL), pre-specified; absolute = its size in either direction, added after the first results (exploratory)")
  it <- out_csv(cfg, "3_galanin_hdl", "agreement_by_HDL.csv")
  if (!is.null(it) && nrow(it)) page_table("3  Does the ELISA-Olink agreement weaken as HDL rises? (exploratory)",
                                           it |> select(with, analysis, n, slope_low, slope_high, interaction, ci_low, ci_high, p),
                                           note = "Olink GAL ~ ELISA x HDL on standardised ranks (+ sex, Olink plate); slope = agreement at HDL -1 SD (low) and +1 SD (high); interaction < 0: weaker where HDL is high")
})
section("Supplementary: proteins vs clinical parameters", {
  if (!is.null(figs[["4_clinical_screen"]]$nsig)) show("4_clinical_screen", "nsig")
  if (!is.null(figs[["4_clinical_screen"]]$heatmap)) show("4_clinical_screen", "heatmap")
  san <- out_csv(cfg, "4_clinical_screen", "sanity_checks.csv")
  if (!is.null(san)) page_table("Sanity checks: associations known from population studies",
                                san |> select(protein, parameter = label, expected, n, rho, p, rho_adj, p_adj, unadjusted, status),
                                note = sprintf("status: adjusted for %s (a covariate is left out when it is the parameter); recovered = p < 0.05 in the expected direction", covs_label(d$cov)))
})
section("Methods", page_text("Methods, caveats and files", c(
  "## Data",
  "LEIP biobank sera (population controls) of the O-MicroAD Olink Explore HT run. The LEIP samples are the rows of the clinical file (Olink_SampleID). Olink values: PC-normalised NPX from the delivered parquet file, LOD from the Olink fixed LOD file (per sample for count-based assays), samples failing Olink QC left out. Clinical data: sheet Key_parameters plus all further variables of sheet All_SORB_parameters; identifiers, log copies, constants, too-small groups and duplicate variables are not tested (0_data/parameters.csv).",
  "## Statistics",
  "Spearman correlation, p-value from the t approximation, 95% CI by Fisher z (Bonett-Wright). Partial Spearman: all variables ranked, covariates regressed out. Olink galanin on its own (2): the ELISA is not used; clinical parameters also within each sex and leaving out one person at a time; proteome axes = principal components of all measurable proteins (scaled; missing values replaced by the protein median). Differences between two correlations with galanin (e.g. HDL vs LDL): bootstrap over persons. Sample identity per person: distance between the lab values and the Olink values of the benchmark proteins (rank-based normal scores, weighted by 1 / (2 (1 - rho))); the own sample should be among the best matches. Olink vs ELISA: also adjusted for sex, within each sex, within plates, above LOD, without flagged samples, leave-one-out, tertile agreement (weighted kappa), NPX per doubling of the ELISA; specificity = rank of Olink GAL among all proteins correlated with the ELISA, and agreement of the two protein-correlation profiles (permutation test). Gene sets: fgsea on the rho ranking.",
  "Added after the first results (exploratory): does the ELISA-Olink agreement weaken as HDL rises - linear model of Olink GAL on ELISA x HDL on standardised ranks (+ sex, Olink plate) - and the size of the disagreement |z(ELISA) - z(Olink GAL)| vs HDL. Sanity checks are judged adjusted for age, sex and Olink plate.",
  "## Caveats",
  "n = 34: exploratory; weak correlations are missed and single p < 0.05 results can be chance (with ~3000 proteins, |rho| of about 0.55 occurs by chance alone). Correlation is not binding: the HDL analysis shows whether the data are compatible with binding; experiments are needed (see the last answer of 3).",
  "Sex matters: women have higher HDL, and galanin may differ by sex; results adjusted for sex and within each sex are shown.",
  "NPX is relative: Olink and the ELISA can only agree in ranking. The Olink assay targets the galanin precursor (UniProt P22466); a galanin ELISA may detect the mature peptide or other fragments. Biobank serum was not collected for peptide measurements (proteolysis).",
  "## Files (output/)",
  "answers.csv | 0_data/ (samples, checks, detection, parameters) | 1_elisa_validation/ | 2_olink_galanin/ | 3_galanin_hdl/ | 4_clinical_screen/ - each with an .xlsx of all its tables and the figures")))
invisible(grDevices::dev.off())
msg("Report: %s", out_file)
