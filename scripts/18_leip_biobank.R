# 18 - LEIP biobank serum only: Olink proteins vs clinical parameters, with a focus on galanin
# Uses only the LEIP biobank samples (population controls) and their clinical data from the SORB
# database (paths$leip_clinical: sheet Key_parameters and, if present, All_SORB_parameters).
#   1. proteins measurable in LEIP serum (>= qc: min_detect_frac of the LEIP samples above LOD)
#   2. every protein vs every clinical parameter: Spearman correlation (primary) and partial
#      Spearman adjusted for age, sex and Olink plate; BH-FDR over the proteins of each parameter
#   3. sanity checks: associations known from population proteomics (e.g. leptin - BMI)
#   4. galanin (Olink assay GAL, pre-specified): is it measurable; do the Olink values agree with
#      the galanin ELISA (benchmarked against the other proteins measured by both the lab and
#      Olink); which clinical parameters and which proteins does galanin correlate with
# Out: output/leip_biobank/ (LEIP_biobank_summary.pdf, leip_biobank.xlsx, answers.csv, figures)
#      output/leip_biobank/galanin/ (galanin.xlsx, figures)

source("R/utils.R")
source("R/leip.R")
source("R/report.R")
cfg     <- load_config()
lc      <- cfg$leip_biobank %||% list()
fdr_cut <- cfg$stats$fdr
meta    <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean   <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
if (is.null(cfg$paths$leip_clinical) || !file.exists(cfg$paths$leip_clinical))
  stop("Step 18 needs the LEIP clinical file (paths: leip_clinical in config.yml): ", cfg$paths$leip_clinical %||% "not set")
clear_outputs(cfg, "leip_biobank")

gal_name  <- lc$galanin_assay %||% "GAL"
adjust    <- unlist(lc$adjust_for %||% c("age", "sex_male", "plate"))
min_n     <- lc$min_n %||% 20
min_group <- lc$min_group %||% 5
fmt_p   <- \(p) ifelse(is.na(p), "n/a", sprintf("%.2g", p))
col_or  <- \(d, nm, default = NA) if (nm %in% names(d)) d[[nm]] else rep(default, nrow(d))
ensure  <- \(d, cols) { for (c in setdiff(cols, names(d))) d[[c]] <- rep(NA_real_, nrow(d)); d }
top_str <- \(name, rho, p, k = 5) {
  o <- head(order(p), k); o <- o[!is.na(p[o])]
  if (!length(o)) "none" else paste(sprintf("%s (%+.2f)", name[o], rho[o]), collapse = ", ")
}
find_assay <- \(det, a) det$OlinkID[str_detect(toupper(det$Assay), paste0("(^|_)", toupper(a), "($|_)"))]
r_crit  <- \(p, n) { t <- qt(1 - p / 2, n - 2); t / sqrt(t^2 + n - 2) }   # |rho| needed for p at n samples

# ---- 1. samples: Olink data and clinical file --------------------------------------------------------------
cl   <- read_leip_clinical(cfg$paths$leip_clinical, lc$galanin_lab %||% "Galanin [pg/mL]")
clin <- cl$data
lp   <- clean |> filter(cohort == "LEIP", matrix == "Serum")
if (!nrow(lp)) stop("No LEIP serum samples in the cleaned data of step 02.")
leip_meta <- meta |> filter(cohort == "LEIP", matrix == "Serum")

# each clinical row must belong to the same person and Olink plate as in the manifest
norm_plate <- \(x) str_remove_all(str_to_lower(as.character(x)), "\\s")
key_vals <- setdiff(intersect(cl$key, names(clin)), leip_technical)
clin_ids <- tibble(SampleID = clin$SampleID, SubjectID_clinical = as.character(col_or(clin, "SubjectID")),
                   plate_clinical = as.character(col_or(clin, "Olink_plate")),
                   n_clinical_values = if (length(key_vals)) rowSums(!is.na(clin[key_vals])) else 0)
id_check <- full_join(leip_meta |> select(SampleID, SubjectID_manifest = SubjectID, plate_manifest = plate),
                      clin_ids, by = "SampleID") |>
  mutate(in_olink_data = SampleID %in% lp$SampleID,
         problem = case_when(
           is.na(SubjectID_manifest) ~ "in the clinical file, but not a LEIP serum sample in the manifest",
           is.na(n_clinical_values) ~ "no row in the clinical file",
           !is.na(SubjectID_clinical) & SubjectID_clinical != SubjectID_manifest ~ "SubjectID differs between manifest and clinical file",
           !is.na(plate_clinical) & norm_plate(plate_clinical) != norm_plate(plate_manifest) ~ "Olink plate differs between manifest and clinical file",
           n_clinical_values == 0 ~ "no clinical values",
           !in_olink_data ~ "not in the Olink data after QC",
           TRUE ~ NA_character_)) |>
  arrange(SampleID)
save_csv(id_check, cfg, "leip_biobank", "sample_check.csv")
if (any(!is.na(id_check$problem))) {
  msg("LEIP sample checks (leip_biobank/sample_check.csv):")
  print(id_check |> filter(!is.na(problem)) |> select(SampleID, SubjectID_manifest, SubjectID_clinical, problem), n = Inf)
}

W <- lp |> select(SampleID, OlinkID, value) |> pivot_wider(names_from = OlinkID, values_from = value) |> arrange(SampleID)
Yall <- as.matrix(W[, -1]); rownames(Yall) <- W$SampleID
S <- W |> select(SampleID) |>
  left_join(leip_meta |> select(SampleID, SubjectID, plate), by = "SampleID") |>
  left_join(clin |> select(-any_of("SubjectID")), by = "SampleID")
sq <- out_csv(cfg, "qc", "sample_qc.csv")

# proteins measurable in LEIP serum: the filter of step 02, applied to the LEIP samples alone
det <- lp |>
  group_by(OlinkID, Assay, across(any_of("UniProt"))) |>
  summarise(n_samples = sum(!is.na(value)),
            frac_above_lod = if (all(is.na(below_lod))) NA_real_ else mean(!below_lod, na.rm = TRUE),
            median_npx = median(value, na.rm = TRUE), median_lod = median(LOD, na.rm = TRUE), .groups = "drop") |>
  mutate(analysed = n_samples >= min_n & coalesce(frac_above_lod >= cfg$qc$min_detect_frac, TRUE),
         protein = if_else(duplicated(Assay) | duplicated(Assay, fromLast = TRUE), paste0(Assay, " (", OlinkID, ")"), Assay)) |>
  arrange(Assay)
save_csv(det, cfg, "leip_biobank", "detection_in_LEIP.csv")
Y <- Yall[, det$OlinkID[det$analysed], drop = FALSE]

# ---- 2. clinical parameters and covariates -----------------------------------------------------------------
pinfo <- select_parameters(clin |> filter(SampleID %in% S$SampleID), cl$key, unlist(lc$exclude_parameters), min_n, min_group)
save_csv(pinfo, cfg, "leip_biobank", "parameters.csv")
params <- pinfo$parameter[pinfo$analysed]
if (!length(params)) stop("No clinical parameter passes the filters - see leip_biobank/parameters.csv")
ptype <- setNames(pinfo$type, pinfo$parameter)

adj <- intersect(adjust, names(S))
if (length(setdiff(adjust, adj)))
  msg("WARNING: covariate(s) %s not found - the adjusted analysis goes without", paste(setdiff(adjust, adj), collapse = ", "))
cov_for <- function(pm = "", drop = character()) {
  z <- S[setdiff(adj, c(pm, drop))]
  if (ncol(z)) z else NULL
}
adj_lab <- \(v) if (length(v)) paste(recode(v, sex_male = "sex", plate = "Olink plate"), collapse = ", ") else "nothing"
adj_txt <- adj_lab(adj)
n_clin  <- sum(id_check$in_olink_data & coalesce(id_check$n_clinical_values, 0) > 0)
msg("LEIP: %d samples in the Olink data (%d with clinical data); %d of %d proteins measurable; %d clinical parameters (%d key)",
    nrow(S), n_clin, ncol(Y), nrow(det), length(params), sum(pinfo$analysed & pinfo$key))

# ---- 3. every protein vs every clinical parameter ------------------------------------------------------------
assoc <- map(params, \(pm) {
  x <- S[[pm]]
  u <- spearman_vs(x, Y, min_n = min_n)
  a <- spearman_vs(x, Y, cov_for(pm), min_n = min_n)
  tibble(parameter = pm, OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
         median_diff = if (ptype[[pm]] == "binary") unname(median_diff(x, Y)) else NA_real_,
         n_adj = a$n, rho_adj = a$rho, p_adj = a$p)
}) |> bind_rows() |>
  group_by(parameter) |>
  mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH")) |>
  ungroup() |>
  mutate(significant = coalesce(fdr < fdr_cut, FALSE), significant_adj = coalesce(fdr_adj < fdr_cut, FALSE)) |>
  left_join(det |> select(OlinkID, Assay, protein), by = "OlinkID") |>
  left_join(pinfo |> select(parameter, label, type, key), by = "parameter") |>
  relocate(parameter, label, type, key, OlinkID, Assay, protein) |>
  arrange(match(parameter, params), p)

nsig <- assoc |>
  group_by(parameter, label, type, key) |>
  summarise(n_samples = max(n, na.rm = TRUE), proteins_tested = sum(!is.na(p)),
            # positive/negative first: summarise() sees a summary as soon as it is created
            positive = sum(significant & rho > 0), negative = sum(significant & rho < 0), significant = sum(significant),
            positive_adj = sum(significant_adj & rho_adj > 0), negative_adj = sum(significant_adj & rho_adj < 0),
            significant_adjusted = sum(significant_adj),
            p_below_05 = sum(p < 0.05, na.rm = TRUE), expected_by_chance = round(0.05 * sum(!is.na(p))),
            top_proteins = top_str(protein, rho, p), .groups = "drop") |>
  relocate(significant, .before = positive) |> relocate(significant_adjusted, .before = positive_adj) |>
  arrange(desc(significant), desc(significant_adjusted))
sig_tab <- assoc |> filter(significant | significant_adj)
top10   <- assoc |> group_by(parameter) |> slice_min(p, n = 10, with_ties = FALSE) |> ungroup()
save_csv(assoc, cfg, "leip_biobank", "associations_all.csv.gz")
save_csv(sig_tab, cfg, "leip_biobank", "associations_significant.csv")
save_csv(nsig, cfg, "leip_biobank", "n_significant_per_parameter.csv")
save_csv(assoc |> filter(key) |> select(OlinkID, protein, label, rho) |> pivot_wider(names_from = label, values_from = rho),
         cfg, "leip_biobank", "rho_key_parameters_wide.csv")
msg("%d protein x parameter tests; %d parameters with proteins at FDR < %g (adjusted: %d)",
    nrow(assoc), sum(nsig$significant > 0), fdr_cut, sum(nsig$significant_adjusted > 0))

# ---- 4. sanity checks: associations known from population studies -------------------------------------------
sanity <- map(lc$expected_associations %||% list(), \(e) {
  row <- tibble(protein = e$protein, parameter = e$parameter, label = param_label(e$parameter), expected = e$direction)
  oid <- find_assay(det, e$protein)[1]
  if (is.na(oid)) return(row |> mutate(status = "protein not measured"))
  if (!e$parameter %in% names(S) || sum(!is.na(S[[e$parameter]])) < min_n) return(row |> mutate(status = "parameter not available"))
  r <- spearman_vs(S[[e$parameter]], Yall[, oid, drop = FALSE], min_n = min_n)
  ok <- isTRUE(r$p < 0.05 && sign(r$rho) == if (e$direction == "negative") -1 else 1)
  row |> mutate(frac_above_lod = det$frac_above_lod[det$OlinkID == oid], n = r$n, rho = r$rho, p = r$p,
                status = if (ok) "recovered" else "not recovered")
}) |> bind_rows()
if (nrow(sanity)) {
  sanity <- ensure(sanity, c("frac_above_lod", "n", "rho", "p")) |> relocate(status, .after = last_col())
  save_csv(sanity, cfg, "leip_biobank", "sanity_checks.csv")
}

# ---- 5. galanin -------------------------------------------------------------------------------------------------
gal_oid <- find_assay(det, gal_name)
if (length(gal_oid) > 1) {
  msg("%s: several assays match (%s) - using %s", gal_name, paste(det$Assay[det$OlinkID %in% gal_oid], collapse = ", "), gal_oid[1])
  gal_oid <- gal_oid[1]
}
has_gal   <- length(gal_oid) == 1
has_elisa <- "galanin_elisa" %in% names(S) && sum(!is.na(S$galanin_elisa)) >= 10
if (!has_gal) msg("WARNING: galanin (%s) is not in the LEIP Olink data - galanin analyses skipped", gal_name)
if (!has_elisa) msg("WARNING: no galanin ELISA values (column '%s') - Olink vs ELISA comparison skipped",
                    lc$galanin_lab %||% "Galanin [pg/mL]")
gd <- S |> select(SampleID, SubjectID, plate, any_of(c("age", "sex_male", "galanin_elisa", "Galanin_ELISA_plate",
                                                       "low_serum", "lipamisch")))
gal_det <- NULL
if (has_gal) {
  gd <- gd |> mutate(npx = Yall[SampleID, gal_oid]) |>
    left_join(lp |> filter(OlinkID == gal_oid) |> select(SampleID, below_lod, LOD), by = "SampleID")
  gal_det <- det |> filter(OlinkID == gal_oid)
  msg("Galanin: %s (%s), %.0f%% of LEIP samples above LOD", gal_det$Assay, gal_oid, 100 * gal_det$frac_above_lod)
}
if (!is.null(sq)) gd <- gd |> left_join(sq |> select(SampleID, SampleQC, qc_outlier = outlier), by = "SampleID")
adj_txt_elisa <- adj_lab(setdiff(adj, "plate"))

# 5a. Olink GAL vs ELISA: agreement, tertiles, sensitivity analyses
agree <- tert_tab <- tert_stats <- NULL
if (has_gal && has_elisa) {
  flagged <- coalesce(col_or(gd, "SampleQC", "PASS") != "PASS", FALSE) | coalesce(as.logical(col_or(gd, "qc_outlier", FALSE)), FALSE) |
    coalesce(col_or(gd, "low_serum", 0) == 1, FALSE) | coalesce(col_or(gd, "lipamisch", 0) == 1, FALSE)
  agree_row <- function(what, keep, covs = NULL, min_pairs = 10) {
    d <- gd[keep & !is.na(gd$npx) & coalesce(gd$galanin_elisa > 0, FALSE), ]
    if (!is.null(covs)) d <- d[stats::complete.cases(d[covs]), ]
    if (nrow(d) < min_pairs) return(tibble(analysis = what, n = nrow(d)))
    r <- spearman_vs(d$galanin_elisa, matrix(d$npx, dimnames = list(NULL, "GAL")), if (!is.null(covs)) d[covs], min_n = min_pairs)
    f <- tryCatch(lm(reformulate(c("log2(galanin_elisa)", covs), "npx"), data = d), error = \(e) NULL)
    ci <- if (!is.null(f)) suppressWarnings(confint(f)[2, ]) else c(NA, NA)
    tibble(analysis = what, n = r$n, rho = r$rho, ci_low = r$ci_low, ci_high = r$ci_high, p = r$p,
           pearson_r_log2 = if (is.null(covs)) cor(log2(d$galanin_elisa), d$npx) else NA_real_,
           npx_per_doubling = if (!is.null(f)) unname(coef(f)[2]) else NA_real_,
           npx_per_doubling_ci_low = unname(ci[1]), npx_per_doubling_ci_high = unname(ci[2]))
  }
  agree <- bind_rows(
    agree_row("all samples", TRUE),
    agree_row("adjusted for Olink plate", TRUE, "plate"),
    agree_row("only samples with GAL above LOD", coalesce(!gd$below_lod, FALSE)),
    agree_row("without flagged samples (QC warning / outlier, low serum, lipaemic)", !flagged),
    map(sort(unique(na.omit(gd$plate))), \(pl) agree_row(paste("within", pl), coalesce(gd$plate == pl, FALSE), min_pairs = 5))) |>
    ensure(c("rho", "ci_low", "ci_high", "p", "pearson_r_log2", "npx_per_doubling"))
  tt <- gd |> filter(!is.na(npx), !is.na(galanin_elisa)) |>
    mutate(ELISA = factor(ntile(galanin_elisa, 3), 1:3, c("low", "middle", "high")),
           Olink = factor(ntile(npx, 3), 1:3, c("low", "middle", "high")))
  tab <- table(ELISA = tt$ELISA, Olink = tt$Olink)
  tert_tab <- as.data.frame.matrix(tab) |> rownames_to_column("ELISA tertile (rows) / Olink tertile (columns)")
  tert_stats <- tibble(n = sum(tab), same_tertile_pct = 100 * sum(diag(tab)) / sum(tab),
                       opposite_tertile_pct = 100 * (tab[1, 3] + tab[3, 1]) / sum(tab), weighted_kappa = weighted_kappa(tab))
}

# 5b. technical factors: were the samples put on the plates by galanin level? does the plate shift GAL?
kw <- function(v, g, what) {
  ok <- !is.na(v) & !is.na(g)
  if (n_distinct(g[ok]) < 2) return(NULL)
  m <- tapply(v[ok], as.character(g[ok]), median)
  tibble(test = what, n = sum(ok), groups = length(m), medians = paste(sprintf("%s: %.3g", names(m), m), collapse = "; "),
         p = kruskal.test(v[ok], factor(g[ok]))$p.value)
}
plate_fx <- bind_rows(
  if (has_elisa) kw(gd$galanin_elisa, gd$plate, "galanin ELISA by Olink plate (sample allocation)"),
  if (has_gal) kw(gd$npx, gd$plate, "Olink GAL by Olink plate"),
  if (has_elisa && "Galanin_ELISA_plate" %in% names(gd)) kw(gd$galanin_elisa, gd$Galanin_ELISA_plate, "galanin ELISA by ELISA plate"))

# 5c. which proteins does the ELISA correlate with, and where does Olink GAL rank among them?
ev <- NULL; gal_rank <- n_rank <- NA
if (has_elisa) {
  Ye <- if (has_gal && !gal_oid %in% colnames(Y)) cbind(Y, Yall[, gal_oid, drop = FALSE]) else Y
  u <- spearman_vs(S$galanin_elisa, Ye, min_n = min_n)
  a <- spearman_vs(S$galanin_elisa, Ye, cov_for("galanin_elisa"), min_n = min_n)
  ev <- tibble(OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
               n_adj = a$n, rho_adj = a$rho, p_adj = a$p) |>
    mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH"), rank_by_rho = rank(-rho, ties.method = "min", na.last = "keep")) |>
    left_join(det |> select(OlinkID, Assay, protein, frac_above_lod, analysed), by = "OlinkID") |>
    relocate(OlinkID, Assay, protein) |> arrange(p)
  n_rank <- sum(!is.na(ev$rho))
  if (has_gal) gal_rank <- ev$rank_by_rho[ev$OlinkID == gal_oid]
}

# 5d. benchmark: every protein measured by both the lab and Olink (how well do lab and Olink agree in general?)
bm_map <- lc$lab_vs_olink %||% setNames(list(gal_name), "galanin_elisa")
bench <- imap(bm_map, \(assay, lab) {
  row <- tibble(lab = lab, label = param_label(lab), olink_assay = assay)
  if (!lab %in% names(S)) return(row |> mutate(status = "not in the clinical file"))
  if (sum(!is.na(S[[lab]])) < min_n) return(row |> mutate(status = sprintf("fewer than %d values", min_n)))
  oid <- find_assay(det, assay)[1]
  if (is.na(oid)) return(row |> mutate(status = "not measured by Olink"))
  allr <- spearman_vs(S[[lab]], if (oid %in% colnames(Y)) Y else cbind(Y, Yall[, oid, drop = FALSE]), min_n = min_n)
  r <- allr |> filter(id == oid)
  row |> mutate(status = "compared", OlinkID = oid, n = r$n, frac_above_lod = det$frac_above_lod[det$OlinkID == oid],
                rho = r$rho, ci_low = r$ci_low, ci_high = r$ci_high, p = r$p,
                rank_among_proteins = rank(-allr$rho, ties.method = "min", na.last = "keep")[allr$id == oid],
                n_proteins = sum(!is.na(allr$rho)))
}) |> bind_rows() |> ensure(c("n", "frac_above_lod", "rho", "ci_low", "ci_high", "p", "rank_among_proteins", "n_proteins")) |>
  arrange(desc(status == "compared"), desc(rho))

# 5e. galanin (Olink and ELISA) vs every clinical parameter; p-value primary (pre-specified protein)
clin_cor <- function(v, measure, drop_cov = character()) {
  V <- matrix(v, dimnames = list(NULL, measure))
  map(setdiff(params, if (measure == "galanin ELISA") "galanin_elisa"), \(pm) {
    u <- spearman_vs(S[[pm]], V, min_n = min_n)
    a <- spearman_vs(S[[pm]], V, cov_for(pm, drop_cov), min_n = min_n)
    tibble(measure = measure, parameter = pm, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
           n_adj = a$n, rho_adj = a$rho, p_adj = a$p)
  }) |> bind_rows() |>
    mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH")) |>
    left_join(pinfo |> select(parameter, label, type, key), by = "parameter") |>
    relocate(measure, parameter, label, type, key)
}
gcl <- bind_rows(if (has_gal) clin_cor(gd$npx, "Olink GAL"), if (has_elisa) clin_cor(gd$galanin_elisa, "galanin ELISA", "plate"))
prof <- NULL; prof_r <- NA
if (has_gal && has_elisa) {
  prof <- gcl |> filter(parameter != "galanin_elisa") |> select(parameter, label, measure, rho) |>
    pivot_wider(names_from = measure, values_from = rho)
  prof_r <- suppressWarnings(cor(prof$`Olink GAL`, prof$`galanin ELISA`, method = "spearman", use = "complete.obs"))
}

# 5f. which proteins does Olink GAL correlate with (+ gene-set enrichment of the ranking)
gp <- gsea <- NULL
if (has_gal) {
  Yp <- Y[, setdiff(colnames(Y), gal_oid), drop = FALSE]
  u <- spearman_vs(gd$npx, Yp, min_n = min_n)
  a <- spearman_vs(gd$npx, Yp, cov_for(), min_n = min_n)
  gp <- tibble(OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p,
               n_adj = a$n, rho_adj = a$rho, p_adj = a$p) |>
    mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH")) |>
    left_join(det |> select(OlinkID, Assay, protein, frac_above_lod), by = "OlinkID") |>
    relocate(OlinkID, Assay, protein) |> arrange(p)
  if (isTRUE(lc$gsea %||% TRUE) && sum(!is.na(gp$rho)) >= 50) {
    gsea <- tryCatch({
      sets <- load_gene_sets(cfg$enrichment$collections)
      st <- gp |> filter(!is.na(rho)) |> transmute(gene = str_split(Assay, "_"), rho) |> unnest(gene) |>
        group_by(gene) |> slice_max(abs(rho), n = 1, with_ties = FALSE) |> ungroup()
      out <- suppressWarnings(fgsea::fgseaMultilevel(sets, setNames(st$rho, st$gene), minSize = cfg$enrichment$min_size,
                                                     maxSize = cfg$enrichment$max_size, eps = 0))
      as_tibble(out) |> mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) |> arrange(pval)
    }, error = \(e) { msg("Gene-set enrichment for galanin skipped: %s", conditionMessage(e)); NULL })
  }
}

# ---- 6. answers -----------------------------------------------------------------------------------------------------
answers <- list()
answer <- function(question, item, verdict, evidence) answers[[length(answers) + 1]] <<- tibble(question, item, verdict, evidence)
if (has_gal) {
  fa <- gal_det$frac_above_lod
  answer("G1 Is galanin (Olink GAL) measurable in LEIP serum?", "Olink GAL",
         case_when(is.na(fa) ~ "unknown (no LOD available)", fa >= 0.9 ~ "yes, in (nearly) all samples",
                   fa >= cfg$qc$min_detect_frac ~ "yes, in most samples", fa >= 0.2 ~ "partly: many samples are below LOD",
                   TRUE ~ "mostly below LOD: the values are dominated by noise"),
         sprintf("%d of %d LEIP samples above LOD (%.0f%%); median NPX %.2f, median LOD %.2f.%s",
                 sum(!gd$below_lod, na.rm = TRUE), sum(!is.na(gd$below_lod)), 100 * fa, gal_det$median_npx, gal_det$median_lod,
                 if (!gal_det$analysed) " GAL fails the detection filter, so it is not in the proteome-wide tables; the galanin analyses use it anyway (pre-specified)." else ""))
}
if (has_gal && has_elisa) {
  a0 <- agree |> filter(analysis == "all samples"); a1 <- agree |> filter(analysis == "adjusted for Olink plate")
  others <- bench |> filter(status == "compared", lab != "galanin_elisa", !is.na(rho))
  pl <- if (nrow(plate_fx)) plate_fx |> filter(str_detect(test, "^galanin ELISA by Olink plate")) else tibble()
  answer("G2 Do the Olink GAL levels correspond to the galanin ELISA?", "Olink GAL vs ELISA",
         case_when(is.na(a0$p) ~ "not enough data",
                   a0$p < 0.05 & a0$rho >= 0.5 ~ "yes: both methods rank the samples similarly",
                   a0$p < 0.05 & a0$rho >= 0.3 ~ "partly: moderate agreement",
                   a0$p < 0.05 & a0$rho < 0 ~ "no: inverse relation (unexpected)",
                   TRUE ~ "no clear agreement"),
         paste0(sprintf("Spearman rho = %.2f (95%% CI %.2f to %.2f), p = %s, n = %d. ", a0$rho, a0$ci_low, a0$ci_high, fmt_p(a0$p), a0$n),
                sprintf("Adjusted for Olink plate: rho = %.2f, p = %s. ", a1$rho, fmt_p(a1$p)),
                sprintf("Same tertile in %.0f%% of persons (chance: 33%%), opposite tertiles in %.0f%%; weighted kappa %.2f. ",
                        tert_stats$same_tertile_pct, tert_stats$opposite_tertile_pct, tert_stats$weighted_kappa),
                if (!is.na(gal_rank)) sprintf("Of %d proteins, GAL ranks %d in correlation with the ELISA. ", n_rank, gal_rank) else "",
                if (nrow(pl) && isTRUE(pl$p[1] < 0.05)) sprintf("Note: the ELISA values differ between the Olink plates (p = %s), so plate and galanin level are partly confounded; see the plate-adjusted result. ", fmt_p(pl$p[1])) else "",
                if (nrow(others)) sprintf("Benchmark - other proteins measured by the lab and by Olink: median rho %.2f (%s).",
                                          median(others$rho), paste(sprintf("%s %.2f", others$olink_assay, others$rho), collapse = ", ")) else ""))
}
for (m in if (nrow(gcl)) unique(gcl$measure)) {
  d <- gcl |> filter(measure == m, parameter != "galanin_elisa") |> arrange(p)
  hits <- d |> filter(p < 0.05)
  answer("G3 Which clinical parameters does galanin correlate with?", m,
         if (nrow(hits)) sprintf("%d of %d parameters at p < 0.05 (%d at FDR < %g across parameters)", nrow(hits), nrow(d),
                                 sum(d$fdr < fdr_cut, na.rm = TRUE), fdr_cut) else sprintf("none of %d parameters at p < 0.05", nrow(d)),
         if (nrow(hits)) paste0(paste(head(sprintf("%s %+.2f (p = %s; adjusted %+.2f, p = %s)", hits$label, hits$rho, fmt_p(hits$p),
                                                  hits$rho_adj, fmt_p(hits$p_adj)), 12), collapse = "; "),
                                sprintf(". Adjusted = partial Spearman for %s. By chance alone about %.1f of %d parameters reach p < 0.05.",
                                        if (m == "Olink GAL") adj_txt else adj_txt_elisa, 0.05 * nrow(d), nrow(d)))
         else sprintf("Strongest: %s.", paste(head(sprintf("%s %+.2f (p = %s)", d$label, d$rho, fmt_p(d$p)), 3), collapse = "; ")))
}
if (!is.na(prof_r))
  answer("G3 Which clinical parameters does galanin correlate with?", "Olink and ELISA: same clinical pattern?",
         sprintf("the two correlation profiles agree with rho = %.2f (%d parameters)", prof_r, nrow(prof)),
         "Descriptive: compares the rho of Olink GAL and of the ELISA with each clinical parameter (the parameters are correlated, so there is no p-value).")
if (has_gal) {
  s <- gp |> filter(fdr < fdr_cut); pos <- gp |> filter(rho > 0); neg <- gp |> filter(rho < 0)
  gs_sig <- if (!is.null(gsea)) gsea |> filter(padj < fdr_cut) else NULL
  answer("G4 Which proteins does galanin correlate with?", "Olink GAL vs all proteins",
         sprintf("%d of %d proteins at FDR < %g (%d positive, %d negative); adjusted for %s: %d", nrow(s), sum(!is.na(gp$p)), fdr_cut,
                 sum(s$rho > 0), sum(s$rho < 0), adj_txt, sum(gp$fdr_adj < fdr_cut, na.rm = TRUE)),
         paste0(sprintf("Strongest positive: %s. Strongest negative: %s.", top_str(pos$protein, pos$rho, pos$p, 10), top_str(neg$protein, neg$rho, neg$p, 5)),
                if (!is.null(gs_sig) && nrow(gs_sig)) sprintf(" Gene sets (GSEA on the correlation ranking, padj < %g): %s.", fdr_cut,
                                                              paste(head(sprintf("%s (NES %+.1f)", gs_sig$pathway, gs_sig$NES), 6), collapse = ", "))
                else if (!is.null(gsea)) " No gene set enriched (GSEA padj >= 0.05)." else ""))
}
if (has_elisa) {
  s <- ev |> filter(fdr < fdr_cut)
  answer("G4 Which proteins does galanin correlate with?", "galanin ELISA vs all proteins",
         sprintf("%d of %d proteins at FDR < %g", nrow(s), n_rank, fdr_cut),
         sprintf("Strongest: %s.", top_str(ev$protein, ev$rho, ev$p, 10)))
}
hit_par <- nsig |> filter(significant > 0) |> arrange(desc(significant))
answer("P1 Which clinical parameters are reflected in the LEIP serum proteome?", "all proteins x all parameters",
       sprintf("%d of %d parameters with proteins at FDR < %g (adjusted for %s: %d)", nrow(hit_par), nrow(nsig), fdr_cut, adj_txt,
               sum(nsig$significant_adjusted > 0)),
       if (nrow(hit_par)) paste(head(sprintf("%s: %d at FDR < %g, strongest %s", hit_par$label, hit_par$significant, fdr_cut, hit_par$top_proteins), 10),
                                collapse = "; ")
       else "No protein passes FDR for any parameter.")
if (nrow(sanity)) {
  tested <- sanity |> filter(status %in% c("recovered", "not recovered"))
  answer("P2 Do known associations show up? (sanity check)", "associations known from population studies",
         sprintf("%d of %d recovered (p < 0.05, expected direction)%s", sum(tested$status == "recovered"), nrow(tested),
                 if (nrow(tested) < nrow(sanity)) sprintf("; %d not testable", nrow(sanity) - nrow(tested)) else ""),
         paste(sprintf("%s ~ %s: %s", sanity$protein, sanity$label,
                       if_else(is.na(sanity$rho), sanity$status, sprintf("rho %+.2f, p = %s", sanity$rho, fmt_p(sanity$p)))), collapse = "; "))
}
answers <- bind_rows(answers)
save_csv(answers, cfg, "leip_biobank", "answers.csv")

# ---- 7. figures --------------------------------------------------------------------------------------------------
figs <- list()
keep_fig <- function(name, p, ..., width = 8, height = 6) { save_plot(p, cfg, ..., width = width, height = height); figs[[name]] <<- p }

nb_par <- nsig |> filter(significant + significant_adjusted > 0 | key)
if (nrow(nb_par)) {
  nb <- nb_par |>
    select(label, `unadjusted|positive` = positive, `unadjusted|negative` = negative,
           `adjusted|positive` = positive_adj, `adjusted|negative` = negative_adj) |>
    pivot_longer(-label, names_to = c("analysis", "direction"), names_sep = "\\|", values_to = "n") |>
    mutate(n = if_else(direction == "negative", -n, n),
           label = factor(label, levels = nb_par |> arrange(significant, significant_adjusted) |> pull(label)),
           analysis = factor(analysis, c("unadjusted", "adjusted"), c("Spearman", paste("adjusted for", adj_txt))))
  p <- ggplot(nb, aes(n, label, fill = direction)) + geom_col() + geom_vline(xintercept = 0, colour = "grey40") +
    facet_wrap(~analysis) + scale_fill_manual(values = c(positive = "firebrick", negative = "steelblue")) +
    labs(title = sprintf("LEIP: proteins correlated with each clinical parameter (FDR < %g)", fdr_cut),
         subtitle = sprintf("%d proteins tested per parameter; key parameters and all parameters with hits", ncol(Y)),
         x = "proteins (negative < 0 < positive)", y = NULL, fill = NULL)
  keep_fig("nsig", p, "leip_biobank", "n_significant_per_parameter.png", width = 11, height = 2 + 0.2 * nrow(nb_par))
}

hk <- assoc |> filter(key, !is.na(p))
if (nrow(hk)) {
  top_hm <- hk |> group_by(OlinkID) |> summarise(best = min(p)) |> slice_min(best, n = 50, with_ties = FALSE) |> pull(OlinkID)
  d <- hk |> filter(OlinkID %in% top_hm)
  m <- d |> select(protein, label, rho) |> pivot_wider(names_from = label, values_from = rho)
  mm <- as.matrix(m[, -1]); mm[is.na(mm)] <- 0
  ord <- if (nrow(mm) > 2) m$protein[hclust(dist(mm))$order] else m$protein
  p <- ggplot(d |> mutate(protein = factor(protein, levels = ord), label = factor(label, levels = unique(pinfo$label[pinfo$key]))),
              aes(label, protein, fill = rho)) + geom_tile() +
    geom_text(aes(label = if_else(significant, "*", "")), size = 4, vjust = 0.75) +
    scale_fill_gradient2(low = "steelblue", high = "firebrick", limits = c(-1, 1)) +
    labs(title = "LEIP: correlation of proteins with the key clinical parameters",
         subtitle = sprintf("the %d proteins with the strongest association (Spearman rho; * FDR < %g). All proteins: associations_all.csv.gz",
                            length(top_hm), fdr_cut), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  keep_fig("heatmap", p, "leip_biobank", "heatmap_key_parameters.png", width = 10, height = 3 + 0.2 * length(top_hm))

  lab_d <- hk |> group_by(label) |> slice_min(p, n = 3, with_ties = FALSE) |> ungroup()
  p <- ggplot(hk, aes(rho, -log10(p))) + geom_point(aes(colour = significant), size = 0.6, alpha = 0.6) +
    geom_text(data = lab_d, aes(label = protein), size = 2.6, vjust = -0.6, check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"),
                        labels = c(`FALSE` = "not significant", `TRUE` = sprintf("FDR < %g", fdr_cut)), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    facet_wrap(~label) + labs(title = "LEIP: every protein vs each key clinical parameter", x = "Spearman rho", y = "-log10 p") +
    theme(legend.position = "bottom")
  keep_fig("volcano", p, "leip_biobank", "volcano_key_parameters.png", width = 13, height = 10)

  tp <- hk |> slice_min(p, n = 12, with_ties = FALSE)
  sc <- map(seq_len(nrow(tp)), \(i) tibble(panel = sprintf("%s vs %s\nrho = %.2f, p = %s", tp$protein[i], tp$label[i], tp$rho[i], fmt_p(tp$p[i])),
                                           x = S[[tp$parameter[i]]], y = Yall[, tp$OlinkID[i]])) |> bind_rows() |>
    mutate(panel = factor(panel, levels = unique(panel))) |> filter(!is.na(x), !is.na(y))
  p <- ggplot(sc, aes(x, y)) + geom_point(size = 1.2) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
    facet_wrap(~panel, scales = "free") +
    labs(title = "LEIP: the strongest protein - key parameter correlations", x = "clinical parameter", y = "NPX")
  keep_fig("top_pairs", p, "leip_biobank", "top_associations.png", width = 12, height = 9)
}

if (has_gal && has_elisa) {
  a0 <- agree |> filter(analysis == "all samples"); a1 <- agree |> filter(analysis == "adjusted for Olink plate")
  d <- gd |> filter(!is.na(npx), coalesce(galanin_elisa > 0, FALSE)) |> mutate(lod = if_else(coalesce(below_lod, FALSE), "below LOD", "above LOD"))
  p <- ggplot(d, aes(galanin_elisa, npx)) +
    (if (!is.na(gal_det$median_lod)) geom_hline(yintercept = gal_det$median_lod, linetype = 2, colour = "grey50")) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", linewidth = 0.6) +
    geom_point(aes(colour = plate, shape = lod), size = 2.4) +
    scale_x_continuous(trans = "log2", breaks = pretty(d$galanin_elisa, n = 6)) +
    scale_shape_manual(values = c(`above LOD` = 16, `below LOD` = 1)) +
    labs(title = "Galanin in LEIP serum: Olink (GAL) vs ELISA",
         subtitle = sprintf("Spearman rho = %.2f (95%% CI %.2f to %.2f), p = %s, n = %d; adjusted for Olink plate: rho = %.2f, p = %s\ndashed line: median Olink LOD",
                            a0$rho, a0$ci_low, a0$ci_high, fmt_p(a0$p), a0$n, a1$rho, fmt_p(a1$p)),
         x = "galanin ELISA (pg/mL, log2 scale)", y = "Olink GAL (NPX, log2 scale)", colour = "Olink plate", shape = "Olink")
  keep_fig("gal_elisa", p, "leip_biobank", "galanin", "GAL_Olink_vs_ELISA.png", width = 8.5, height = 6)
}

b <- bench |> filter(status == "compared", !is.na(rho))
if (nrow(b)) {
  b <- b |> mutate(name = sprintf("%s  vs  Olink %s", label, olink_assay), galanin = lab == "galanin_elisa",
                   info = sprintf("%s above LOD, rank %d of %d", if_else(is.na(frac_above_lod), "?", sprintf("%.0f%%", 100 * frac_above_lod)),
                                  rank_among_proteins, n_proteins))
  p <- ggplot(b, aes(rho, reorder(name, rho), colour = galanin)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.4) +
    geom_text(aes(x = 1.05, label = info), hjust = 0, size = 3, colour = "grey30") +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "firebrick"), guide = "none") +
    scale_x_continuous(breaks = seq(-1, 1, 0.5)) + coord_cartesian(xlim = c(-1, 2)) +
    labs(title = "Lab measurement vs Olink for the same protein (LEIP serum)",
         subtitle = "Spearman rho with 95% CI. Rank = place of the matching Olink assay among all proteins correlated with the lab value (1 = best)",
         x = "Spearman rho", y = NULL)
  keep_fig("bench", p, "leip_biobank", "galanin", "lab_vs_Olink_benchmark.png", width = 11, height = 2 + 0.35 * nrow(b))
}

gal_clin_plot <- function(d, title_extra = "") {
  first_m <- if (has_gal) "Olink GAL" else "galanin ELISA"
  lev <- d |> filter(measure == first_m) |> arrange(rho) |> pull(label)
  lev <- c(lev, setdiff(unique(d$label), lev))
  d <- d |> mutate(y = match(label, lev) + if_else(measure == "Olink GAL", 0.18, -0.18), sig = coalesce(p < 0.05, FALSE))
  ggplot(d, aes(rho, y, colour = measure)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), alpha = 0.45) + geom_point(aes(shape = sig), size = 2) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    scale_colour_manual(values = c(`Olink GAL` = "firebrick", `galanin ELISA` = "steelblue"), name = NULL) +
    scale_y_continuous(breaks = seq_along(lev), labels = lev, expand = expansion(add = 0.6)) +
    labs(title = paste0("What does galanin correlate with? Clinical parameters in LEIP", title_extra),
         subtitle = sprintf("Spearman rho with 95%% CI; filled = p < 0.05 (single tests)%s",
                            if (!is.na(prof_r)) sprintf("; agreement of the Olink and ELISA profiles: rho = %.2f", prof_r) else ""),
         x = "Spearman rho", y = NULL) + theme(legend.position = "bottom")
}
if (nrow(gcl)) {
  p <- gal_clin_plot(gcl)
  keep_fig("gal_clin", p, "leip_biobank", "galanin", "galanin_vs_clinical.png", width = 9, height = 2.5 + 0.2 * n_distinct(gcl$label))
  sel <- gcl |> group_by(parameter) |> filter(any(coalesce(p < 0.05, FALSE)) | any(key)) |> ungroup()
  figs$gal_clin_short <- gal_clin_plot(sel, " (key parameters and p < 0.05)")
}

if (has_gal) {
  tc <- gcl |> filter(measure == "Olink GAL", parameter != "galanin_elisa", !is.na(p)) |> slice_min(p, n = 6, with_ties = FALSE)
  if (nrow(tc)) {
    sc <- map(seq_len(nrow(tc)), \(i) tibble(panel = sprintf("%s: rho = %.2f, p = %s", tc$label[i], tc$rho[i], fmt_p(tc$p[i])),
                                             x = S[[tc$parameter[i]]], y = gd$npx)) |> bind_rows() |>
      mutate(panel = factor(panel, levels = unique(panel))) |> filter(!is.na(x), !is.na(y))
    p <- ggplot(sc, aes(x, y)) + geom_point(size = 1.5) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
      facet_wrap(~panel, scales = "free_x") +
      labs(title = "Olink GAL vs its most strongly correlated clinical parameters (LEIP)", x = NULL, y = "Olink GAL (NPX)")
    keep_fig("gal_clin_scatter", p, "leip_biobank", "galanin", "GAL_top_clinical_scatter.png", width = 10, height = 6.5)
  }
  gv <- gp |> filter(!is.na(p))
  fdr_line <- if (any(gv$fdr < fdr_cut)) max(gv$p[gv$fdr < fdr_cut]) else NA
  p <- ggplot(gv, aes(rho, -log10(p))) +
    (if (!is.na(fdr_line)) geom_hline(yintercept = -log10(fdr_line), linetype = 2, colour = "grey50")) +
    geom_point(aes(colour = fdr < fdr_cut), size = 0.8, alpha = 0.7) +
    geom_text(data = head(gv, 20), aes(label = protein), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    scale_colour_manual(values = c(`FALSE` = "grey65", `TRUE` = "firebrick"),
                        labels = c(`FALSE` = "not significant", `TRUE` = sprintf("FDR < %g", fdr_cut)), name = NULL) +
    labs(title = "Which proteins correlate with galanin (Olink GAL) in LEIP serum?",
         subtitle = sprintf("Spearman rho vs -log10 p for %d proteins; dashed line: FDR cutoff", nrow(gv)), x = "Spearman rho", y = "-log10 p") +
    theme(legend.position = "bottom")
  keep_fig("gal_volcano", p, "leip_biobank", "galanin", "GAL_vs_proteins_volcano.png", width = 9, height = 7)
  tq <- head(gv, 9)
  sc <- map(seq_len(nrow(tq)), \(i) tibble(panel = sprintf("%s: rho = %.2f", tq$protein[i], tq$rho[i]), x = Yall[, tq$OlinkID[i]], y = gd$npx)) |>
    bind_rows() |> mutate(panel = factor(panel, levels = unique(panel))) |> filter(!is.na(x), !is.na(y))
  p <- ggplot(sc, aes(x, y)) + geom_point(size = 1.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
    facet_wrap(~panel, scales = "free_x") +
    labs(title = "Olink GAL vs its most strongly correlated proteins (LEIP)", x = "NPX of the other protein", y = "Olink GAL (NPX)")
  keep_fig("gal_partners", p, "leip_biobank", "galanin", "GAL_top_proteins_scatter.png", width = 10, height = 8)
  if (!is.null(gsea) && nrow(gsea)) {
    g <- gsea |> slice_min(pval, n = 20, with_ties = FALSE)
    p <- ggplot(g, aes(NES, reorder(str_trunc(pathway, 60), NES), size = size, colour = padj < fdr_cut)) + geom_point() +
      geom_vline(xintercept = 0, colour = "grey60") +
      scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "firebrick"), labels = c(`FALSE` = "padj >= 0.05", `TRUE` = "padj < 0.05"), name = NULL) +
      labs(title = "Gene sets among the proteins correlated with galanin (GSEA)",
           subtitle = "NES > 0: gene set enriched among proteins positively correlated with GAL; top 20 by p-value", x = "normalised enrichment score", y = NULL)
    keep_fig("gal_gsea", p, "leip_biobank", "galanin", "GAL_gene_sets.png", width = 10, height = 7)
  }
}
if (has_elisa) {
  e2 <- ev |> filter(!is.na(p)) |> mutate(is_gal = OlinkID %in% gal_oid)
  p <- ggplot(e2, aes(rho, -log10(p))) + geom_point(colour = "grey60", size = 0.8, alpha = 0.7) +
    geom_text(data = head(e2, 12), aes(label = protein), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
    geom_point(data = e2 |> filter(is_gal), colour = "firebrick", size = 3) +
    geom_text(data = e2 |> filter(is_gal), aes(label = protein), colour = "firebrick", vjust = 1.8, fontface = "bold") +
    labs(title = "Which Olink proteins correlate with the galanin ELISA?",
         subtitle = if (!is.na(gal_rank)) sprintf("red: Olink GAL, rank %d of %d proteins by correlation", gal_rank, n_rank) else "Olink GAL not measured",
         x = "Spearman rho", y = "-log10 p")
  keep_fig("elisa_volcano", p, "leip_biobank", "galanin", "ELISA_vs_proteins_volcano.png", width = 9, height = 7)
}

# ---- 8. tables -----------------------------------------------------------------------------------------------------
samples_out <- S |> select(SampleID, SubjectID, plate, any_of(setdiff(cl$key, leip_technical))) |>
  left_join(id_check |> select(SampleID, problem), by = "SampleID")
if (has_gal) samples_out <- samples_out |> mutate(GAL_NPX = gd$npx, GAL_below_LOD = gd$below_lod)
save_csv(samples_out, cfg, "leip_biobank", "leip_samples.csv")
sheets <- list(answers = answers, parameters = pinfo, samples = samples_out, sample_check = id_check, detection = det,
               n_significant = nsig, significant = sig_tab, top10_per_parameter = top10, sanity_checks = sanity)
writexl::write_xlsx(keep(sheets, \(d) is.data.frame(d) && ncol(d) > 0), out_path(cfg, "leip_biobank", "leip_biobank.xlsx"))
if (has_gal || has_elisa) {
  gal_sheets <- list(answers = answers |> filter(str_detect(question, "^G")), detection = gal_det, agreement = agree,
                     tertiles = tert_tab, tertile_agreement = tert_stats, plate_effects = plate_fx, lab_vs_olink = bench,
                     galanin_vs_clinical = gcl, clinical_profiles = prof, GAL_vs_proteins = gp, ELISA_vs_proteins = ev,
                     GAL_gene_sets = gsea, sample_values = gd)
  writexl::write_xlsx(keep(compact(gal_sheets), \(d) is.data.frame(d) && ncol(d) > 0), out_path(cfg, "leip_biobank", "galanin", "galanin.xlsx"))
  if (!is.null(agree)) save_csv(agree, cfg, "leip_biobank", "galanin", "olink_vs_elisa.csv")
  save_csv(bench, cfg, "leip_biobank", "galanin", "lab_vs_olink_benchmark.csv")
  if (nrow(gcl)) save_csv(gcl, cfg, "leip_biobank", "galanin", "galanin_vs_clinical.csv")
  if (!is.null(gp)) save_csv(gp, cfg, "leip_biobank", "galanin", "GAL_vs_proteins.csv")
  if (!is.null(ev)) save_csv(ev, cfg, "leip_biobank", "galanin", "ELISA_vs_proteins.csv")
}

# ---- 9. PDF summary ----------------------------------------------------------------------------------------------------
pdf_file <- out_path(cfg, "leip_biobank", "LEIP_biobank_summary.pdf")
pdf_open(pdf_file)
section("Overview", {
  wc <- S |> filter(SampleID %in% id_check$SampleID[coalesce(id_check$n_clinical_values, 0) > 0])
  q <- \(v, d = 1) if (all(is.na(v))) "n/a" else
    sprintf(paste0("%.", d, "f (IQR %.", d, "f-%.", d, "f)"), median(v, na.rm = TRUE), quantile(v, 0.25, na.rm = TRUE), quantile(v, 0.75, na.rm = TRUE))
  sx <- col_or(wc, "sex_male")
  left_out <- pinfo |> filter(!analysed) |> mutate(r = str_replace(reason, "^same ranks as.*", "same ranks as another variable")) |> count(r)
  n_med <- if (n_clin) n_clin else nrow(S)
  probs <- id_check |> filter(!is.na(problem)) |> group_by(problem) |>
    summarise(which = if (n() <= 5) paste(SampleID, collapse = ", ") else sprintf("%d samples", n()), .groups = "drop")
  page_text("LEIP biobank serum: proteins, clinical parameters and galanin", c(
    "## Samples",
    sprintf("%d LEIP serum samples in the Olink data after QC, %d of them with clinical data (SORB database).", nrow(S), nrow(wc)),
    sprintf("Age %s years; %d men, %d women; BMI %s; galanin ELISA %s pg/mL (medians).", q(col_or(wc, "age")), sum(sx == 1, na.rm = TRUE),
            sum(sx == 0, na.rm = TRUE), q(col_or(wc, "BMI")), q(col_or(wc, "galanin_elisa"), 0)),
    if (nrow(probs)) paste0("Sample checks: ", paste(sprintf("%s (%s)", probs$problem, probs$which), collapse = "; "), ". Details: sample_check.csv.")
    else "Manifest and clinical file agree for every sample (SubjectID, Olink plate).",
    "## Proteins",
    sprintf("%d of %d proteins are measurable in LEIP serum (>= %.0f%% of the LEIP samples above LOD) and are analysed.",
            ncol(Y), nrow(det), 100 * cfg$qc$min_detect_frac),
    "## Clinical parameters",
    sprintf("%d of %d variables of the clinical file are analysed, %d of them key parameters (sheet Key_parameters). Left out: %s. Full list with reasons: parameters.csv.",
            length(params), nrow(pinfo), sum(pinfo$analysed & pinfo$key), paste(sprintf("%s (%d)", left_out$r, left_out$n), collapse = "; ")),
    "## Statistics",
    sprintf("Spearman rank correlation: robust to outliers and skewed values, no transformation needed. Adjusted: partial Spearman correlation for %s (all variables ranked, covariates regressed out). Proteome-wide: Benjamini-Hochberg FDR over the %d proteins of each parameter, significant = FDR < %g. Galanin is pre-specified: its single-test p-values are primary, the FDR is shown alongside.",
            adj_txt, ncol(Y), fdr_cut),
    sprintf("Power: with n = %d, a correlation needs about |rho| >= %.2f for p < 0.05 and about |rho| >= %.2f to pass the proteome-wide FDR. Weaker true correlations are missed; 'not significant' does not mean 'no correlation'.",
            n_med, r_crit(0.05, n_med), r_crit(0.05 / max(ncol(Y), 1), n_med))),
    subtitle = "Olink Explore HT, PC-normalised NPX; clinical data from the LEIP clinical file")
})
section("Answers", {
  items <- character()
  for (qq in unique(answers$question)) {
    a <- answers |> filter(question == qq)
    items <- c(items, paste("##", qq), sprintf("%s: %s. %s", a$item, a$verdict, a$evidence))
  }
  page_text("Answers", items, subtitle = "automatically derived from the tests; evidence in leip_biobank.xlsx and galanin/galanin.xlsx", size = 9.5)
})
section("Galanin: Olink vs ELISA", {
  if (!has_gal || !has_elisa) stop("Olink GAL or the galanin ELISA is not available")
  page_plot(figs$gal_elisa)
  page_table("Galanin: agreement of Olink GAL with the ELISA",
             agree |> select(analysis, n, rho, ci_low, ci_high, p, npx_per_doubling),
             note = "rho = Spearman; npx_per_doubling = change of Olink NPX per doubling of the ELISA value (1 = same fold-change)")
  page_table("Galanin tertiles: ELISA (rows) vs Olink (columns)", tert_tab,
             note = sprintf("same tertile: %.0f%% (chance 33%%); weighted kappa %.2f", tert_stats$same_tertile_pct, tert_stats$weighted_kappa))
  if (nrow(plate_fx)) page_table("Galanin and plates (Kruskal-Wallis)", plate_fx |> mutate(medians = str_trunc(medians, 38)),
                                 note = "full medians: galanin/galanin.xlsx, sheet plate_effects")
})
section("Lab vs Olink", {
  if (is.null(figs$bench)) stop("no protein is measured by both the lab and Olink")
  page_plot(figs$bench)
})
section("Galanin vs clinical parameters", {
  if (is.null(figs$gal_clin_short)) stop("no galanin measurement available")
  page_plot(figs$gal_clin_short)
  page_table("Galanin vs clinical parameters (p < 0.05)", gcl |> filter(p < 0.05) |> arrange(measure, p) |>
               select(measure, parameter = label, n, rho, p, fdr, rho_adj, p_adj),
             note = sprintf("adjusted: Olink GAL for %s, ELISA for %s. All parameters: galanin/galanin_vs_clinical.csv", adj_txt, adj_txt_elisa))
})
section("Galanin vs proteins", {
  if (!has_gal) stop("Olink GAL is not available")
  page_plot(figs$gal_volcano)
  page_table("Proteins most strongly correlated with Olink GAL", gp |> head(30) |> select(protein, n, rho, p, fdr, rho_adj, p_adj),
             note = "all proteins: galanin/GAL_vs_proteins.csv")
  if (!is.null(gsea) && nrow(gsea)) page_table("Gene sets (GSEA on the GAL correlation ranking)",
                                               gsea |> head(25) |> transmute(pathway = str_trunc(pathway, 60), size, NES, pval, padj))
})
section("Proteome vs clinical parameters", {
  if (!is.null(figs$nsig)) page_plot(figs$nsig)
  if (!is.null(figs$heatmap)) page_plot(figs$heatmap)
  page_table("Proteins per clinical parameter", nsig |> head(40) |>
               select(parameter = label, n = n_samples, `FDR sig.` = significant, positive, negative, adjusted = significant_adjusted,
                      `p < 0.05` = p_below_05, `by chance` = expected_by_chance),
             note = "by chance = number of proteins expected at p < 0.05 without any true correlation; all parameters: n_significant_per_parameter.csv")
})
section("Sanity checks", {
  if (!nrow(sanity)) stop("no expected associations configured (leip_biobank: expected_associations)")
  page_table("Sanity checks: associations known from population studies",
             sanity |> select(protein, parameter = label, expected, n, rho, p, status),
             note = "recovered = p < 0.05 in the expected direction; with n = 34 weaker true associations can be missed")
})
section("Methods", page_text("Methods and caveats", c(
  "## Data",
  "LEIP biobank serum samples only (population controls), Olink Explore HT, PC-normalised NPX from step 02 (samples failing Olink QC removed). Clinical data from the LEIP clinical file: sheet Key_parameters plus all further variables of sheet All_SORB_parameters, linked by Olink sample ID and checked against the manifest (SubjectID, Olink plate).",
  "Variables left out: identifiers, log copies (same ranks as the raw value), variables without variation or with too few values, coarsened copies listed in config.yml, and variables with the same ranks as another one (e.g. glucose in mmol/l and mg/dl). Column names are converted to plain ASCII (the micro sign becomes u, umlauts lose their dots).",
  "## Statistics",
  "Spearman rank correlation; p-value from the t approximation (as cor.test with exact = FALSE); 95% CI by Fisher z (Bonett-Wright variance). Partial Spearman: all variables ranked, covariates regressed out, Pearson correlation of the residuals (as in ppcor). Binary variables (e.g. sex) are coded 0/1; for them the median NPX difference is given as well.",
  "Olink vs ELISA: Spearman correlation, also adjusted for Olink plate, within each plate, only above LOD and without flagged samples; slope of NPX on log2(ELISA) (1 = the same fold-change); agreement of tertiles with linear-weighted kappa. Benchmark: the same comparison for every protein measured by the lab and by Olink (config: lab_vs_olink), with the rank of the matching Olink assay among all proteins.",
  "## Caveats",
  "n = 34: exploratory; only strong correlations reach the proteome-wide FDR. Correlation is not causation, and many clinical parameters are correlated with each other (e.g. BMI, waist, insulin), so their protein lists overlap.",
  "NPX is relative: Olink and ELISA can only agree in ranking, not in absolute values. Olink GAL and the ELISA may detect different forms of galanin (precursor vs mature peptide), and values below the Olink LOD are mostly noise.",
  "The LEIP samples were collected and stored differently from the study samples (see the main pipeline); this analysis compares LEIP samples only with each other.",
  "## Files (output/leip_biobank/)",
  "answers.csv - the answers above | leip_biobank.xlsx - parameters, samples, detection, significant associations, top 10 per parameter, sanity checks | associations_all.csv.gz - every protein x parameter | galanin/ - galanin.xlsx and figures")))
invisible(grDevices::dev.off())
msg("Done: %s", pdf_file)
for (i in seq_len(nrow(answers))) msg("%s | %s: %s", str_extract(answers$question[i], "^[A-Z][0-9]"), answers$item[i], answers$verdict[i])
