# Helpers of the LEIP galanin study: configuration, reading the Olink serum file and the LEIP
# clinical file, choosing the clinical variables, and fast (partial) Spearman correlations.
# Uses the repository's generic helpers: R/utils.R (config, output paths, logging), R/qc.R (Olink LOD).

leip_config <- function() load_config(Sys.getenv("LEIP_CONFIG", "leip_galanin/config.yml"))

fmt_p   <- \(p) ifelse(is.na(p), "n/a", sprintf("%.2g", p))
col_or  <- \(d, nm, default = NA) if (nm %in% names(d)) d[[nm]] else rep(default, nrow(d))
ensure  <- \(d, cols) { for (c in setdiff(cols, names(d))) d[[c]] <- rep(NA_real_, nrow(d)); d }
zscore  <- \(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
r_crit  <- \(p, n) { t <- qt(1 - p / 2, n - 2); t / sqrt(t^2 + n - 2) }      # |rho| needed for p at n samples
top_str <- \(name, rho, p, k = 5) {
  o <- head(order(p), k); o <- o[!is.na(p[o])]
  if (!length(o)) "none" else paste(sprintf("%s (%+.2f)", name[o], rho[o]), collapse = ", ")
}
#' OlinkIDs of an assay given as gene symbol (also inside combined names such as "IL12A_IL12B").
find_assay <- \(det, a) det$OlinkID[str_detect(toupper(det$Assay), paste0("(^|_)", toupper(a), "($|_)"))]

#' The data of step 01, plus: M = measurable proteins (samples x proteins), gal = OlinkID of galanin,
#' npx = its values, elisa = the ELISA values, cov = covariates of the adjusted analyses.
leip_load <- function(cfg) {
  d <- read_step(cfg, "0_data", "leip_data.rds", step = "leip_galanin/scripts/01_data.R")
  d$M <- d$Y[, d$det$OlinkID[d$det$measurable], drop = FALSE]
  g <- find_assay(d$det, cfg$galanin$olink_assay %||% "GAL")
  d$gal <- if (length(g)) g[1] else NA_character_
  d$npx <- if (!is.na(d$gal)) unname(d$Y[, d$gal]) else rep(NA_real_, nrow(d$Y))
  d$elisa <- col_or(d$S, "galanin_elisa", NA_real_)
  d$has_gal <- !is.na(d$gal); d$has_elisa <- sum(!is.na(d$elisa)) >= 10
  d$cov <- intersect(unlist(cfg$covariates %||% c("age", "sex_male", "plate")), names(d$S))
  d
}
covs_label <- \(v) if (length(v)) paste(recode(v, sex_male = "sex", plate = "Olink plate"), collapse = ", ") else "nothing"

#' Figures of one analysis step: fig() saves the PNG and keeps the plot, figs_save() stores all plots
#' for the PDF report (plots must only use columns of their own data, no outside variables).
figs_new  <- function() { e <- new.env(); e$list <- list(); e }
fig       <- function(f, name, p, cfg, ..., width = 8, height = 6) {
  save_plot(p, cfg, ..., width = width, height = height); f$list[[name]] <- p; invisible(p)
}
figs_save <- function(f, cfg, ...) saveRDS(f$list, out_path(cfg, ...))

#' Answers of one analysis step: add rows with answer(), write them with answers_save().
answers_new <- function() { e <- new.env(); e$rows <- list(); e }
answer <- function(a, question, item, verdict, evidence)
  a$rows[[length(a$rows) + 1]] <- tibble(question = question, item = item, verdict = verdict, evidence = evidence)
answers_save <- function(a, cfg, ...) { d <- bind_rows(a$rows); save_csv(d, cfg, ...); d }

# ---- Olink data ------------------------------------------------------------------------------------------------

#' Olink rows of the LEIP samples: reads every NPX parquet file that contains one of `ids`, adds the
#' LOD (R/qc.R add_lod: Olink fixed LOD file, per sample for count-based assays), and returns the
#' sample x assay rows of these samples. Control IDs are made unique per plate, as in step 02 of
#' the dISF/serum pipeline.
read_leip_npx <- function(npx_dir, ids, fixed_lod = NULL, value_col = "PCNormalizedNPX") {
  files <- list.files(npx_dir, pattern = "\\.parquet$", full.names = TRUE)
  if (!length(files)) stop("No .parquet files in ", npx_dir)
  d <- map(files, \(f) {
    x <- OlinkAnalyze::read_npx(f) |> as_tibble()
    if (!any(x$SampleID %in% ids)) { msg("%s: no LEIP sample - skipped", basename(f)); return(NULL) }
    msg("%s: %d LEIP samples", basename(f), n_distinct(x$SampleID[x$SampleID %in% ids]))
    req <- c("SampleID", "SampleType", "PlateID", "OlinkID", "Assay", "AssayType", value_col, "SampleQC", "DataAnalysisRefID")
    miss <- setdiff(req, names(x))
    if (length(miss)) stop(basename(f), " is missing columns: ", paste(miss, collapse = ", "))
    x <- x |> mutate(SampleID = if_else(SampleType == "SAMPLE", SampleID, paste(SampleID, PlateID, sep = "@")))
    add_lod(x, fixed_lod, "auto", value_col = value_col) |> mutate(source_file = basename(f))
  }) |> compact()
  if (!length(d)) stop("None of the LEIP samples of the clinical file is in the NPX files in ", npx_dir)
  bind_rows(d) |>
    filter(SampleType == "SAMPLE", AssayType == "assay", SampleID %in% ids) |>
    mutate(value = .data[[value_col]], below_lod = if_else(is.na(LOD), NA, value < LOD))
}

# ---- clinical data -----------------------------------------------------------------------------------------------

#' Plain-ASCII names (micro sign -> u, umlauts -> a/o/u, sharp s -> ss): the same on every system.
ascii_names <- function(x) {
  x <- stringi::stri_replace_all_regex(x, paste0("[", intToUtf8(c(0xb5, 0x3bc)), "]"), "u")
  stringi::stri_trans_general(x, "Latin-ASCII")
}

#' LEIP clinical file: sheet Key_parameters plus every further variable of sheet All_SORB_parameters
#' (if present). One row per Olink sample. Column names become ASCII; the galanin ELISA column is
#' renamed to galanin_elisa and sex_MF becomes sex_male (1 = M, 0 = F).
#' Returns list(data, key = names of the Key_parameters columns).
read_leip_clinical <- function(path, galanin_col = "Galanin [pg/mL]") {
  sheets <- readxl::excel_sheets(path)
  if (!"Key_parameters" %in% sheets) stop("Sheet 'Key_parameters' not found in ", path)
  gal <- ascii_names(galanin_col)
  prep <- function(d) {
    names(d) <- ascii_names(names(d))
    if (!"Olink_SampleID" %in% names(d)) stop("Column 'Olink_SampleID' not found in ", path)
    d <- d |> rename(SampleID = Olink_SampleID) |> filter(!is.na(SampleID))
    if (gal %in% names(d)) d <- d |> rename(galanin_elisa = all_of(gal))
    if ("sex_MF" %in% names(d))
      d <- d |> mutate(sex_male = case_when(toupper(sex_MF) == "M" ~ 1, toupper(sex_MF) == "F" ~ 0), .after = sex_MF)
    d |> select(-any_of(c("sex_MF", "sex")))            # sex is kept once, as sex_male
  }
  key <- prep(readxl::read_excel(path, sheet = "Key_parameters"))
  key_cols <- setdiff(names(key), "SampleID")
  if ("All_SORB_parameters" %in% sheets) {
    all <- prep(readxl::read_excel(path, sheet = "All_SORB_parameters"))
    key <- key |> left_join(all |> select(SampleID, all_of(setdiff(names(all), names(key)))), by = "SampleID")
  }
  list(data = key, key = key_cols)
}

# identifiers and technical columns of the clinical file (never correlated)
leip_technical <- c("SubjectID", "SORB_barcode", "Olink_plate", "Olink_well", "Galanin_ELISA_plate",
                    "Galanin_ELISA_well", "Ifd.Nr_Uli")

# readable names for the SORB variables that are unambiguous; all others keep their SORB name
leip_labels <- c(
  age = "age", sex_male = "sex (male = 1)", BMI = "BMI", WHR = "waist-hip ratio", WtHR = "waist-height ratio",
  c_fett = "body fat %", c_gewich = "body weight", c_groess = "height", c_bauch = "waist circumference",
  c_huefte = "hip circumference", c_mager = "lean mass", c_water = "body water",
  Gluc0_mg_dl = "fasting glucose", Gluc30_mg_dl = "glucose 30 min OGTT", Gluc120_mg_dl = "glucose 120 min OGTT",
  GLUK_30 = "glucose 30 min OGTT", GLUK_120 = "glucose 120 min OGTT",
  Ins0_uU_ml = "fasting insulin", Ins30_uU_ml = "insulin 30 min OGTT", Ins120_uU_ml = "insulin 120 min OGTT",
  ins30 = "insulin 30 min OGTT", ins120 = "insulin 120 min OGTT",
  HOMA_IR = "HOMA-IR", HOMA_B = "HOMA-B", c_CRP = "CRP (lab)", C_CHOL = "total cholesterol", C_HDL = "HDL cholesterol",
  C_LDL = "LDL cholesterol", C_TRIGLY = "triglycerides", c_apo = "ApoA-I (lab)", C_APO_B = "ApoB (lab)",
  C_LIPO = "Lp(a) (lab)", MDRD_kurz = "eGFR (MDRD)", galanin_elisa = "galanin ELISA",
  c_adiponectin = "adiponectin (lab)", IL10 = "IL-10 (lab)", serum_irisin = "irisin (lab)", Progranulin = "progranulin (lab)",
  Vaspin0min = "vaspin fasting (lab)", c_chemerin = "chemerin (lab)", BMP2 = "BMP2 (lab)", c_FGF21 = "FGF21 (lab)",
  c_IGF1 = "IGF-1 (lab)", c_AFABP4 = "A-FABP / FABP4 (lab)", c_AGF = "AGF / ANGPTL6 (lab)", c_alat = "ALT", c_asat = "AST",
  c_ggt = "GGT", c_bili = "bilirubin", c_tsh = "TSH", c_FT3 = "free T3", c_FT4 = "free T4", Harnsaure_im_Serum = "uric acid",
  C_KREAT = "creatinine", HST = "urea", C_G_EIW = "total protein", SYS_MW = "systolic blood pressure",
  c_DIA_MW = "diastolic blood pressure", IMD_L_MW = "intima-media thickness left", IMD_R_MW = "intima-media thickness right",
  RESTRAINT = "eating restraint (TFEQ)", c_DISINHIB = "disinhibition (TFEQ)", HUNGER = "hunger (TFEQ)",
  cfDNA = "cell-free DNA")

param_label <- function(x) unname(if_else(x %in% names(leip_labels), leip_labels[x], x))

#' Which clinical variables are analysed? One row per column of `clin`; `reason` says why a
#' variable is left out (NA = analysed). Rules, in this order: identifiers / technical, excluded in
#' the config, not numeric, log copies (ln_/lg_: same ranks as the raw value), no variation, fewer
#' than min_n values, fewer than min_group samples outside the most common value, and variables
#' whose ranks are (nearly) identical to an earlier variable (e.g. the same value in other units).
select_parameters <- function(clin, key_cols, exclude = character(), min_n = 20, min_group = 5, dup_rho = 0.99) {
  cols <- setdiff(names(clin), "SampleID")
  cols <- c(intersect(key_cols, cols), setdiff(cols, key_cols))      # key parameters first: kept among duplicates
  info <- tibble(
    parameter = cols, label = param_label(cols), key = cols %in% key_cols,
    numeric = map_lgl(cols, \(c) is.numeric(clin[[c]]) || is.logical(clin[[c]])),
    n = map_int(cols, \(c) sum(!is.na(clin[[c]]))),
    n_values = map_int(cols, \(c) n_distinct(na.omit(clin[[c]]))),
    n_outside_mode = map_int(cols, \(c) { t <- table(clin[[c]]); if (length(t)) sum(t) - max(t) else 0L })
  ) |>
    mutate(type = case_when(!numeric ~ "text", n_values == 2 ~ "binary", n_values <= 4 ~ "ordinal", TRUE ~ "continuous"),
           reason = case_when(
             parameter %in% leip_technical ~ "identifier / technical",
             parameter %in% exclude ~ "excluded in the config (clinical: exclude)",
             !numeric ~ "not numeric",
             str_detect(parameter, "^(ln|lg)_|_lg$|_lg_") ~ "log copy of another variable (same ranks)",
             n_values < 2 ~ "no variation",
             n < min_n ~ sprintf("too few values (< %d)", min_n),
             n_outside_mode < min_group ~ sprintf("fewer than %d samples differ from the most common value", min_group),
             TRUE ~ NA_character_))
  # rank duplicates: e.g. glucose in mmol/l and mg/dl, HOMA-IS = 1 / HOMA-IR, the same AUC in two units
  cand <- info$parameter[is.na(info$reason)]
  if (length(cand) > 1) {
    r <- suppressWarnings(cor(as.matrix(clin[cand]), method = "spearman", use = "pairwise.complete.obs"))
    shared <- crossprod(!is.na(as.matrix(clin[cand])))
    kept <- character()
    for (c in cand) {
      dup <- kept[abs(r[c, kept]) >= dup_rho & shared[c, kept] >= min_n]
      dup <- dup[!is.na(dup)]
      if (length(dup)) info$reason[info$parameter == c] <- sprintf("same ranks as %s (rho = %.3f)", dup[1], r[c, dup[1]])
      else kept <- c(kept, c)
    }
  }
  info |> mutate(analysed = is.na(reason)) |> select(-numeric) |>
    group_by(analysed, label) |>                                     # labels must be unique among analysed variables
    mutate(label = if (n() > 1) paste0(label, " [", parameter, "]") else label) |>
    ungroup()
}

# ---- statistics ------------------------------------------------------------------------------------------------

#' Covariate matrix for partial correlations: intercept, ranked numeric covariates, dummy-coded
#' factors. Covariates without variation in these rows are dropped.
covariate_matrix <- function(Z) {
  cols <- map(names(Z), \(nm) {
    z <- Z[[nm]]
    if (is.numeric(z) || is.logical(z)) {
      if (n_distinct(z) < 2) return(NULL)
      return(matrix(rank(z), dimnames = list(NULL, nm)))
    }
    f <- factor(z)
    if (nlevels(f) < 2) return(NULL)
    m <- stats::model.matrix(~f)[, -1, drop = FALSE]
    colnames(m) <- paste0(nm, levels(f)[-1])
    m
  }) |> compact()
  do.call(cbind, c(list(`(Intercept)` = rep(1, nrow(Z))), cols))
}

#' Spearman correlation of x with every column of Y (samples in rows). With covariates Z, a partial
#' Spearman correlation: all variables are ranked and the covariates regressed out (as in ppcor).
#' Samples with a missing value are left out per protein. p-value from the t approximation, as in
#' cor.test(method = "spearman", exact = FALSE); 95% CI by Fisher z with the Bonett-Wright variance.
spearman_vs <- function(x, Y, Z = NULL, min_n = 10) {
  Y <- as.matrix(Y)
  ok <- !is.na(x)
  if (!is.null(Z) && ncol(Z)) ok <- ok & stats::complete.cases(Z) else Z <- NULL
  pattern <- apply(is.na(Y) & ok, 2, \(v) paste(which(v), collapse = ","))    # proteins with the same missing samples
  n <- k <- rho <- rep(NA_real_, ncol(Y))
  for (pt in unique(pattern)) {
    j <- which(pattern == pt)
    rows <- setdiff(which(ok), as.integer(strsplit(pt, ",")[[1]]))
    if (length(rows) < min_n) next
    rx <- rank(x[rows])
    RY <- apply(Y[rows, j, drop = FALSE], 2, rank)
    kk <- 0
    if (!is.null(Z)) {
      q <- qr(covariate_matrix(Z[rows, , drop = FALSE]))
      kk <- q$rank - 1
      rx <- qr.resid(q, rx); RY <- qr.resid(q, RY)
    }
    rho[j] <- suppressWarnings(stats::cor(rx, RY))[1, ]
    n[j] <- length(rows); k[j] <- kk
  }
  df <- n - 2 - k
  rho[!is.na(rho) & df < 2] <- NA
  tt <- rho * sqrt(df / pmax(1 - rho^2, 1e-12))
  se <- sqrt((1 + rho^2 / 2) / (df - 1))
  tibble(id = colnames(Y), n = as.integer(n), n_covariates = as.integer(k), rho = rho,
         ci_low = tanh(atanh(rho) - 1.96 * se), ci_high = tanh(atanh(rho) + 1.96 * se),
         p = 2 * stats::pt(-abs(tt), df))
}

#' One Spearman (or partial Spearman) correlation: n, rho, 95% CI, p.
spearman1 <- function(x, y, Z = NULL, min_n = 5)
  spearman_vs(x, matrix(y, dimnames = list(NULL, "y")), Z, min_n) |> select(-id)

#' Bootstrap comparison of two correlations that share x: rho(x, a) - rho(x, b), with covariates Z
#' as partial Spearman correlations; percentile 95% CI and a two-sided bootstrap p-value.
boot_rho_diff <- function(x, a, b, Z = NULL, B = 2000) {
  dat <- bind_cols(tibble(x = x, a = a, b = b), if (!is.null(Z)) as_tibble(Z))
  dat <- dat[stats::complete.cases(dat), ]
  zc <- if (is.null(Z)) NULL else names(Z)
  rd <- function(s) {                                    # the same partial Spearman as spearman_vs(), for two pairs
    r <- cbind(rank(s$x), rank(s$a), rank(s$b))
    if (!is.null(zc)) r <- qr.resid(qr(covariate_matrix(s[zc])), r)
    suppressWarnings(cor(r[, 1], r[, 2]) - cor(r[, 1], r[, 3]))
  }
  bs <- replicate(B, rd(dat[sample.int(nrow(dat), replace = TRUE), ]))
  bs <- bs[is.finite(bs)]
  tibble(n = nrow(dat), difference = rd(dat), ci_low = unname(quantile(bs, 0.025)), ci_high = unname(quantile(bs, 0.975)),
         p_boot = max(1 / length(bs), min(1, 2 * min(mean(bs <= 0), mean(bs >= 0)))))   # never below 1 / resamples
}

#' Median NPX difference between the two values of a binary parameter (higher code minus lower).
median_diff <- function(x, Y) {
  ok <- !is.na(x)
  hi <- ok & x == max(x, na.rm = TRUE); lo <- ok & !hi
  apply(Y[hi, , drop = FALSE], 2, median, na.rm = TRUE) - apply(Y[lo, , drop = FALSE], 2, median, na.rm = TRUE)
}

#' Linear-weighted Cohen's kappa of a square table of ordered categories.
weighted_kappa <- function(tab) {
  p <- tab / sum(tab); k <- nrow(tab)
  w <- 1 - abs(outer(seq_len(k), seq_len(k), "-")) / (k - 1)
  po <- sum(w * p); pe <- sum(w * outer(rowSums(p), colSums(p)))
  (po - pe) / (1 - pe)
}

#' MSigDB gene sets (msigdbr) as a named list of gene symbols; collections like "H" or "C2:CP:REACTOME".
load_gene_sets <- function(collections) {
  map(collections, \(cl) {
    parts <- str_split_fixed(cl, ":", 2)
    g <- if (parts[2] == "") msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1])
         else msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = parts[2])
    split(g$gene_symbol, g$gs_name)
  }) |> unlist(recursive = FALSE)
}

#' GSEA (fgsea) of a protein ranking; `stat` is named by Olink assay (combined assays count per gene).
rank_gsea <- function(assay, stat, collections, min_size = 10, max_size = 500) {
  sets <- load_gene_sets(collections)
  st <- tibble(gene = str_split(assay, "_"), stat) |> filter(!is.na(stat)) |> unnest(gene) |>
    group_by(gene) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |> ungroup()
  out <- suppressWarnings(fgsea::fgseaMultilevel(sets, setNames(st$stat, st$gene), minSize = min_size, maxSize = max_size, eps = 0))
  as_tibble(out) |> mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";")) |> arrange(pval)
}
