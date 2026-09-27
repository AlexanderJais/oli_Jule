# Helpers for the LEIP biobank analysis (step 18): reading the clinical file, choosing the
# clinical parameters, and fast (partial) Spearman correlations of one variable with many proteins.

#' Plain-ASCII names (micro sign -> u, umlauts -> a/o/u, sharp s -> ss): the same on every system.
ascii_names <- function(x) {
  x <- stringi::stri_replace_all_regex(x, "[µμ]", "u")
  stringi::stri_trans_general(x, "Latin-ASCII")
}

#' LEIP clinical file: sheet Key_parameters, plus every further variable of sheet All_SORB_parameters
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
             parameter %in% exclude ~ "excluded in config.yml (leip_biobank: exclude_parameters)",
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
