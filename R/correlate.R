# Helpers for the correlation analysis (step 16) and the serum vs dISF signature analysis (step 17).

#' OlinkID of a protein given as gene symbol / assay name, OlinkID or UniProt (NA if not on the panel).
#' Also finds the protein inside combined assay names such as "IL12A_IL12B".
find_assay <- function(clean, name) {
  N <- toupper(name)
  a <- clean |> distinct(OlinkID, Assay, UniProt)
  hit <- a |> filter(toupper(Assay) == N | OlinkID == name | UniProt == name)
  if (!nrow(hit)) hit <- a |> filter(str_detect(toupper(Assay), paste0("(^|_)", N, "($|_)")))
  if (nrow(hit)) hit$OlinkID[1] else NA_character_
}

#' Spearman correlation with n; NA when fewer than min_n complete pairs.
spearman_n <- function(x, y, min_n = 5) {
  ok <- !is.na(x) & !is.na(y)
  n <- sum(ok)
  if (n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(tibble(n = n, rho = NA_real_, p = NA_real_))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

#' Repeated-measures correlation (within subjects); NA if too few repeated subjects.
rmcorr_n <- function(x, y, subj, min_n = 5) {
  ok <- !is.na(x) & !is.na(y) & !is.na(subj)
  ok <- ok & subj %in% names(which(table(subj[ok]) >= 2))
  if (sum(ok) < min_n || n_distinct(subj[ok]) < 2) return(tibble(n_rm = sum(ok), subjects_rm = n_distinct(subj[ok]), r_within = NA_real_, p_within = NA_real_))
  r <- tryCatch(suppressWarnings(rmcorr::rmcorr(participant = s, measure1 = x, measure2 = y,
                                                dataset = data.frame(s = factor(subj[ok]), x = x[ok], y = y[ok]))),
                error = \(e) NULL)
  tibble(n_rm = sum(ok), subjects_rm = n_distinct(subj[ok]), r_within = r$r %||% NA_real_, p_within = r$p %||% NA_real_)
}

#' Full-rank covariate design (drops constant factors and aliased columns).
covariate_design <- function(info, covars) {
  covars <- covars[vapply(covars, \(v) n_distinct(info[[v]], na.rm = TRUE) > 1, logical(1))]
  X <- if (length(covars)) model.matrix(as.formula(paste("~", paste(covars, collapse = " + "))), info) else
    matrix(1, nrow(info), 1, dimnames = list(NULL, "(Intercept)"))
  q <- qr(X)
  X[, q$pivot[seq_len(q$rank)], drop = FALSE]
}

#' Partial Spearman correlation of `anchor` with every row of `expr`, adjusted for the covariates:
#' both are rank-transformed, residualised on the covariates, then correlated.
#' p from the t distribution with n - 2 - (number of covariate columns - 1) degrees of freedom.
partial_spearman <- function(anchor, expr, info, covars, min_n = 8) {
  map(rownames(expr), \(a) {
    y <- expr[a, ]
    ok <- !is.na(y) & !is.na(anchor) & stats::complete.cases(info[, covars, drop = FALSE])
    n <- sum(ok)
    out <- tibble(OlinkID = a, n = n, rho_partial = NA_real_, p_partial = NA_real_)
    if (n < min_n) return(out)
    X <- covariate_design(info[ok, , drop = FALSE], covars)
    if (n - ncol(X) - 1 < 3) return(out)
    rx <- lm.fit(X, rank(anchor[ok]))$residuals
    ry <- lm.fit(X, rank(y[ok]))$residuals
    if (sd(rx) == 0 || sd(ry) == 0) return(out)
    r <- cor(rx, ry)
    dfree <- n - 2 - (ncol(X) - 1)
    out$rho_partial <- r
    out$p_partial <- 2 * pt(-abs(r * sqrt(dfree / max(1e-12, 1 - r^2))), dfree)
    out
  }) |> bind_rows()
}

#' limma (with duplicateCorrelation for repeated subjects) returning per-contrast effects and, for each
#' named set of contrasts, a joint F-test. Used where a joint test is needed (group x visit interaction).
#' @param form   fixed-effect formula, e.g. ~ 0 + gv + plate
#' @param block  column with the subject ID (NULL = no repeated measures)
#' @param contrasts named character vector of contrasts
#' @param f_sets named list of character vectors (names of `contrasts`) to test jointly
fit_limma_f <- function(expr, info, form, block = NULL, contrasts, f_sets = list(), min_group_n = 3) {
  info <- as.data.frame(info)
  vars <- all.vars(form)
  info <- info[stats::complete.cases(info[, vars, drop = FALSE]) & info$SampleID %in% colnames(expr), , drop = FALSE]
  for (v in vars) if (is.character(info[[v]])) info[[v]] <- factor(info[[v]])
  info <- droplevels(info)
  if (nrow(info) < 4) return(NULL)
  # drop factors that are constant or make the design rank-deficient (e.g. plate confounded with visit)
  terms_f <- attr(terms(form), "term.labels")
  keep_terms <- terms_f[vapply(terms_f, \(v) !(v %in% names(info) && is.factor(info[[v]]) && nlevels(info[[v]]) < 2), logical(1))]
  form <- as.formula(paste("~", if (attr(terms(form), "intercept") == 0) "0 +" else "", paste(keep_terms, collapse = " + ")))
  X <- model.matrix(form, info)
  if (qr(X)$rank < ncol(X) && "plate" %in% keep_terms) {
    form <- update(form, ~ . - plate); X <- model.matrix(form, info)
  }
  if (qr(X)$rank < ncol(X)) return(NULL)
  colnames(X) <- make.names(colnames(X))
  n_col <- colSums(X != 0)
  ok <- vapply(contrasts, \(ct) { cols <- all.vars(parse(text = ct)[[1]]); all(cols %in% colnames(X)) && all(n_col[cols] >= min_group_n) }, logical(1))
  contrasts <- contrasts[ok]
  if (!length(contrasts)) return(NULL)
  e <- expr[, info$SampleID, drop = FALSE]
  e <- e[rowMeans(is.na(e)) <= 0.2, , drop = FALSE]
  if (anyNA(e)) e <- t(apply(e, 1, \(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))
  cm <- limma::makeContrasts(contrasts = contrasts, levels = X); colnames(cm) <- names(contrasts)
  if (!is.null(block) && anyDuplicated(info[[block]])) {
    cf <- limma::duplicateCorrelation(e, X, block = info[[block]])
    lf <- limma::lmFit(e, X, block = info[[block]], correlation = cf$consensus.correlation)
  } else lf <- limma::lmFit(e, X)
  fit <- limma::eBayes(limma::contrasts.fit(lf, cm))
  per <- map(names(contrasts), \(ct) {
    tt <- limma::topTable(fit, coef = ct, number = Inf, sort.by = "none")
    tibble(OlinkID = rownames(tt), contrast = ct, logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
  }) |> bind_rows()
  ftests <- imap(f_sets, \(set, nm) {
    set <- intersect(set, names(contrasts))
    if (length(set) < 1) return(NULL)
    tt <- limma::topTable(fit, coef = set, number = Inf, sort.by = "none")
    tibble(OlinkID = rownames(tt), test = nm, df1 = length(set), F = if (length(set) > 1) tt$F else tt$t^2,
           P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
  }) |> bind_rows()
  list(effects = per, f = ftests, n_samples = nrow(info), n_subjects = n_distinct(info$SubjectID),
       formula = paste(deparse(form), collapse = ""))
}

#' Tissue / cell origin from the Human Protein Atlas download (proteinatlas.tsv), if available.
#' Returns one row per gene with the HPA specificity columns and a simple origin label.
hpa_annotate <- function(genes, path) {
  empty <- tibble(gene = character(), hpa_origin = character())
  h <- read_hpa(path)
  if (is.null(h)) return(empty)
  col <- \(nm) if (nm %in% names(h)) h[[nm]] else rep(NA_character_, nrow(h))
  h <- tibble(gene = h$Gene, tissue_specificity = col("RNA tissue specificity"), tissue_specific = col("RNA tissue specific nTPM"),
              cell_type_specific = col("RNA single cell type specific nCPM"), blood_cell_specific = col("RNA blood cell specific nTPM"),
              secretome = col("Secretome location")) |>
    filter(gene %in% genes) |> distinct(gene, .keep_all = TRUE)
  immune <- "T-cell|B-cell|NK-cell|[Mm]acrophage|[Mm]onocyte|[Nn]eutrophil|[Dd]endritic|cDC|pDC|[Mm]ast cell|[Pp]lasma cell|[Gg]ranulocyte|[Bb]asophil|[Ee]osinophil|Langerhans|[Ll]ymphoid|[Kk]upffer|[Hh]ofbauer|[Mm]icroglia"
  h |> mutate(
    skin = coalesce(str_detect(tissue_specific, "(^|;)\\s*skin"), FALSE),
    keratinocyte = coalesce(str_detect(cell_type_specific, "[Kk]eratinocyte"), FALSE),
    immune = coalesce(str_detect(cell_type_specific, immune), FALSE) | coalesce(str_detect(tissue_specific, "lymphoid tissue|bone marrow"), FALSE) |
      !is.na(blood_cell_specific),
    liver = coalesce(str_detect(tissue_specific, "(^|;)\\s*liver"), FALSE) | coalesce(str_detect(cell_type_specific, "[Hh]epatocyte"), FALSE),
    hpa_origin = pmap_chr(list(skin, keratinocyte, immune, liver, tissue_specificity), \(s, k, i, l, ts) {
      o <- c(if (s) "skin", if (k) "keratinocyte", if (i) "immune", if (l) "liver")
      if (length(o)) paste(o, collapse = " + ") else if (identical(ts, "Low tissue specificity")) "not tissue-specific" else "other"
    }))
}

#' Annotate a table that has an Assay column (gene symbols; combined assays use the first gene).
add_hpa <- function(df, hpa) {
  if (!nrow(hpa)) return(df |> mutate(hpa_origin = NA_character_))
  df |> mutate(.gene = str_extract(Assay, "^[^_]+")) |>
    left_join(hpa |> select(.gene = gene, hpa_origin, hpa_secretome = secretome), by = ".gene") |> select(-.gene)
}
