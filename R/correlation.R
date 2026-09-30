# Helpers for correlation analyses (steps 18-20): assay lookup, Spearman / repeated-measures
# correlation with a fixed output format, below-LOD handling, and a limma model with subject blocking.

#' OlinkID of an assay given by gene symbol, OlinkID or UniProt ("TPSAB1" also finds "TPSAB1_TPSB2").
find_assay <- function(assay_map, name) {
  N <- toupper(name)
  hit <- assay_map |>
    filter(toupper(Assay) == N | str_detect(toupper(Assay), paste0("(^|_)", N, "($|_)")) | OlinkID == name |
             (if ("UniProt" %in% names(assay_map)) coalesce(UniProt == name, FALSE) else FALSE))
  if (nrow(hit)) hit$OlinkID[1] else NA_character_
}

#' Spearman correlation as a one-row tibble (n, rho, p); NA below 5 complete pairs.
spearman_row <- function(x, y, min_n = 5) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < min_n || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2)
    return(tibble(n = sum(ok), rho = NA_real_, p = NA_real_))
  s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = sum(ok), rho = unname(s$estimate), p = s$p.value)
}

#' Repeated-measures correlation (Bakdash & Marusich 2017; = correlation of subject-centred values,
#' i.e. the within-subject relationship). Uses subjects with >= 2 complete observations.
rmcorr_row <- function(x, y, subj, min_n = 5) {
  ok <- !is.na(x) & !is.na(y)
  multi <- ok & subj %in% names(which(table(subj[ok]) >= 2))
  empty <- tibble(n_within = sum(multi), n_subjects_within = n_distinct(subj[multi]), r_within = NA_real_, p_within = NA_real_)
  if (sum(multi) < min_n || n_distinct(subj[multi]) < 2) return(empty)
  r <- tryCatch(suppressWarnings(rmcorr::rmcorr(participant = subj, measure1 = x, measure2 = y,
                                                dataset = data.frame(subj = factor(subj[multi]), x = x[multi], y = y[multi]))),
                error = \(e) NULL)
  if (is.null(r)) return(empty)
  empty |> mutate(r_within = r$r, p_within = r$p)
}

#' Values below LOD set to the LOD (censoring): for rank correlation all of them tie at the bottom.
censor_at_lod <- function(value, lod) if_else(!is.na(lod) & value < lod, lod, value)

#' Drop aliased design columns (e.g. plate levels that coincide with a compartment).
estimable_design <- function(X) {
  q <- qr(X)
  keep <- sort(q$pivot[seq_len(q$rank)])
  X[, keep, drop = FALSE]
}

#' limma with subject blocking (duplicateCorrelation) and optional variance weights per group of samples
#' (arrayWeights, e.g. per compartment when dISF and serum are stacked in one model).
#' @param contrasts named character vector of contrasts on the design columns (make.names applied)
#' @param ftests named list of character vectors: contrast names tested jointly (moderated F)
#' @return list(contrasts = tidy per-contrast results, ftests = tidy F results, info)
limma_block <- function(expr, info, form, contrasts, ftests = list(), block = "SubjectID", weights_by = NULL,
                        min_group_n = 3, label = "") {
  info <- as.data.frame(info)
  vars <- unique(c(all.vars(form), block, weights_by))
  info <- info[stats::complete.cases(info[, vars, drop = FALSE]) & info$SampleID %in% colnames(expr), , drop = FALSE]
  for (v in all.vars(form)) if (is.character(info[[v]])) info[[v]] <- factor(info[[v]])
  info <- droplevels(info)
  if (nrow(info) < 6) { msg("  %s: skipped (only %d samples)", label, nrow(info)); return(NULL) }
  single <- all.vars(form)[vapply(all.vars(form), \(v) is.factor(info[[v]]) && nlevels(info[[v]]) < 2, logical(1))]
  if (length(single)) form <- update(form, as.formula(paste("~ . -", paste(single, collapse = " - "))))
  X <- estimable_design(model.matrix(form, info))
  colnames(X) <- make.names(colnames(X))
  n_col <- colSums(X != 0)
  ok <- vapply(contrasts, \(ct) {
    cols <- all.vars(parse(text = ct)[[1]])
    all(cols %in% colnames(X)) && all(n_col[cols] >= min_group_n)
  }, logical(1))
  contrasts <- contrasts[ok]
  if (!length(contrasts)) { msg("  %s: no estimable contrast", label); return(NULL) }
  ex <- expr[, info$SampleID, drop = FALSE]
  ex <- ex[rowMeans(is.na(ex)) <= 0.2, , drop = FALSE]
  if (anyNA(ex)) ex <- t(apply(ex, 1, \(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))
  w <- if (!is.null(weights_by) && n_distinct(info[[weights_by]]) > 1)
    limma::arrayWeights(ex, X, var.group = info[[weights_by]]) else NULL
  bl <- info[[block]]
  cor <- if (any(duplicated(bl))) limma::duplicateCorrelation(ex, X, block = bl, weights = w)$consensus.correlation else NA_real_
  fit <- if (is.na(cor)) limma::lmFit(ex, X, weights = w) else limma::lmFit(ex, X, block = bl, correlation = cor, weights = w)
  cm <- limma::makeContrasts(contrasts = contrasts, levels = X)
  colnames(cm) <- names(contrasts)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  meta_cols <- tibble(model = label, n_samples = nrow(info), n_subjects = n_distinct(bl), block_correlation = cor)
  res <- map(names(contrasts), \(ct) {
    tt <- limma::topTable(fit2, coef = ct, number = Inf, sort.by = "none")
    tibble(OlinkID = rownames(tt), contrast = ct, logFC = tt$logFC, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
  }) |> bind_rows() |> mutate(!!!as.list(meta_cols))
  fres <- imap(ftests, \(cts, nm) {
    cts <- intersect(cts, names(contrasts))
    if (length(cts) < 1) return(NULL)
    tt <- limma::topTable(fit2, coef = cts, number = Inf, sort.by = "none")
    tibble(OlinkID = rownames(tt), ftest = nm, n_contrasts = length(cts),
           stat = if (length(cts) > 1) tt$F else tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val)
  }) |> bind_rows()
  if (nrow(fres)) fres <- fres |> mutate(!!!as.list(meta_cols))
  list(contrasts = res, ftests = fres, info = as_tibble(info))
}

#' Major-axis (orthogonal) slope of y on x: symmetric in x and y, so it is not attenuated by noise in
#' x the way an OLS slope is. Used to ask whether serum effects are a scaled-down copy of dISF effects.
ma_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  v <- stats::cov(cbind(x[ok], y[ok]))
  e <- eigen(v, symmetric = TRUE)$vectors[, 1]
  e[2] / e[1]
}
