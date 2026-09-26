# Per-protein differential models.
# Formulas with a random effect, e.g. (1|SubjectID), are fitted with variancePartition::dream
# (limma with a linear mixed model and empirical Bayes moderation); formulas without one use limma.

# lme4 >= 1.1-36 moved nobars() to the reformulas package
nobars <- function(f) if (requireNamespace("reformulas", quietly = TRUE)) reformulas::nobars(f) else lme4::nobars(f)

has_random <- function(form) any(grepl("|", deparse(form), fixed = TRUE))

#' Fit one model and return tidy results for the requested contrasts.
#' @param expr   assays x samples matrix (rownames = OlinkID)
#' @param info   sample table; must contain column SampleID matching colnames(expr)
#' @param form   model formula (right-hand side only)
#' @param contrasts named character vector, e.g. c(L_vs_NL = "condAD_L - condAD_NL");
#'        for a numeric covariate give the coefficient name, e.g. c(per_week = "weeks")
#' @param model  label written into the results
fit_contrasts <- function(expr, info, form, contrasts, model, min_group_n = 3) {
  info <- as.data.frame(info)
  vars <- all.vars(form)
  info <- info[stats::complete.cases(info[, vars, drop = FALSE]), , drop = FALSE]
  info <- info[info$SampleID %in% colnames(expr), , drop = FALSE]
  for (v in vars) if (is.character(info[[v]])) info[[v]] <- factor(info[[v]])
  info <- droplevels(info)
  expr <- expr[, info$SampleID, drop = FALSE]
  rownames(info) <- info$SampleID
  if (ncol(expr) < 4) { msg("  %s: skipped (only %d samples)", model, ncol(expr)); return(NULL) }

  # drop constant factors (e.g. a single plate in a subset)
  fixed_terms <- attr(terms(nobars(form)), "term.labels")
  single <- fixed_terms[vapply(fixed_terms, \(v) v %in% names(info) && is.factor(info[[v]]) && nlevels(info[[v]]) < 2, logical(1))]
  if (length(single)) {
    form <- update(form, as.formula(paste("~ . -", paste(single, collapse = " - "))))
    msg("  %s: dropped constant term(s) %s", model, paste(single, collapse = ", "))
  }

  X <- model.matrix(nobars(form), info)
  n_per_col <- colSums(X != 0)
  ok <- vapply(contrasts, \(ct) {
    cols <- all.vars(parse(text = ct)[[1]])
    all(cols %in% colnames(X)) && all(n_per_col[cols] >= min_group_n)
  }, logical(1))
  if (!all(ok)) {
    msg("  %s: contrast(s) not estimable or group < %d samples, skipped: %s",
        model, min_group_n, paste(names(contrasts)[!ok], collapse = ", "))
    contrasts <- contrasts[ok]
  }
  if (!length(contrasts)) return(NULL)

  # assays with too many missing values in this subset cannot be modelled
  expr <- expr[rowMeans(is.na(expr)) <= 0.2, , drop = FALSE]
  if (anyNA(expr)) expr <- t(apply(expr, 1, \(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x }))

  if (has_random(form)) {
    L <- suppressWarnings(variancePartition::makeContrastsDream(form, info, contrasts = contrasts))
    fit <- suppressMessages(suppressWarnings(
      variancePartition::dream(expr, form, info, L, BPPARAM = bpparam_cores(), quiet = TRUE)))
    fit <- variancePartition::eBayes(fit)
    top <- \(ct) variancePartition::topTable(fit, coef = ct, number = Inf, sort.by = "none")
  } else {
    design <- model.matrix(form, info)
    colnames(design) <- make.names(colnames(design))   # e.g. "platePlate 3" -> "platePlate.3"
    cm <- limma::makeContrasts(contrasts = contrasts, levels = design)
    colnames(cm) <- names(contrasts)
    fit <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr, design), cm))
    top <- \(ct) limma::topTable(fit, coef = ct, number = Inf, sort.by = "none")
  }

  map(names(contrasts), \(ct) {
    tt <- top(ct)
    tibble(model = model, contrast = ct, OlinkID = rownames(tt), logFC = tt$logFC,
           AveExpr = tt$AveExpr, t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
           n_samples = ncol(expr), n_subjects = n_distinct(info$SubjectID),
           formula = paste(deparse(form), collapse = ""))
  }) |> bind_rows()
}

#' Add assay names and a significance flag.
annotate_results <- function(res, assay_map, fdr) {
  res |>
    left_join(assay_map, by = "OlinkID") |>
    relocate(Assay, .after = OlinkID) |>
    mutate(significant = adj.P.Val < fdr)
}

volcano <- function(res, title) {
  top <- res |> group_by(contrast) |> slice_min(P.Value, n = 8, with_ties = FALSE)
  ggplot(res, aes(logFC, -log10(P.Value), colour = significant)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text(data = top, aes(label = Assay), size = 2.6, vjust = -0.6, colour = "black", check_overlap = TRUE) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) +
    facet_wrap(~contrast, scales = "free") +
    labs(title = title, x = "difference in NPX (log2)", y = "-log10 p") +
    theme(legend.position = "bottom")
}

#' Run a list of model specifications and write results + volcano plots.
#' Each spec: list(name, samples (character SampleIDs), formula, contrasts)
run_model_specs <- function(specs, expr, info, cfg, prefix, assay_map) {
  res <- map(specs, \(sp) {
    msg("Model %s: %d samples", sp$name, length(sp$samples))
    fit_contrasts(expr, info |> filter(SampleID %in% sp$samples), sp$formula, sp$contrasts,
                  sp$name, cfg$stats$min_group_n)
  }) |> bind_rows()
  if (!nrow(res)) return(invisible(res))
  res <- annotate_results(res, assay_map, cfg$stats$fdr)
  save_csv(res, cfg, "models", paste0(prefix, "_results.csv"))
  summ <- res |> group_by(model, contrast, n_samples, n_subjects) |>
    summarise(n_assays = n(), n_sig = sum(significant), n_up = sum(significant & logFC > 0),
              n_down = sum(significant & logFC < 0), .groups = "drop")
  save_csv(summ, cfg, "models", paste0(prefix, "_summary.csv"))
  print(as.data.frame(summ))
  for (mdl in unique(res$model)) {
    r <- res |> filter(model == mdl)
    save_plot(volcano(r, paste(prefix, mdl)), cfg, "models", "volcano", paste0(prefix, "_", mdl, ".png"),
              width = 4 + 3.5 * min(3, n_distinct(r$contrast)), height = 4 + 3.5 * (n_distinct(r$contrast) > 3))
  }
  invisible(res)
}

#' Pre-specified test for ONE protein (focus proteins, step 12): same formula and contrasts as the
#' proteome-wide models, fitted with lmerTest (Satterthwaite df) or lm, without empirical Bayes
#' moderation. Returns one row per contrast; p is the unadjusted single-protein p-value.
#' @param df one row per sample with columns `value` and the model variables
test_single <- function(df, form, contrasts, model, min_group_n = 3) {
  df <- as.data.frame(df)
  vars <- c("value", all.vars(form))
  df <- df[stats::complete.cases(df[, intersect(vars, names(df)), drop = FALSE]), , drop = FALSE]
  for (v in all.vars(form)) if (is.character(df[[v]])) df[[v]] <- factor(df[[v]])
  df <- droplevels(df)
  fixed_terms <- attr(terms(nobars(form)), "term.labels")
  single <- fixed_terms[vapply(fixed_terms, \(v) v %in% names(df) && is.factor(df[[v]]) && nlevels(df[[v]]) < 2, logical(1))]
  if (length(single)) form <- update(form, as.formula(paste("~ . -", paste(single, collapse = " - "))))
  if (nrow(df) < 4) return(NULL)
  X <- model.matrix(nobars(form), df)
  full <- as.formula(paste("value", paste(deparse(form), collapse = "")))
  random <- has_random(form) && n_distinct(df$SubjectID) < nrow(df)
  fit <- tryCatch(
    if (random) suppressMessages(lmerTest::lmer(full, data = df)) else lm(nobars(full), data = df),
    error = \(e) NULL)
  if (is.null(fit)) return(NULL)
  b <- if (random) lme4::fixef(fit) else coef(fit)
  n_obs <- nrow(df); n_subj <- n_distinct(df$SubjectID)
  map(names(contrasts), \(nm) {
    ct <- contrasts[[nm]]
    cols <- all.vars(parse(text = ct)[[1]])
    if (!all(cols %in% names(b)) || any(is.na(b[cols])) || any(colSums(X[, cols, drop = FALSE] != 0) < min_group_n))
      return(NULL)
    # contrast vector: evaluate the (linear) contrast expression with unit vectors
    L <- vapply(names(b), \(cl) eval(parse(text = ct), envir = as.list(setNames(as.numeric(names(b) == cl), names(b)))),
                numeric(1))
    L[is.na(b)] <- 0
    if (random) {
      r <- lmerTest::contest1D(fit, L)
      est <- r$Estimate; se <- r$`Std. Error`; dfree <- r$df; p <- r$`Pr(>|t|)`
    } else {
      bb <- b; bb[is.na(bb)] <- 0
      V <- vcov(fit); V[is.na(V)] <- 0
      est <- sum(L * bb); se <- sqrt(drop(t(L) %*% V %*% L)); dfree <- df.residual(fit)
      p <- 2 * pt(-abs(est / se), dfree)
    }
    tibble(model = model, contrast = nm, estimate = est, se = se, df_resid = dfree, t = est / se, p = p,
           ci_low = est - qt(0.975, dfree) * se, ci_high = est + qt(0.975, dfree) * se,
           n_samples = n_obs, n_subjects = n_subj,
           method = if (random) "lmer (Satterthwaite)" else "lm")
  }) |> bind_rows()
}
