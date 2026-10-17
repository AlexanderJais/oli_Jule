# Helpers for the V1 relapse screen (run_relapse_v1.R): which dISF proteins are already higher at baseline in patients
# who later relapse? Built for very small groups (e.g. 5 vs 5): exact relabelling check, leave-one-out, AUC.

#' Proteins x patients matrix of one sample set (one sample per patient).
patient_matrix <- function(d, samples, tested) {
  w <- d |> filter(SampleID %in% samples$SampleID, OlinkID %in% .env$tested) |>
    group_by(OlinkID, SubjectID) |> summarise(v = mean(value, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = SubjectID, values_from = v)
  M <- as.matrix(w[, -1]); rownames(M) <- w$OlinkID; M
}

#' limma t / p / logFC for every protein (rows of M) for one labelling (TRUE = relapse).
limma_t <- function(M, is_rel) {
  design <- cbind(1, as.numeric(is_rel))
  fit <- limma::eBayes(limma::lmFit(M, design))
  tibble(OlinkID = rownames(M), logFC = fit$coefficients[, 2], t = fit$t[, 2], p = fit$p.value[, 2])
}

#' Exact relabelling check: every way of choosing n_rel of the patients as "relapse". Per protein the exact one-sided
#' (higher in relapse) and two-sided p; per labelling the number of nominal hits (p < 0.05, higher in relapse).
relabel_check <- function(M, is_rel) {
  n <- length(is_rel); combs <- utils::combn(n, sum(is_rel))
  obs <- limma_t(M, is_rel)
  ts <- matrix(NA_real_, nrow(M), ncol(combs)); hits <- integer(ncol(combs))
  up_hits <- down_hits <- matrix(FALSE, nrow(M), ncol(combs), dimnames = list(rownames(M), NULL))
  for (j in seq_len(ncol(combs))) {
    lab <- seq_len(n) %in% combs[, j]
    r <- limma_t(M, lab)
    ts[, j] <- r$t; hits[j] <- sum(r$p < 0.05 & r$logFC > 0, na.rm = TRUE)
    up_hits[, j] <- coalesce(r$p < 0.05 & r$logFC > 0, FALSE); down_hits[, j] <- coalesce(r$p < 0.05 & r$logFC < 0, FALSE)
  }
  list(per_protein = obs |> mutate(p_exact_up = rowMeans(ts >= t - 1e-12, na.rm = TRUE), p_exact_down = rowMeans(ts <= t + 1e-12, na.rm = TRUE),
                                   p_exact_two_sided = rowMeans(abs(ts) >= abs(t) - 1e-12, na.rm = TRUE)),
       hits = hits, observed_hits = sum(obs$p < 0.05 & obs$logFC > 0, na.rm = TRUE), n_labellings = ncol(combs),
       up_hits = up_hits, down_hits = down_hits, combs = combs)   # per labelling (same order for every matrix of the same patients)
}

#' Leave-one-patient-out: is the protein higher in relapse in every refit, and the largest p across refits.
leave_one_out <- function(M, is_rel) {
  fits <- map(seq_along(is_rel), \(i) limma_t(M[, -i, drop = FALSE], is_rel[-i]))
  tibble(OlinkID = rownames(M),
         loo_always_higher = Reduce(`&`, map(fits, \(f) coalesce(f$logFC > 0, FALSE))),
         loo_max_p = do.call(pmax, map(fits, "p")))
}

#' AUC per protein: share of (relapse, non-relapse) pairs in which the relapse patient has the higher value.
auc_rows <- function(M, is_rel) {
  apply(M, 1, \(v) {
    a <- v[is_rel]; b <- v[!is_rel]; a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (!length(a) || !length(b)) return(NA_real_)
    mean(outer(a, b, ">") + 0.5 * outer(a, b, "=="))
  })
}

#' One marker as a predictor, leave one patient out: the cut-off (midpoint of the two group means) and the direction
#' are set on the other patients, then the left-out patient is classified. Returns the cut-off on all patients too.
loo_predict_marker <- function(v, is_rel) {
  fit <- function(x, y) { a <- mean(x[y], na.rm = TRUE); b <- mean(x[!y], na.rm = TRUE); list(cut = (a + b) / 2, up = a > b) }
  pred <- map_lgl(seq_along(v), \(i) { f <- fit(v[-i], is_rel[-i]); if (is.na(v[i])) NA else if (f$up) v[i] > f$cut else v[i] < f$cut })
  f <- fit(v, is_rel)
  tibble(cutoff_NPX = f$cut, loo_correct = sum(pred == is_rel, na.rm = TRUE), loo_n = sum(!is.na(pred)),
         loo_sensitivity = sum(pred & is_rel, na.rm = TRUE) / sum(is_rel), loo_specificity = sum(!pred & !is_rel, na.rm = TRUE) / sum(!is_rel))
}

#' Can V1 predict relapse at all? The whole procedure inside leave-one-patient-out: on the other patients pick the k
#' proteins with the largest |limma t|, score = mean of their signed z values, cut-off = midpoint of the group means,
#' classify the left-out patient. Returns the number of correctly predicted patients.
nested_loo_correct <- function(M, is_rel, k = 1) {
  sum(map_lgl(seq_along(is_rel), \(i) {
    tr <- M[, -i, drop = FALSE]; y <- is_rel[-i]
    r <- limma_t(tr, y); top <- head(order(-abs(r$t)), k); sg <- sign(r$t[top])
    mu <- rowMeans(tr[top, , drop = FALSE], na.rm = TRUE); sdv <- apply(tr[top, , drop = FALSE], 1, sd, na.rm = TRUE)
    sc <- \(x) colMeans(sg * (x[top, , drop = FALSE] - mu) / sdv, na.rm = TRUE)
    s_tr <- sc(tr); cut <- (mean(s_tr[y]) + mean(s_tr[!y])) / 2
    isTRUE((sc(M[, i, drop = FALSE]) > cut) == is_rel[i])          # no value for the left-out patient = not predicted
  }))
}

#' Exact AUC of every protein for every labelling (ranks do not depend on the labels): number of proteins that
#' separate the groups perfectly (AUC 1 or 0) per labelling.
perfect_separators <- function(M, combs) {
  R <- t(apply(M, 1, rank, na.last = "keep")); n1 <- nrow(combs); n0 <- ncol(M) - n1
  apply(combs, 2, \(cm) { a <- (rowSums(R[, cm, drop = FALSE]) - n1 * (n1 + 1) / 2) / (n1 * n0); sum(a %in% c(0, 1), na.rm = TRUE) })
}
