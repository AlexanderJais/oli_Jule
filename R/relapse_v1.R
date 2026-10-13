# Helpers for the V1 relapse screen (run_relapse_v1.R): which dISF proteins are already higher at baseline in patients
# who later relapse? Built for very small groups (e.g. 4 vs 6): exact relabelling check, leave-one-out, AUC.

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
  for (j in seq_len(ncol(combs))) {
    lab <- seq_len(n) %in% combs[, j]
    r <- limma_t(M, lab)
    ts[, j] <- r$t; hits[j] <- sum(r$p < 0.05 & r$logFC > 0, na.rm = TRUE)
  }
  list(per_protein = obs |> mutate(p_exact_up = rowMeans(ts >= t - 1e-12, na.rm = TRUE), p_exact_down = rowMeans(ts <= t + 1e-12, na.rm = TRUE),
                                   p_exact_two_sided = rowMeans(abs(ts) >= abs(t) - 1e-12, na.rm = TRUE)),
       hits = hits, observed_hits = sum(obs$p < 0.05 & obs$logFC > 0, na.rm = TRUE), n_labellings = ncol(combs))
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
