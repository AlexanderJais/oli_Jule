# Single-case analysis of one Leipzig biobank sample against the other Leipzig samples (step 08b):
# which proteins differ, and what the proteome says about the person's clinical values.
# Only Leipzig (LEIP) serum samples are used.

clinical_cols <- c("age", "sex", "BMI", "WHR", "c_fett", "HOMA_IR", "c_CRP", "C_CHOL", "C_HDL", "C_LDL",
                   "C_TRIGLY", "c_apo", "MDRD_kurz", "Gluc0_mg_dl")

#' The Leipzig sample to profile: "auto" = the one Leipzig sample without clinical data; otherwise a SampleID or SubjectID.
#' Returns NULL (step skipped) if "auto" finds no such sample or several.
resolve_case <- function(samples, case = "auto") {
  num <- intersect(setdiff(clinical_cols, "sex"), names(samples))
  if (identical(case, "auto") || is.null(case)) {
    if (!length(num)) { msg("No LEIP clinical data (paths$leip_clinical) - step 08b skipped; set leip_case$case to profile a sample."); return(NULL) }
    no_clin <- samples$SampleID[rowSums(!is.na(samples[num])) == 0]
    if (length(no_clin) != 1) {
      msg("leip_case$case = auto: %d Leipzig samples without clinical data (%s) - step 08b skipped; set leip_case$case in config.yml.",
          length(no_clin), paste(no_clin, collapse = ", "))
      return(NULL)
    }
    return(no_clin)
  }
  hit <- samples$SampleID[samples$SampleID == case | samples$SubjectID == case]
  if (length(hit) != 1) stop(sprintf("leip_case$case = '%s' matches %d Leipzig serum samples (need exactly 1).", case, length(hit)), call. = FALSE)
  hit
}

#' Clinical values as numbers: "<0.3" -> 0.15 (half the limit), ">60" -> 60, decimal commas allowed.
parse_clinical <- function(v) {
  if (is.numeric(v)) return(v)
  v <- str_replace(str_trim(as.character(v)), ",", ".")
  num <- suppressWarnings(as.numeric(str_remove(v, "^[<>]=?\\s*")))
  if_else(str_detect(coalesce(v, ""), "^<"), num / 2, num)
}

#' Sex as "M" / "F" (also accepts m/f, w, male/female, maennlich/weiblich); anything else becomes NA.
parse_sex <- function(v) {
  v <- str_to_lower(str_trim(as.character(v)))
  case_when(v %in% c("m", "male", "man", "mann", "maennlich", "männlich") ~ "M",
            v %in% c("f", "w", "female", "woman", "frau", "weiblich") ~ "F", TRUE ~ NA_character_)
}

#' Samples x proteins matrix (one column per Assay name; duplicate assays averaged).
case_matrix <- function(d, col = "value") {
  w <- d |> group_by(SampleID, Assay) |> summarise(v = mean(.data[[col]], na.rm = TRUE), .groups = "drop") |>
    mutate(v = if_else(is.nan(v), NA_real_, v)) |>
    pivot_wider(names_from = Assay, values_from = v)
  m <- as.matrix(w[, -1]); rownames(m) <- w$SampleID; m
}

#' Crawford & Howell (1998) single-case t-test of one value x against reference matrix R (samples x proteins).
crawford_howell <- function(x, R) {
  n <- colSums(!is.na(R)); m <- colMeans(R, na.rm = TRUE); s <- apply(R, 2, sd, na.rm = TRUE)
  med <- apply(R, 2, median, na.rm = TRUE); mad <- apply(R, 2, mad, na.rm = TRUE)
  t <- (x - m) / (s * sqrt((n + 1) / n))
  tibble(Assay = colnames(R), n_ref = n, case_value = x, ref_mean = m, ref_sd = s, ref_median = med,
         ref_min = apply(R, 2, min, na.rm = TRUE), ref_max = apply(R, 2, max, na.rm = TRUE),
         diff = x - med, z = (x - m) / s, robust_z = (x - med) / mad, t = t,
         p = 2 * pt(-abs(t), df = n - 1),
         percentile = 100 * (colSums(R < x, na.rm = TRUE) + 0.5 * colSums(R == x, na.rm = TRUE)) / n)
}

#' Test one case against the reference: proteins detected in >= min_detect of the reference are tested.
case_test <- function(x, R, below_R, min_detect = 0.5, fdr = 0.05, min_diff = 0.5) {
  tested <- colMeans(!below_R, na.rm = TRUE) >= min_detect & colSums(!is.na(R)) >= 5 & !is.na(x)
  crawford_howell(x[tested], R[, tested, drop = FALSE]) |>
    mutate(q = p.adjust(p, "BH"), outside_range = case_value > ref_max | case_value < ref_min,
           hit = q < fdr & abs(diff) >= min_diff & outside_range)
}

#' Correlation (rank based) of every protein with each clinical value in the reference samples.
clinical_assoc <- function(R, clin) {
  vars <- intersect(setdiff(clinical_cols, "sex"), names(clin))
  out <- map(vars, \(v) {
    y <- clin[[v]]; ok <- !is.na(y)
    if (sum(ok) < 10) return(NULL)
    rho <- suppressWarnings(cor(rank(y[ok]), apply(R[ok, , drop = FALSE], 2, rank, na.last = "keep"), use = "pairwise.complete.obs"))[1, ]
    n <- sum(ok); tt <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-12))
    tibble(Assay = colnames(R), parameter = v, rho = rho, p = 2 * pt(-abs(tt), n - 2))
  })
  if ("sex" %in% names(clin) && length(unique(na.omit(clin$sex))) == 2) {
    ok <- !is.na(clin$sex)
    out <- c(out, list(tibble(Assay = colnames(R), parameter = "sex (M - F)",
      rho = apply(R[ok, , drop = FALSE], 2, \(v) median(v[clin$sex[ok] == "M"], na.rm = TRUE) - median(v[clin$sex[ok] == "F"], na.rm = TRUE)),
      p = apply(R[ok, , drop = FALSE], 2, \(v) tryCatch(suppressWarnings(wilcox.test(v ~ clin$sex[ok], exact = FALSE)$p.value), error = \(e) NA_real_)))))
  }
  bind_rows(out) |> group_by(parameter) |> mutate(fdr = p.adjust(p, "BH")) |> ungroup()
}

# ---- clinical profile -----------------------------------------------------------------------------------------------

#' z-scores against the reference samples (mean / SD of the reference only).
ref_z <- function(X, ref_ids) {
  m <- colMeans(X[ref_ids, , drop = FALSE], na.rm = TRUE); s <- apply(X[ref_ids, , drop = FALSE], 2, sd, na.rm = TRUE)
  s[!is.na(s) & s == 0] <- NA
  sweep(sweep(X, 2, m), 2, s, "/")
}

#' Weighted marker score: sum(w * z) / sum(|w|) over the available markers (NA-skipping); clock = TRUE gives sum(w * z).
marker_score <- function(Z, w, clock = FALSE) {
  a <- intersect(names(w), colnames(Z))
  if (!length(a)) return(rep(NA_real_, nrow(Z)))
  Zs <- Z[, a, drop = FALSE]; W <- matrix(w[a], nrow(Zs), length(a), byrow = TRUE) * !is.na(Zs)
  Zs[is.na(Zs)] <- 0
  s <- rowSums(Zs * W)
  if (clock) s else s / rowSums(abs(W))
}

#' Sex from pre-specified sex-specific proteins, calibrated on the reference samples (leave-one-out).
predict_sex <- function(Z, below, ref_ids, case_id, sex, w, min_detect = 0.5, cap = c(0.03, 0.97)) {
  usable <- function(ids) {          # marker detected in >= min_detect of the sex it is expected in
    a <- intersect(names(w), colnames(Z))
    a[map_lgl(a, \(m) { pos <- ids[sex[ids] == if (w[m] > 0) "M" else "F"]; length(pos) > 0 && mean(!below[pos, m], na.rm = TRUE) >= min_detect })]
  }
  fit <- function(ids) {
    a <- usable(ids)
    if (length(a) == 0 || length(unique(sex[ids])) < 2) return(NULL)
    sc <- marker_score(Z[ids, a, drop = FALSE], w[a])
    mM <- mean(sc[sex[ids] == "M"]); mF <- mean(sc[sex[ids] == "F"])
    sd_p <- sqrt(sum(c(sc[sex[ids] == "M"] - mM, sc[sex[ids] == "F"] - mF)^2) / (length(ids) - 2))
    list(markers = a, mM = mM, mF = mF, sd = max(sd_p, 1e-6))
  }
  ref_ids <- ref_ids[!is.na(sex[ref_ids])]
  loo <- map_chr(ref_ids, \(i) {
    f <- fit(setdiff(ref_ids, i)); if (is.null(f)) return(NA_character_)
    s <- marker_score(Z[i, f$markers, drop = FALSE], w[f$markers]); if (s > (f$mM + f$mF) / 2) "M" else "F"
  })
  f <- fit(ref_ids)
  if (is.null(f)) return(list(call = "indeterminate", p_male = 0.5, reason = "no usable sex markers", markers = tibble(), loo_errors = NA, n_ref = length(ref_ids), ref_scores = tibble()))
  sc <- marker_score(Z[case_id, f$markers, drop = FALSE], w[f$markers])
  lr <- dnorm(sc, f$mM, f$sd) / dnorm(sc, f$mF, f$sd)
  p_male <- min(max(lr / (1 + lr), cap[1]), cap[2])
  mk <- tibble(Assay = f$markers, expected_higher_in = if_else(w[f$markers] > 0, "male", "female"),
               case_z = Z[case_id, f$markers],
               mean_z_male = colMeans(Z[ref_ids[sex[ref_ids] == "M"], f$markers, drop = FALSE], na.rm = TRUE),
               mean_z_female = colMeans(Z[ref_ids[sex[ref_ids] == "F"], f$markers, drop = FALSE], na.rm = TRUE),
               case_below_LOD = below[case_id, f$markers]) |>
    mutate(vote = if_else(abs(case_z - mean_z_male) < abs(case_z - mean_z_female), "male", "female"))
  errors <- sum(loo != sex[ref_ids], na.rm = TRUE); net <- sum(mk$vote == "male") - sum(mk$vote == "female")
  call <- if (errors <= 1 && p_male >= 0.9 && net >= min(4, nrow(mk))) "male"
          else if (errors <= 1 && p_male <= 0.1 && -net >= min(4, nrow(mk))) "female" else "indeterminate"
  reason <- if (call != "indeterminate") "" else if (errors > 1) sprintf("%d of %d reference samples misclassified", errors, length(ref_ids))
            else "markers do not agree clearly"
  list(call = call, p_male = p_male, reason = reason, markers = mk, loo_errors = errors, n_ref = length(ref_ids),
       ref_scores = tibble(SampleID = c(ref_ids, case_id), sex = c(sex[ref_ids], "case"),
                           score = marker_score(Z[c(ref_ids, case_id), f$markers, drop = FALSE], w[f$markers]),
                           loo_call = c(loo, NA)))
}

#' Fit y ~ score (+ sex) on the reference samples with leave-one-out evaluation, a permutation test and a
#' jackknife+ 80% interval for the case. log = TRUE fits log(y). Returns LOO predictions and the case prediction.
fit_profile <- function(y, score, sex = NULL, case_score, log = FALSE, n_perm = 999, level = 0.8) {
  ok <- !is.na(y) & !is.na(score) & (if (is.null(sex)) TRUE else !is.na(sex))
  if (log) ok <- ok & y > 0
  yt <- if (log) base::log(y[ok]) else y[ok]
  X <- cbind(1, score[ok], if (!is.null(sex)) as.numeric(sex[ok] == "M"))
  n <- length(yt)
  if (n < 10 || qr(X)$rank < ncol(X)) return(NULL)
  loo_pred <- function(yy) {                       # exact leave-one-out predictions of a linear model
    f <- lm.fit(X, yy); h <- rowSums((X %*% solve(crossprod(X))) * X)
    yy - f$residuals / (1 - h)
  }
  r2cv <- function(yy) { p <- loo_pred(yy); b <- (sum(yy) - yy) / (n - 1); 1 - sum((yy - p)^2) / sum((yy - b)^2) }
  pl <- loo_pred(yt); r2 <- r2cv(yt)
  null <- replicate(n_perm, r2cv(sample(yt)))
  back <- if (log) exp else identity
  base_loo <- (sum(yt) - yt) / (n - 1)
  # jackknife+ (Barber et al. 2021): models without sample i, evaluated at the case
  pred_case <- function(cs) {
    xc <- c(1, case_score, if (!is.null(sex)) as.numeric(cs == "M"))
    mu <- map_dbl(seq_len(n), \(i) sum(lm.fit(X[-i, , drop = FALSE], yt[-i])$coefficients * xc))
    r <- abs(yt - pl)
    k_lo <- floor((1 - level) * (n + 1)); k_hi <- ceiling(level * (n + 1))      # |residuals| make each side use 1 - level
    full <- sum(lm.fit(X, yt)$coefficients * xc)
    list(est = full, lo = sort(mu - r)[max(k_lo, 1)], hi = sort(mu + r)[min(k_hi, n)], dist = mu + (yt - pl))
  }
  cases <- if (is.null(sex)) list(pred_case(NULL)) else set_names(map(c("M", "F"), pred_case), c("M", "F"))
  list(n = n, r2cv = r2, p_perm = (1 + sum(null >= r2)) / (n_perm + 1),
       mae = mean(abs(back(yt) - back(pl))), mae_baseline = mean(abs(back(yt) - back(base_loo))),
       loo = tibble(SampleID = names(y)[ok], observed = back(yt), predicted = back(pl)), cases = cases, back = back,
       score_pct = 100 * mean(score[ok] < case_score), score_outside = case_score < min(score[ok]) | case_score > max(score[ok]))
}

#' Combine the case prediction over the sex call: called sex, or both sexes weighted by P(male) if indeterminate.
case_prediction <- function(fp, sex_call, p_male) {
  b <- fp$back
  if (length(fp$cases) == 1 || sex_call != "indeterminate") {
    cp <- if (length(fp$cases) == 1) fp$cases[[1]] else fp$cases[[if (sex_call == "male") "M" else "F"]]
    return(list(est = b(cp$est), lo = b(cp$lo), hi = b(cp$hi), dist = b(cp$dist), w = rep(1, length(cp$dist))))
  }
  m <- fp$cases$M; f <- fp$cases$F
  list(est = b(p_male * m$est + (1 - p_male) * f$est), lo = b(min(m$lo, f$lo)), hi = b(max(m$hi, f$hi)),
       dist = b(c(m$dist, f$dist)), w = c(rep(p_male, length(m$dist)), rep(1 - p_male, length(f$dist))))
}

#' Probability of each band (e.g. BMI classes) from the jackknife+ predictive distribution.
band_probs <- function(dist, w, breaks, labels) {
  b <- cut(dist, c(-Inf, breaks, Inf), labels = labels, right = FALSE)
  p <- tapply(w, b, sum); p[is.na(p)] <- 0; p <- p / sum(p)
  paste(sprintf("%s %.0f%%", names(p), 100 * p), collapse = "; ")
}

#' Which apolipoprotein is c_apo? From its correlation with HDL and LDL and its skew in the reference (no proteins used).
apo_analyte <- function(clin, setting = "auto") {
  if (!identical(setting, "auto")) return(setting)
  if (!"c_apo" %in% names(clin) || sum(!is.na(clin$c_apo)) < 10) return(NA_character_)
  r_hdl <- suppressWarnings(cor(clin$c_apo, clin$C_HDL, use = "complete.obs", method = "spearman"))
  r_ldl <- suppressWarnings(cor(clin$c_apo, clin$C_LDL, use = "complete.obs", method = "spearman"))
  x <- na.omit(clin$c_apo); skew <- mean((x - mean(x))^3) / sd(x)^3
  if (isTRUE(r_hdl >= 0.5 && r_hdl > r_ldl)) "ApoA1" else if (isTRUE(r_ldl >= 0.5)) "ApoB" else if (skew > 1) "Lpa" else NA_character_
}
