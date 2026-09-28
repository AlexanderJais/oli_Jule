# 02 - Aim 1: does Olink confirm the galanin ELISA?
# Galanin is measured twice in the same LEIP sera: by ELISA (pg/mL) and by Olink (assay GAL; NPX,
# a relative log2 value). The two can only agree in ranking the samples, not in absolute values.
#   - is Olink GAL above its LOD in LEIP serum?
#   - agreement: Spearman (all samples; adjusted for Olink plate; adjusted for sex and plate; women and
#     men separately; within each plate; above LOD only; without flagged samples; without extreme ELISA
#     values), leave-one-out range, tertile agreement (weighted kappa), and the change in NPX per
#     doubling of the ELISA value (1 = same fold-change)
#   - specificity: the rank of Olink GAL among all proteins correlated with the ELISA, and whether the
#     ELISA and Olink GAL correlate with the same proteins (protein profiles, permutation test)
#   - benchmark: the same comparison for the other proteins measured by the lab and by Olink; it also
#     shows whether the clinical file and the Olink samples are correctly matched - overall, and per
#     person (does the profile of lab values match the person's own Olink sample best?)
#   - technical factors: Olink plate, ELISA plate
# Out: output/1_elisa_validation/

source("R/utils.R")
source("R/leip.R")
cfg <- load_config()
set.seed(cfg$seed %||% 1)
d <- leip_load(cfg)
clear_outputs(cfg, "1_elisa_validation")
Q <- "1 Does Olink confirm the galanin ELISA?"
ans <- answers_new()
if (!d$has_gal || !d$has_elisa) {
  answer(ans, Q, "Olink GAL vs ELISA", "not possible",
         if (!d$has_gal) "Galanin (GAL) is not in the Olink data." else "The clinical file has no galanin ELISA values.")
  answers_save(ans, cfg, "1_elisa_validation", "answers.csv")
  quit(save = "no")
}
S <- d$S
gdet <- d$det |> filter(OlinkID == d$gal)
gd <- S |>
  select(SampleID, SubjectID, plate, any_of(c("SampleQC", "qc_outlier", "sex_male", "age", "C_HDL", "c_apo",
                                              "low_serum", "lipamisch", "Galanin_ELISA_plate"))) |>
  mutate(elisa = d$elisa, npx = d$npx) |>
  left_join(d$lod |> filter(OlinkID == d$gal) |> select(SampleID, LOD, below_lod), by = "SampleID")

# ---- is GAL measurable? -----------------------------------------------------------------------------------------------------
fa <- gdet$frac_above_lod
answer(ans, Q, "Is Olink GAL above its LOD in LEIP serum?",
       case_when(is.na(fa) ~ "unknown (no LOD)", fa >= 0.9 ~ "yes, in (nearly) all samples", fa >= 0.5 ~ "yes, in most samples",
                 fa >= 0.2 ~ "partly: many samples are below LOD", TRUE ~ "mostly below LOD: the Olink values are mostly noise"),
       sprintf("%d of %d samples above LOD (%.0f%%); median NPX %.2f, median LOD %.2f. The assay targets UniProt %s (the galanin precursor).",
               sum(!gd$below_lod, na.rm = TRUE), sum(!is.na(gd$below_lod)), 100 * fa, gdet$median_npx, gdet$median_lod,
               col_or(gdet, "UniProt", "P22466")))

# ---- agreement ---------------------------------------------------------------------------------------------------------------------
pair <- gd |> filter(!is.na(npx), coalesce(elisa > 0, FALSE)) |>
  mutate(log2_elisa = log2(elisa), discordance = zscore(log2_elisa) - zscore(npx))
agree_row <- function(what, keep = rep(TRUE, nrow(pair)), covs = NULL, min_pairs = 10) {
  x <- pair[keep, ]
  if (!is.null(covs)) x <- x[stats::complete.cases(x[covs]), ]
  if (nrow(x) < min_pairs) return(tibble(analysis = what, n = nrow(x)))
  r <- spearman1(x$elisa, x$npx, if (!is.null(covs)) x[covs], min_n = min_pairs)
  f <- tryCatch(lm(reformulate(c("log2_elisa", covs), "npx"), data = x), error = \(e) NULL)
  ci <- if (!is.null(f)) suppressWarnings(confint(f)["log2_elisa", ]) else c(NA, NA)
  tibble(analysis = what, n = r$n, rho = r$rho, ci_low = r$ci_low, ci_high = r$ci_high, p = r$p,
         pearson_r_log2 = if (is.null(covs)) cor(x$log2_elisa, x$npx) else NA_real_,
         npx_per_doubling = if (!is.null(f)) unname(coef(f)["log2_elisa"]) else NA_real_,
         npx_per_doubling_ci_low = unname(ci[1]), npx_per_doubling_ci_high = unname(ci[2]))
}
flagged <- coalesce(pair$SampleQC != "PASS", FALSE) | coalesce(pair$qc_outlier, FALSE) |
  coalesce(col_or(pair, "low_serum", 0) == 1, FALSE) | coalesce(col_or(pair, "lipamisch", 0) == 1, FALSE)
has_sex <- "sex_male" %in% names(pair) && n_distinct(na.omit(pair$sex_male)) == 2
agree <- bind_rows(
  agree_row("all samples"),
  agree_row("adjusted for Olink plate", covs = "plate"),
  if (has_sex) list(agree_row("adjusted for sex and Olink plate", covs = c("sex_male", "plate")),
                    agree_row("women only", coalesce(pair$sex_male == 0, FALSE), min_pairs = 5),
                    agree_row("men only", coalesce(pair$sex_male == 1, FALSE), min_pairs = 5)),
  agree_row("only samples with GAL above LOD", coalesce(!pair$below_lod, FALSE)),
  agree_row("without flagged samples (QC warning / outlier, low serum, lipaemic)", !flagged),
  agree_row("without extreme ELISA values (|z| > 3, log scale)", abs(zscore(pair$log2_elisa)) <= 3),
  map(sort(unique(na.omit(pair$plate))), \(pl) agree_row(paste("within", pl), coalesce(pair$plate == pl, FALSE), min_pairs = 5))) |>
  ensure(c("rho", "ci_low", "ci_high", "p", "pearson_r_log2", "npx_per_doubling"))
a0 <- agree |> filter(analysis == "all samples"); a1 <- agree |> filter(analysis == "adjusted for Olink plate")
loo <- tibble(left_out = pair$SubjectID,
              rho = map_dbl(seq_len(nrow(pair)), \(i) cor(pair$elisa[-i], pair$npx[-i], method = "spearman"))) |>
  mutate(change = rho - a0$rho) |> arrange(desc(abs(change)))
tt <- pair |> mutate(ELISA = factor(ntile(elisa, 3), 1:3, c("low", "middle", "high")), Olink = factor(ntile(npx, 3), 1:3, c("low", "middle", "high")))
tab <- table(ELISA = tt$ELISA, Olink = tt$Olink)
tert_tab <- as.data.frame.matrix(tab) |> rownames_to_column("ELISA tertile (rows) / Olink tertile (columns)")
tert <- tibble(n = sum(tab), same_tertile_pct = 100 * sum(diag(tab)) / sum(tab),
               opposite_tertile_pct = 100 * (tab[1, 3] + tab[3, 1]) / sum(tab), weighted_kappa = weighted_kappa(tab))

# ---- technical factors -------------------------------------------------------------------------------------------------------------
kw <- function(v, g, what) {
  ok <- !is.na(v) & !is.na(g)
  if (n_distinct(g[ok]) < 2) return(NULL)
  m <- tapply(v[ok], as.character(g[ok]), median)
  m <- m[str_order(names(m), numeric = TRUE)]                  # plate 2 before plate 10
  tibble(test = what, n = sum(ok), groups = length(m), medians = paste(sprintf("%s: %.3g", names(m), m), collapse = "; "),
         p = kruskal.test(v[ok], factor(g[ok]))$p.value)
}
plate_fx <- bind_rows(kw(gd$elisa, gd$plate, "galanin ELISA by Olink plate (were samples put on plates by level?)"),
                      kw(gd$npx, gd$plate, "Olink GAL by Olink plate"),
                      if ("Galanin_ELISA_plate" %in% names(gd)) kw(gd$elisa, gd$Galanin_ELISA_plate, "galanin ELISA by ELISA plate"))

# ---- specificity: which proteins does the ELISA follow? ---------------------------------------------------------------------------
Ye <- if (d$gal %in% colnames(d$M)) d$M else cbind(d$M, d$Y[, d$gal, drop = FALSE])
u <- spearman_vs(d$elisa, Ye)
a <- spearman_vs(d$elisa, Ye, S[d$cov])
ev <- tibble(OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p, rho_adj = a$rho, p_adj = a$p) |>
  mutate(fdr = p.adjust(p, "BH"), fdr_adj = p.adjust(p_adj, "BH"), rank_by_rho = rank(-rho, ties.method = "min", na.last = "keep"),
         is_galanin = OlinkID == d$gal) |>
  left_join(d$det |> select(OlinkID, Assay, protein, frac_above_lod), by = "OlinkID") |>
  relocate(OlinkID, Assay, protein) |> arrange(p)
gal_rank <- ev$rank_by_rho[ev$is_galanin]; n_rank <- sum(!is.na(ev$rho))

# do Olink GAL and the ELISA correlate with the same proteins? (same samples; complete proteins; permutation of the ELISA)
ok <- !is.na(d$elisa) & !is.na(d$npx)
Mx <- d$M[ok, setdiff(colnames(d$M), d$gal), drop = FALSE]
Mx <- Mx[, colSums(is.na(Mx)) == 0 & apply(Mx, 2, sd) > 0, drop = FALSE]
R <- apply(Mx, 2, rank)
prof <- tibble(OlinkID = colnames(Mx), rho_GAL = as.numeric(cor(rank(d$npx[ok]), R)), rho_ELISA = as.numeric(cor(rank(d$elisa[ok]), R))) |>
  left_join(d$det |> select(OlinkID, protein), by = "OlinkID")
prof_r <- cor(prof$rho_GAL, prof$rho_ELISA)
B <- cfg$permutations %||% 1000
E <- replicate(B, sample(rank(d$elisa[ok])))
null_r <- as.numeric(cor(t(cor(E, R)), prof$rho_GAL))
prof_p <- (1 + sum(null_r >= prof_r)) / (B + 1)

# ---- benchmark: lab vs Olink for other proteins ----------------------------------------------------------------------------------
bench <- imap(cfg$lab_vs_olink %||% list(galanin_elisa = "GAL"), \(assay, lab) {
  row <- tibble(lab = lab, label = param_label(lab), olink_assay = assay)
  if (!lab %in% names(S)) return(row |> mutate(status = "not in the clinical file"))
  if (sum(!is.na(S[[lab]])) < 10) return(row |> mutate(status = "fewer than 10 values"))
  oid <- find_assay(d$det, assay)[1]
  if (is.na(oid)) return(row |> mutate(status = "not measured by Olink"))
  allr <- spearman_vs(S[[lab]], if (oid %in% colnames(d$M)) d$M else cbind(d$M, d$Y[, oid, drop = FALSE]))
  r <- allr |> filter(id == oid)
  row |> mutate(status = "compared", OlinkID = oid, n = r$n, frac_above_lod = d$det$frac_above_lod[d$det$OlinkID == oid],
                rho = r$rho, ci_low = r$ci_low, ci_high = r$ci_high, p = r$p,
                rank_among_proteins = rank(-allr$rho, ties.method = "min", na.last = "keep")[allr$id == oid],
                n_proteins = sum(!is.na(allr$rho)))
}) |> bind_rows() |>
  ensure(c("n", "frac_above_lod", "rho", "ci_low", "ci_high", "p", "rank_among_proteins", "n_proteins")) |>
  arrange(desc(status == "compared"), desc(rho))
others <- bench |> filter(status == "compared", lab != "galanin_elisa", !is.na(rho))

# ---- per person: is each Olink sample the right one? ------------------------------------------------------------------------------
# The lab values of the benchmark proteins that agree with Olink (rho >= 0.5) form a profile per person; so do their
# Olink values. Distance between the lab profile of person i and the Olink profile of sample j: rank-based normal
# scores, mean squared difference weighted by how well each protein agrees (1 / (2 (1 - rho)), the expected squared
# difference of a correct pair). A person's own sample should be among the best matches; a swapped or mislabelled
# sample fits no better than a stranger's - also where a single swap barely changes the correlations over all persons.
nscore <- \(v) { r <- rep(NA_real_, length(v)); k <- !is.na(v); r[k] <- qnorm((rank(v[k]) - 0.5) / sum(k)); r }
idp <- others |> filter(rho >= 0.5)
ident <- tibble(SampleID = character(), SubjectID = character(), own_rank = integer(), of = integer(), own_distance = numeric(),
                best_match = character(), fit = character())
if (nrow(idp) >= 3) {
  L <- sapply(idp$lab, \(l) nscore(S[[l]]))                     # persons x proteins: lab
  O <- sapply(idp$OlinkID, \(o) nscore(d$Y[, o]))              # samples x proteins: Olink
  w <- 1 / (2 * (1 - pmin(idp$rho, 0.95)))
  D <- sapply(seq_len(nrow(O)), \(j) {                         # D[i, j]: lab of person i vs Olink of sample j
    sq <- sweep(L, 2, O[j, ])^2; ok <- !is.na(sq)
    sqrt(as.vector(ifelse(ok, sq, 0) %*% w) / as.vector(ok %*% w))
  })
  D[rowSums(!is.na(L)) < 3, ] <- NA
  best <- apply(D, 1, \(r) if (all(is.na(r))) NA_integer_ else which.min(r))
  ident <- tibble(SampleID = S$SampleID, SubjectID = S$SubjectID, plate = S$plate, proteins = rowSums(!is.na(L)),
                  own_distance = diag(D),
                  own_rank = map_int(seq_len(nrow(D)), \(i) if (is.na(D[i, i])) NA_integer_ else as.integer(sum(D[i, ] < D[i, i], na.rm = TRUE) + 1)),
                  of = colSums(!is.na(t(D))), best_match = S$SubjectID[best], best_distance = D[cbind(seq_len(nrow(D)), best)]) |>
    left_join(pair |> select(SampleID, galanin_discordance = discordance), by = "SampleID") |>
    filter(!is.na(own_rank)) |>
    mutate(fit = case_when(own_rank == 1 ~ "best match", own_rank <= of / 4 ~ "among the best 25%",
                           TRUE ~ "possible swap (no better than others)")) |>
    arrange(desc(own_rank), desc(own_distance))
}

# ---- answers -------------------------------------------------------------------------------------------------------------------------
verdict <- case_when(is.na(a0$p) ~ "not enough data",
                     a0$p < 0.05 & a0$rho >= 0.5 ~ "yes: both methods rank the samples similarly",
                     a0$p < 0.05 & a0$rho >= 0.3 ~ "partly: moderate agreement",
                     a0$p < 0.05 & a0$rho < 0 ~ "no: opposite ranking (unexpected)",
                     TRUE ~ "no: the Olink values do not follow the ELISA")
pl <- if (nrow(plate_fx)) plate_fx |> filter(str_detect(test, "^galanin ELISA by Olink plate")) else tibble()
ro <- \(a) agree |> filter(analysis == a)
sx_txt <- if (has_sex) with(list(s1 = ro("adjusted for sex and Olink plate"), w = ro("women only"), m = ro("men only")),
  sprintf("Adjusted for sex and Olink plate rho = %.2f (p = %s); women rho = %.2f (p = %s, n = %d), men rho = %.2f (p = %s, n = %d). ",
          s1$rho, fmt_p(s1$p), w$rho, fmt_p(w$p), w$n, m$rho, fmt_p(m$p), m$n)) else ""
answer(ans, Q, "Olink GAL vs ELISA", verdict,
       paste0(sprintf("Spearman rho = %.2f (95%% CI %.2f to %.2f), p = %s, n = %d; adjusted for Olink plate rho = %.2f (p = %s). ",
                      a0$rho, a0$ci_low, a0$ci_high, fmt_p(a0$p), a0$n, a1$rho, fmt_p(a1$p)),
              sx_txt,
              sprintf("Leaving out one person gives rho %.2f to %.2f. ", min(loo$rho), max(loo$rho)),
              sprintf("Same tertile in %.0f%% of persons (chance 33%%), opposite tertiles in %.0f%%; weighted kappa %.2f. ",
                      tert$same_tertile_pct, tert$opposite_tertile_pct, tert$weighted_kappa),
              sprintf("Olink NPX change per doubling of the ELISA value: %.2f (95%% CI %.2f to %.2f; 1 = same fold-change). ",
                      a0$npx_per_doubling, a0$npx_per_doubling_ci_low, a0$npx_per_doubling_ci_high),
              if (nrow(pl) && isTRUE(pl$p[1] < 0.05)) sprintf("The ELISA values differ between the Olink plates (p = %s): see the plate-adjusted value. ", fmt_p(pl$p[1])) else ""))
answer(ans, Q, "Specificity: does the ELISA follow Olink GAL more than other proteins?",
       if (isTRUE(gal_rank <= 3)) "yes: Olink GAL is among the proteins the ELISA follows most"
       else if (isTRUE(gal_rank / n_rank <= 0.05)) "partly: Olink GAL is in the top 5% of the proteins the ELISA follows"
       else "no: the ELISA follows other proteins more than Olink GAL",
       paste0(sprintf("Olink GAL ranks %d of %d proteins by correlation with the ELISA. The ELISA correlates most with %s. ",
                      gal_rank, n_rank, top_str(ev$protein, ev$rho, ev$p, 6)),
              sprintf("Protein profiles: across %d proteins, the correlations of Olink GAL and of the ELISA agree with r = %.2f (permutation p = %s).",
                      nrow(prof), prof_r, fmt_p(prof_p))))
answer(ans, Q, "Benchmark: how well do lab assays and Olink agree for other proteins?",
       if (!nrow(others)) "no other protein measured by both" else sprintf("median rho %.2f over %d proteins (galanin: %.2f)", median(others$rho), nrow(others), a0$rho),
       if (!nrow(others)) "" else paste(sprintf("%s vs Olink %s: rho %.2f (rank %d of %d)", others$label, others$olink_assay, others$rho,
                                                others$rank_among_proteins, others$n_proteins), collapse = "; "))
good <- others |> filter(rho >= 0.7, rank_among_proteins <= 5)
answer(ans, Q, "Are the clinical data and the Olink data of the same persons correctly matched?",
       if (!nrow(others)) "cannot be checked: no other protein measured by the lab and by Olink"
       else if (nrow(good) >= 3) sprintf("yes: for %d lab assays the matching Olink protein is the best or near-best match among all proteins", nrow(good))
       else "not confirmed: fewer than 3 lab assays agree clearly with Olink",
       paste0(if (nrow(good)) paste0(paste(sprintf("%s rho %.2f (rank %d)", good$label, good$rho, good$rank_among_proteins), collapse = "; "), ". ") else "",
              "A general mix-up between the rows of the clinical file and the Olink samples would destroy these correlations; a single swap would not (see the next answer). ",
              "It does not check the galanin ELISA values themselves: a mix-up on the ELISA plates, or an ELISA run on another aliquot or blood draw, would affect galanin only."))
top_g <- pair |> slice_max(abs(discordance), n = 3, with_ties = FALSE) |> left_join(ident |> select(SampleID, own_rank), by = "SampleID")
mism <- ident |> filter(own_rank > of / 4)
answer(ans, Q, "Per person: does each Olink sample belong to the right person?",
       if (!nrow(ident)) "cannot be checked: fewer than 3 lab assays agree with Olink (rho >= 0.5)"
       else if (!nrow(mism)) sprintf("yes: for all %d persons their own Olink sample is among the best matches of their lab values (best match: %d)",
                                     nrow(ident), sum(ident$own_rank == 1))
       else sprintf("check %d of %d persons: their own Olink sample fits their lab values no better than other persons' samples", nrow(mism), nrow(ident)),
       if (!nrow(ident)) "" else paste0(
         sprintf("Profile of %d proteins measured by the lab and by Olink (%s); own sample the best match for %d of %d persons. ", nrow(idp),
                 paste(idp$olink_assay, collapse = ", "), sum(ident$own_rank == 1), nrow(ident)),
         if (nrow(mism)) paste0("Possible swap or mix-up: ", paste(sprintf("%s (own sample rank %d of %d; best match: Olink sample of %s)", mism$SubjectID,
                                                                        mism$own_rank, mism$of, mism$best_match), collapse = "; "), ". ") else "",
         "The 3 persons where ELISA and Olink galanin disagree most: ",
         paste(sprintf("%s (own sample rank %s)", top_g$SubjectID, coalesce(as.character(top_g$own_rank), "n/a")), collapse = ", "),
         ". If their own samples fit well, the galanin disagreement is not a sample swap on the Olink side; it then lies in the galanin measurements (or a mix-up of the ELISA values only)."))
interp <- case_when(
  is.na(fa) | is.na(a0$rho) ~ "not enough data to judge",
  fa < 0.5 ~ "Olink GAL is mostly below LOD, so Olink cannot validate the ELISA here",
  a0$p < 0.05 & a0$rho >= 0.5 ~ "the ELISA is supported by an independent method",
  nrow(others) > 0 & median(others$rho) >= 0.5 & a0$rho < 0.3 ~ "other lab assays agree well with Olink, galanin does not: the disagreement is specific to galanin - different forms measured (precursor vs mature peptide), degradation of the peptide, the ELISA itself (specificity, matrix effects, plate-to-plate variation) or binding to HDL (see 3)",
  nrow(others) > 0 & median(others$rho) < 0.3 ~ "lab assays and Olink agree poorly in general in these samples (pre-analytics, sample age?), so galanin cannot be judged alone",
  TRUE ~ "partial agreement; see the details")
answer(ans, Q, "Interpretation", interp,
       "Olink measures relative amounts (NPX) and the assay targets the galanin precursor (UniProt P22466); a galanin ELISA may detect the mature peptide. Agreement can only be judged as ranking.")
answers <- answers_save(ans, cfg, "1_elisa_validation", "answers.csv")

# ---- figures --------------------------------------------------------------------------------------------------------------------------
F <- figs_new()
disc <- pair |> slice_max(abs(discordance), n = 3, with_ties = FALSE)
p1 <- ggplot(pair, aes(elisa, npx)) +
  (if (!is.na(gdet$median_lod)) geom_hline(yintercept = gdet$median_lod, linetype = 2, colour = "grey50")) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", linewidth = 0.6) +
  geom_point(aes(colour = plate, shape = if_else(coalesce(below_lod, FALSE), "below LOD", "above LOD")), size = 2.4) +
  geom_text(data = disc, aes(label = SubjectID), size = 2.8, vjust = -0.9, colour = "grey25") +
  scale_x_continuous(trans = "log2", breaks = pretty(pair$elisa, n = 6)) +
  scale_shape_manual(values = c(`above LOD` = 16, `below LOD` = 1)) +
  labs(title = "Galanin in LEIP serum: Olink (GAL) vs ELISA",
       subtitle = sprintf("Spearman rho = %.2f (95%% CI %.2f to %.2f), p = %s, n = %d; adjusted for Olink plate rho = %.2f\ndashed: median Olink LOD; labelled: the 3 persons where the methods disagree most",
                          a0$rho, a0$ci_low, a0$ci_high, fmt_p(a0$p), a0$n, a1$rho),
       x = "galanin ELISA (pg/mL, log2 scale)", y = "Olink GAL (NPX, log2)", colour = "Olink plate", shape = "Olink")
fig(F, "scatter", p1, cfg, "1_elisa_validation", "GAL_Olink_vs_ELISA.png", width = 8.5, height = 6)

zp <- pair |> mutate(mean_z = (zscore(log2_elisa) + zscore(npx)) / 2, HDL = col_or(pair, "C_HDL", NA_real_))
sd_d <- sd(zp$discordance)
p2 <- ggplot(zp, aes(mean_z, discordance)) +
  geom_hline(yintercept = c(-1.96, 1.96) * sd_d, linetype = 2, colour = "grey55") + geom_hline(yintercept = 0, colour = "grey40") +
  geom_point(aes(colour = HDL), size = 2.4) +
  geom_text(data = zp |> slice_max(abs(discordance), n = 3, with_ties = FALSE), aes(label = SubjectID), size = 2.8, vjust = -0.9) +
  scale_colour_gradient(low = "grey80", high = "firebrick", na.value = "grey50") +
  labs(title = "Where do ELISA and Olink disagree? (agreement on the z-score scale)",
       subtitle = "y > 0: ELISA higher than Olink relative to the other persons; dashed: +/- 1.96 SD; colour: HDL cholesterol (see 3)",
       x = "mean of the two z-scores", y = "z(ELISA) - z(Olink GAL)", colour = "HDL (mmol/l)")
fig(F, "z_agreement", p2, cfg, "1_elisa_validation", "agreement_z_scores.png", width = 8.5, height = 6)

b <- bench |> filter(status == "compared", !is.na(rho))
if (nrow(b)) {
  b <- b |> mutate(name = sprintf("%s  vs  Olink %s", label, olink_assay), galanin = lab == "galanin_elisa",
                   info = sprintf("%s above LOD, rank %d of %d", if_else(is.na(frac_above_lod), "?", sprintf("%.0f%%", 100 * frac_above_lod)),
                                  rank_among_proteins, n_proteins))
  p3 <- ggplot(b, aes(rho, reorder(name, rho), colour = galanin)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = ci_low, xmax = ci_high), size = 0.4) +
    geom_text(aes(x = 1.05, label = info), hjust = 0, size = 3, colour = "grey30") +
    scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "firebrick"), guide = "none") +
    scale_x_continuous(breaks = seq(-1, 1, 0.5)) + coord_cartesian(xlim = c(-1, 2)) +
    labs(title = "Benchmark: lab assay vs Olink for the same protein (LEIP serum)",
         subtitle = "Spearman rho with 95% CI; rank = place of the matching Olink assay among all proteins correlated with the lab value (1 = best)",
         x = "Spearman rho", y = NULL)
  fig(F, "benchmark", p3, cfg, "1_elisa_validation", "lab_vs_Olink_benchmark.png", width = 11, height = 2 + 0.35 * nrow(b))
}

e2 <- ev |> filter(!is.na(p))
p4 <- ggplot(e2, aes(rho, -log10(p))) + geom_point(colour = "grey60", size = 0.8, alpha = 0.7) +
  geom_text(data = head(e2, 12), aes(label = protein), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
  geom_point(data = e2 |> filter(is_galanin), colour = "firebrick", size = 3) +
  geom_text(data = e2 |> filter(is_galanin), aes(label = protein), colour = "firebrick", vjust = 1.8, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(title = "Which Olink proteins does the galanin ELISA follow?",
       subtitle = sprintf("red: Olink GAL, rank %d of %d proteins by correlation with the ELISA", gal_rank, n_rank),
       x = "Spearman rho with the ELISA", y = "-log10 p")
fig(F, "elisa_volcano", p4, cfg, "1_elisa_validation", "ELISA_vs_all_proteins.png", width = 9, height = 7)

lab_pr <- prof |> mutate(s = abs(rho_GAL) + abs(rho_ELISA)) |> slice_max(s, n = 8, with_ties = FALSE)
p5 <- ggplot(prof, aes(rho_GAL, rho_ELISA)) + geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(colour = "grey45", size = 0.8, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "firebrick", linewidth = 0.6) +
  geom_text(data = lab_pr, aes(label = protein), size = 2.7, vjust = -0.6, check_overlap = TRUE) +
  labs(title = "Do Olink GAL and the ELISA correlate with the same proteins?",
       subtitle = sprintf("each point = one protein; r = %.2f over %d proteins, permutation p = %s (%d permutations of the ELISA)",
                          prof_r, nrow(prof), fmt_p(prof_p), B),
       x = "rho with Olink GAL", y = "rho with the galanin ELISA")
fig(F, "profiles", p5, cfg, "1_elisa_validation", "protein_profiles_GAL_vs_ELISA.png", width = 8, height = 7)
figs_save(F, cfg, "1_elisa_validation", "figures.rds")

# ---- tables ---------------------------------------------------------------------------------------------------------------------------
samples <- pair |> select(SampleID, SubjectID, plate, elisa, log2_elisa, npx, LOD, below_lod, discordance, any_of(c("sex_male", "C_HDL", "c_apo")))
save_csv(agree, cfg, "1_elisa_validation", "agreement.csv")
if (nrow(plate_fx)) save_csv(plate_fx, cfg, "1_elisa_validation", "plate_effects.csv")
save_csv(bench, cfg, "1_elisa_validation", "lab_vs_Olink_benchmark.csv")
if (nrow(ident)) save_csv(ident, cfg, "1_elisa_validation", "sample_identity.csv")
save_csv(ev, cfg, "1_elisa_validation", "ELISA_vs_all_proteins.csv")
save_csv(samples, cfg, "1_elisa_validation", "galanin_values_per_sample.csv")
writexl::write_xlsx(list(answers = answers, GAL_detection = gdet, agreement = agree, leave_one_out = loo, tertiles = tert_tab,
                         tertile_agreement = tert, plate_effects = plate_fx, lab_vs_Olink = bench, sample_identity = ident, ELISA_vs_proteins = ev,
                         protein_profiles = prof |> arrange(desc(abs(rho_GAL) + abs(rho_ELISA))),
                         profile_test = tibble(proteins = nrow(prof), r = prof_r, permutations = B, p = prof_p),
                         samples = samples) |> keep(\(x) is.data.frame(x) && ncol(x) > 0),
                    out_path(cfg, "1_elisa_validation", "elisa_validation.xlsx"))
for (i in seq_len(nrow(answers))) msg("%s: %s", answers$item[i], answers$verdict[i])
