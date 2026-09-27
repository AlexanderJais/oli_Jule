# 04 - Aim 3: galanin and HDL. Hypothesis: galanin binds to HDL particles in the circulation.
# Correlations cannot prove binding, but they show whether the data fit the idea:
#   a) galanin vs the lab lipids (HDL-C, ApoA-I, LDL-C, total cholesterol, triglycerides, ApoB, Lp(a)):
#      all persons; adjusted for sex (women have higher HDL); adjusted for sex, age and BMI; women and
#      men separately; effect size from a linear model (% galanin per SD of HDL)
#   b) specificity: is galanin closer to HDL than to LDL, triglycerides or ApoB? (bootstrap)
#   c) HDL proteins measured by Olink (APOA1, APOA2, APOM, LCAT, PON1 ...): do they track the lab HDL,
#      and does galanin follow them (each protein, and an "HDL protein score")?
#   d) where does Olink GAL rank among all proteins correlated with HDL-C?
#   e) do the ELISA and Olink disagree more when HDL is high? (expected if one assay does not see
#      HDL-bound galanin)
# Out: output/3_galanin_hdl/

source("R/utils.R")
source("leip_galanin/R/leip.R")
cfg <- leip_config()
set.seed(cfg$seed %||% 1)
d <- leip_load(cfg)
S <- d$S
h <- cfg$hdl %||% list()
B <- cfg$bootstrap %||% 2000
clear_outputs(cfg, "3_galanin_hdl")
Q <- "3 Galanin and HDL"
ans <- answers_new()
hdl_col <- h$hdl %||% "C_HDL"; apo_col <- h$apoa1 %||% "c_apo"
lipids <- intersect(unlist(h$lipids %||% c("C_HDL", "c_apo", "C_LDL", "C_CHOL", "C_TRIGLY", "C_APO_B", "C_LIPO")), names(S))
measures <- compact(list(`galanin ELISA` = if (d$has_elisa) log2(d$elisa), `Olink GAL` = if (d$has_gal) d$npx))
if (!length(measures) || !hdl_col %in% names(S)) {
  answer(ans, Q, "galanin vs HDL", "not possible", if (!length(measures)) "No galanin measurement." else sprintf("No HDL column '%s'.", hdl_col))
  answers_save(ans, cfg, "3_galanin_hdl", "answers.csv")
  quit(save = "no")
}
has_sex <- "sex_male" %in% names(S) && n_distinct(na.omit(S$sex_male)) == 2
plate_for <- \(m) if (m == "Olink GAL" && "plate" %in% names(S)) "plate" else character()    # Olink values: also adjusted for plate
Zof <- \(cv) { cv <- intersect(cv, names(S)); if (length(cv)) S[cv] }

# ---- a) galanin vs lab lipids ---------------------------------------------------------------------------------------------------
sets <- list(`all persons` = character(), `adjusted for sex` = "sex_male",
             `adjusted for sex, age, BMI` = c("sex_male", unlist(h$adjust_extra %||% c("age", "BMI"))))
if (!has_sex) sets <- sets[1]
lip <- imap(measures, \(v, m) map(lipids, \(lp) {
  rows <- imap(sets, \(cv, nm) spearman1(v, S[[lp]], Zof(c(cv, if (length(cv)) plate_for(m))), min_n = 8) |> mutate(analysis = nm))
  if (has_sex) rows <- c(rows, list(
    spearman1(v[S$sex_male %in% 0], S[[lp]][S$sex_male %in% 0], min_n = 8) |> mutate(analysis = "women only"),
    spearman1(v[S$sex_male %in% 1], S[[lp]][S$sex_male %in% 1], min_n = 8) |> mutate(analysis = "men only")))
  bind_rows(rows) |> mutate(measure = m, lipid = lp)
}) |> bind_rows()) |> bind_rows() |>
  mutate(label = param_label(lipid), analysis = factor(analysis, c(names(sets), "women only", "men only"))) |>
  relocate(measure, lipid, label, analysis)
ga <- \(m, lp, an) lip |> filter(measure == m, lipid == lp, analysis == an)

# effect size: linear model of galanin (log2 ELISA or NPX) on the lipid (per SD), with covariates
lm_tab <- imap(measures, \(v, m) map(intersect(c(hdl_col, apo_col), names(S)), \(lp) map(names(sets), \(nm) {
  cv <- intersect(c(sets[[nm]], if (length(sets[[nm]])) plate_for(m)), names(S))
  x <- S |> mutate(y = v, lipid_sd = zscore(.data[[lp]])) |> select(y, lipid_sd, all_of(cv)) |> na.omit()
  f <- lm(reformulate(c("lipid_sd", cv), "y"), data = x)
  ci <- suppressWarnings(confint(f)["lipid_sd", ]); b <- coef(f)[["lipid_sd"]]
  tibble(measure = m, lipid = lp, label = param_label(lp), model = nm, n = nrow(x), estimate_per_SD = b, ci_low = ci[[1]], ci_high = ci[[2]],
         p = summary(f)$coefficients["lipid_sd", 4],
         meaning = if (m == "galanin ELISA") sprintf("%+.0f%% galanin per SD (95%% CI %+.0f to %+.0f%%)", 100 * (2^b - 1), 100 * (2^ci[[1]] - 1), 100 * (2^ci[[2]] - 1))
                   else sprintf("%+.2f NPX per SD", b))
}) |> bind_rows()) |> bind_rows()) |> bind_rows()

# ---- b) specificity for HDL ---------------------------------------------------------------------------------------------------------------
comps <- keep(h$comparisons %||% list(c("C_HDL", "C_LDL"), c("C_HDL", "C_TRIGLY"), c("c_apo", "C_APO_B"), c("c_apo", "C_HDL")),
              \(cp) all(unlist(cp) %in% names(S)))
spec <- imap(measures, \(v, m) map(comps, \(cp) {
  cp <- unlist(cp)
  bind_rows(boot_rho_diff(v, S[[cp[1]]], S[[cp[2]]], B = B) |> mutate(analysis = "all persons"),
            if (has_sex) boot_rho_diff(v, S[[cp[1]]], S[[cp[2]]], Zof(c("sex_male", plate_for(m))), B = B) |> mutate(analysis = "adjusted for sex")) |>
    mutate(measure = m, comparison = sprintf("rho(%s) - rho(%s)", param_label(cp[1]), param_label(cp[2])))
}) |> bind_rows()) |> bind_rows() |> relocate(measure, comparison, analysis)

# ---- c) HDL proteins measured by Olink ------------------------------------------------------------------------------------------------------
hp <- bind_rows(tibble(gene = unlist(h$hdl_proteins), group = "HDL particle"),
                tibble(gene = unlist(h$other_lipoprotein_proteins), group = "other lipoproteins")) |>
  mutate(OlinkID = map_chr(gene, \(g) find_assay(d$det, g)[1])) |> distinct(gene, .keep_all = TRUE)
hp_found <- hp |> filter(!is.na(OlinkID)) |> left_join(d$det |> select(OlinkID, frac_above_lod, measurable), by = "OlinkID")
score_ids <- hp_found |> filter(gene %in% unlist(h$score_proteins %||% c("APOA1", "APOA2", "APOM", "LCAT", "PON1", "PON3")), measurable) |> pull(OlinkID)
hdl_score <- if (length(score_ids) >= 2) rowMeans(apply(d$Y[, score_ids, drop = FALSE], 2, zscore), na.rm = TRUE) else NULL
Xp <- cbind(d$Y[, hp_found$OlinkID, drop = FALSE], if (!is.null(hdl_score)) cbind(`HDL protein score` = hdl_score))
labs_p <- c(hp_found$gene, if (!is.null(hdl_score)) "HDL protein score")
grp_p  <- c(hp_found$group, if (!is.null(hdl_score)) "score")
hp_tab <- if (ncol(Xp)) {
  valid <- bind_rows(
    spearman_vs(S[[hdl_col]], Xp, min_n = 8) |> mutate(with = "lab HDL-C"),
    if (apo_col %in% names(S)) spearman_vs(S[[apo_col]], Xp, min_n = 8) |> mutate(with = "lab ApoA-I")) |>
    mutate(protein = labs_p[match(id, colnames(Xp))]) |> select(protein, with, rho) |>
    pivot_wider(names_from = with, values_from = rho, names_prefix = "rho with ")
  imap(measures, \(v, m) bind_rows(
    spearman_vs(v, Xp, min_n = 8) |> mutate(analysis = "all persons"),
    if (has_sex) spearman_vs(v, Xp, Zof(c("sex_male", plate_for(m))), min_n = 8) |> mutate(analysis = "adjusted for sex")) |>
      mutate(measure = m)) |> bind_rows() |>
    mutate(protein = labs_p[match(id, colnames(Xp))], group = grp_p[match(id, colnames(Xp))]) |>
    left_join(valid, by = "protein") |>
    left_join(hp_found |> select(OlinkID, frac_above_lod), by = c("id" = "OlinkID")) |>
    select(measure, analysis, protein, group, OlinkID = id, n, rho, ci_low, ci_high, p, starts_with("rho with"), frac_above_lod)
} else tibble()

# ---- d) HDL-C vs all proteins: does GAL behave like an HDL-associated protein? -------------------------------------------------------------
Yh <- if (d$has_gal && !d$gal %in% colnames(d$M)) cbind(d$M, d$Y[, d$gal, drop = FALSE]) else d$M
hv <- bind_rows(spearman_vs(S[[hdl_col]], Yh) |> mutate(analysis = "all persons"),
                if (has_sex) spearman_vs(S[[hdl_col]], Yh, Zof(c("sex_male", "plate"))) |> mutate(analysis = "adjusted for sex, Olink plate")) |>
  group_by(analysis) |> mutate(fdr = p.adjust(p, "BH"), rank_by_rho = rank(-rho, ties.method = "min", na.last = "keep")) |> ungroup() |>
  rename(OlinkID = id) |> left_join(d$det |> select(OlinkID, protein), by = "OlinkID") |>
  mutate(is_galanin = OlinkID %in% d$gal, hdl_protein = OlinkID %in% hp_found$OlinkID[hp_found$group == "HDL particle"]) |>
  relocate(analysis, OlinkID, protein) |> arrange(analysis, p)
gal_hdl_rank <- hv |> filter(is_galanin)

# ---- e) do the two assays disagree more when HDL is high? -----------------------------------------------------------------------------------
disc <- NULL; disc_tab <- tibble(); split_tab <- tibble()
if (d$has_gal && d$has_elisa) {
  ok <- !is.na(d$npx) & coalesce(d$elisa > 0, FALSE)
  disc <- rep(NA_real_, nrow(S)); disc[ok] <- zscore(log2(d$elisa[ok])) - zscore(d$npx[ok])
  targets <- compact(list(`HDL cholesterol` = S[[hdl_col]], `ApoA-I (lab)` = if (apo_col %in% names(S)) S[[apo_col]], `HDL protein score` = hdl_score))
  disc_tab <- imap(targets, \(t, nm) bind_rows(spearman1(disc, t, min_n = 8) |> mutate(analysis = "all persons"),
                                               if (has_sex) spearman1(disc, t, Zof(c("sex_male", "plate")), min_n = 8) |> mutate(analysis = "adjusted for sex, Olink plate")) |>
                     mutate(with = nm)) |> bind_rows() |> relocate(with, analysis)
  # ELISA-Olink agreement below and above the median HDL; bootstrap over persons for the difference
  hi  <- S[[hdl_col]] > median(S[[hdl_col]], na.rm = TRUE)
  idx <- which(ok & !is.na(hi))
  rho_in <- \(x, grp) suppressWarnings(cor(d$elisa[x][grp], d$npx[x][grp], method = "spearman"))
  low <- rho_in(idx, !hi[idx]); high <- rho_in(idx, hi[idx])
  bs <- replicate(B, { x <- sample(idx, replace = TRUE); rho_in(x, hi[x]) - rho_in(x, !hi[x]) })
  bs <- bs[is.finite(bs)]
  split_tab <- tibble(HDL = c("below the median", "above the median", "difference (above - below)"),
                      n = c(sum(!hi[idx]), sum(hi[idx]), NA),
                      rho_ELISA_vs_Olink = c(low, high, high - low),
                      ci_low = c(NA, NA, unname(quantile(bs, 0.025))), ci_high = c(NA, NA, unname(quantile(bs, 0.975))))
}

# ---- answers --------------------------------------------------------------------------------------------------------------------------------------
for (m in names(measures)) {
  a0 <- ga(m, hdl_col, "all persons"); a1 <- ga(m, hdl_col, "adjusted for sex"); a2 <- ga(m, hdl_col, "adjusted for sex, age, BMI")
  w <- ga(m, hdl_col, "women only"); mn <- ga(m, hdl_col, "men only")
  eff <- lm_tab |> filter(measure == m, lipid == hdl_col, model == if (has_sex) "adjusted for sex" else "all persons")
  sig0 <- isTRUE(a0$p < 0.05); sig1 <- isTRUE(a1$p < 0.05); pos0 <- isTRUE(a0$rho > 0); pos1 <- isTRUE(a1$rho > 0)
  verdict <- if (!nrow(a0) || is.na(a0$p)) "not enough data"
    else if (sig0 && pos0 && has_sex && sig1 && pos1) "yes, positive - also within sex (not explained by sex)"
    else if (sig0 && pos0 && has_sex) "yes, positive - but largely explained by sex (women have higher HDL and higher galanin)"
    else if (sig0 && pos0) "yes, positive"
    else if (sig0) "negative correlation"
    else if (has_sex && sig1 && pos1) "yes, positive after adjusting for sex"
    else "no clear correlation"
  answer(ans, Q, sprintf("Does the %s correlate with HDL cholesterol?", m), verdict,
         paste0(sprintf("All persons rho %+.2f (p = %s, n = %d)", a0$rho, fmt_p(a0$p), a0$n),
                if (has_sex) sprintf("; adjusted for sex%s %+.2f (p = %s); sex, age, BMI %+.2f (p = %s); women %+.2f (p = %s, n = %d); men %+.2f (p = %s, n = %d)",
                                     if (m == "Olink GAL") " (and plate)" else "", a1$rho, fmt_p(a1$p), a2$rho, fmt_p(a2$p),
                                     w$rho, fmt_p(w$p), w$n, mn$rho, fmt_p(mn$p), mn$n) else "",
                if (nrow(eff)) sprintf(". Linear model (%s): %s, p = %s.", eff$model, eff$meaning, fmt_p(eff$p)) else "."))
  oth <- lip |> filter(measure == m, analysis == "all persons", lipid != hdl_col)
  sp <- spec |> filter(measure == m)
  answer(ans, Q, sprintf("Is the %s correlation specific to HDL?", m),
         if (!nrow(oth)) "no other lipids" else if (all(abs(oth$rho[oth$lipid != apo_col]) < abs(a0$rho), na.rm = TRUE) && isTRUE(a0$p < 0.05))
           "yes: HDL is the lipid most closely related to galanin" else "not clearly: other lipids relate as closely",
         paste0(paste(sprintf("%s %+.2f (p = %s)", oth$label, oth$rho, fmt_p(oth$p)), collapse = "; "),
                if (nrow(sp)) paste0(". Bootstrap differences: ", paste(sprintf("%s [%s] %+.2f (95%% CI %+.2f to %+.2f)", sp$comparison, sp$analysis,
                                                                          sp$difference, sp$ci_low, sp$ci_high), collapse = "; ")) else "", "."))
}
if (nrow(hp_tab)) {
  sc <- hp_tab |> filter(protein == "HDL protein score")
  hdlp <- hp_tab |> filter(group == "HDL particle", analysis == if (has_sex) "adjusted for sex" else "all persons") |> arrange(p)
  answer(ans, Q, "Olink HDL proteins: does galanin follow them?",
         if (nrow(sc)) paste(sprintf("%s vs HDL protein score: rho %+.2f (p = %s, %s)", sc$measure, sc$rho, fmt_p(sc$p), sc$analysis), collapse = "; ")
         else "no HDL protein score (fewer than 2 score proteins measured)",
         paste0(if (nrow(sc)) sprintf("The score (mean z of %s) tracks the lab HDL-C with rho %.2f. ", paste(hp_found$gene[hp_found$OlinkID %in% score_ids], collapse = ", "),
                                      sc$`rho with lab HDL-C`[1]) else "",
                "Strongest HDL proteins (adjusted for sex): ", paste(head(sprintf("%s ~ %s %+.2f (p = %s)", hdlp$measure, hdlp$protein, hdlp$rho, fmt_p(hdlp$p)), 10), collapse = "; "),
                if (length(setdiff(hp$gene, hp_found$gene))) sprintf(". Not on the Olink panel: %s.", paste(setdiff(hp$gene, hp_found$gene), collapse = ", ")) else "."))
}
if (nrow(gal_hdl_rank)) {
  g0 <- gal_hdl_rank |> filter(analysis == "all persons"); n_hv <- sum(!is.na(hv$rho[hv$analysis == "all persons"]))
  top_h <- hv |> filter(analysis == "all persons")
  answer(ans, Q, "Does Olink GAL behave like an HDL-associated protein?",
         if (isTRUE(g0$rank_by_rho / n_hv <= 0.05)) sprintf("yes: GAL is in the top 5%% of proteins correlated with HDL-C (rank %d of %d)", g0$rank_by_rho, n_hv)
         else sprintf("no: GAL ranks %d of %d proteins by correlation with HDL-C", g0$rank_by_rho, n_hv),
         sprintf("HDL-C vs Olink GAL rho %+.2f (p = %s). The proteins most correlated with HDL-C: %s.", g0$rho, fmt_p(g0$p), top_str(top_h$protein, top_h$rho, top_h$p, 8)))
}
if (nrow(disc_tab)) {
  d0 <- disc_tab |> filter(with == "HDL cholesterol", analysis == "all persons")
  d1 <- disc_tab |> filter(with == "HDL cholesterol", analysis != "all persons")
  answer(ans, Q, "Do the ELISA and Olink disagree more when HDL is high?",
         case_when(is.na(d0$p) ~ "not enough data", d0$p < 0.05 & d0$rho > 0 ~ "yes: the ELISA reads relatively higher than Olink when HDL is high",
                   d0$p < 0.05 ~ "yes: the ELISA reads relatively lower than Olink when HDL is high", TRUE ~ "no clear relation"),
         paste0(sprintf("Discordance = z(ELISA) - z(Olink GAL) vs HDL-C: rho %+.2f (p = %s)", d0$rho, fmt_p(d0$p)),
                if (nrow(d1)) sprintf("; %s: %+.2f (p = %s)", d1$analysis, d1$rho, fmt_p(d1$p)) else "",
                if (nrow(split_tab)) sprintf(". ELISA vs Olink: rho %.2f below and %.2f above the median HDL (difference %+.2f, 95%% CI %+.2f to %+.2f)",
                                             split_tab$rho_ELISA_vs_Olink[1], split_tab$rho_ELISA_vs_Olink[2], split_tab$rho_ELISA_vs_Olink[3],
                                             split_tab$ci_low[3], split_tab$ci_high[3]) else "",
                ". If galanin binds HDL and one assay does not detect the bound form (e.g. a hidden epitope), the disagreement grows with HDL."))
}
answer(ans, Q, "What would confirm binding to HDL?", "needs experiments - correlations cannot show binding",
       "Measure galanin in lipoprotein fractions (ultracentrifugation or size-exclusion chromatography: HDL vs LDL/VLDL vs lipoprotein-free); precipitate ApoA-I (anti-ApoA-I beads) and measure the co-precipitated galanin; compare the ELISA and Olink in HDL-rich vs HDL-depleted serum or after delipidation (spike-in recovery).")
answers <- answers_save(ans, cfg, "3_galanin_hdl", "answers.csv")

# ---- figures -----------------------------------------------------------------------------------------------------------------------------------------
F <- figs_new()
sexlab <- if (has_sex) if_else(S$sex_male == 1, "men", "women") else rep("all", nrow(S))
sc_d <- imap(measures, \(v, m) tibble(measure = m, hdl = S[[hdl_col]], y = v, sex = sexlab)) |> bind_rows() |> filter(!is.na(hdl), !is.na(y), !is.na(sex)) |>
  mutate(measure = factor(measure, names(measures), c(`galanin ELISA` = "galanin ELISA (log2 pg/mL)", `Olink GAL` = "Olink GAL (NPX)")[names(measures)]))
sub <- map_chr(names(measures), \(m) { a0 <- ga(m, hdl_col, "all persons"); a1 <- ga(m, hdl_col, "adjusted for sex")
  sprintf("%s: rho %+.2f (p = %s)%s", m, a0$rho, fmt_p(a0$p), if (nrow(a1)) sprintf(", adjusted for sex %+.2f (p = %s)", a1$rho, fmt_p(a1$p)) else "") })
p1 <- ggplot(sc_d, aes(hdl, y)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey30", linetype = 2, linewidth = 0.5) +
  geom_smooth(aes(colour = sex), method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
  geom_point(aes(colour = sex), size = 2.2) + facet_wrap(~measure, scales = "free_y") +
  scale_colour_manual(values = c(women = "#c0392b", men = "#2471a3", all = "grey30")) +
  labs(title = "Galanin vs HDL cholesterol (LEIP)", subtitle = paste(c(sub, "dashed: all persons; coloured lines: within sex"), collapse = "\n"),
       x = "HDL cholesterol (lab)", y = NULL, colour = NULL)
fig(F, "hdl_by_sex", p1, cfg, "3_galanin_hdl", "galanin_vs_HDL_by_sex.png", width = 11, height = 5.5)

an_lev <- levels(lip$analysis)
lp_d <- lip |> mutate(lab = factor(label, rev(unique(param_label(lipids)))),
                      y = as.numeric(lab) + (as.numeric(analysis) - (length(an_lev) + 1) / 2) * 0.14, sig = coalesce(p < 0.05, FALSE))
p2 <- ggplot(lp_d, aes(rho, y, colour = analysis)) + geom_vline(xintercept = 0, colour = "grey60") +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), alpha = 0.4) + geom_point(aes(shape = sig), size = 1.9) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
  scale_colour_manual(values = setNames(c("black", "#8e44ad", "#16a085", "#c0392b", "#2471a3")[seq_along(an_lev)], an_lev), name = NULL) +
  scale_y_continuous(breaks = seq_along(levels(lp_d$lab)), labels = levels(lp_d$lab)) + facet_wrap(~measure) +
  labs(title = "Galanin vs the lab lipids: is it HDL, and is it sex?", subtitle = "Spearman / partial Spearman with 95% CI (Olink GAL: adjusted models also for plate)",
       x = "rho", y = NULL) + theme(legend.position = "bottom")
fig(F, "lipids", p2, cfg, "3_galanin_hdl", "galanin_vs_lipids.png", width = 11, height = 6.5)

if (nrow(hp_tab)) {
  x <- hp_tab |> filter(analysis == if (has_sex) "adjusted for sex" else "all persons")
  lev <- x |> group_by(protein) |> summarise(v = mean(`rho with lab HDL-C`, na.rm = TRUE)) |> arrange(v) |> pull(protein)
  x <- x |> mutate(y = match(protein, lev) + if_else(measure == "Olink GAL", 0.15, -0.15), sig = coalesce(p < 0.05, FALSE),
                   ylab = sprintf("%s (HDL-C rho %.2f)", protein, `rho with lab HDL-C`))
  p3 <- ggplot(x, aes(rho, y, colour = measure)) + geom_vline(xintercept = 0, colour = "grey60") +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), alpha = 0.45) + geom_point(aes(shape = sig), size = 2.1) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"), name = NULL) +
    scale_colour_manual(values = c(`Olink GAL` = "firebrick", `galanin ELISA` = "steelblue"), name = NULL) +
    scale_y_continuous(breaks = seq_along(lev), labels = x$ylab[match(lev, x$protein)], expand = expansion(add = 0.6)) +
    labs(title = "Galanin vs lipoprotein proteins measured by Olink", subtitle = sprintf("rho %s; in brackets: how closely the protein tracks the lab HDL-C",
                                                                                       if (has_sex) "adjusted for sex (Olink GAL also for plate)" else "(unadjusted)"),
         x = "rho with galanin", y = NULL) + theme(legend.position = "bottom")
  fig(F, "hdl_proteins", p3, cfg, "3_galanin_hdl", "galanin_vs_HDL_proteins.png", width = 9, height = 2.5 + 0.3 * length(lev))
}
v0 <- hv |> filter(analysis == "all persons", !is.na(p))
p4 <- ggplot(v0, aes(rho, -log10(p))) + geom_point(colour = "grey65", size = 0.8, alpha = 0.7) +
  geom_point(data = v0 |> filter(hdl_protein), colour = "#16a085", size = 1.8) +
  geom_text(data = bind_rows(head(v0, 10), v0 |> filter(hdl_protein)) |> distinct(OlinkID, .keep_all = TRUE), aes(label = protein),
            size = 2.6, vjust = -0.6, check_overlap = TRUE) +
  geom_point(data = v0 |> filter(is_galanin), colour = "firebrick", size = 3) +
  geom_text(data = v0 |> filter(is_galanin), aes(label = protein), colour = "firebrick", vjust = 1.8, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  labs(title = "Which Olink proteins track HDL cholesterol? (LEIP)",
       subtitle = sprintf("green: HDL-particle proteins; red: Olink GAL%s", if (nrow(gal_hdl_rank)) sprintf(" (rank %d)", gal_hdl_rank$rank_by_rho[gal_hdl_rank$analysis == "all persons"]) else ""),
       x = "Spearman rho with HDL-C", y = "-log10 p")
fig(F, "hdl_volcano", p4, cfg, "3_galanin_hdl", "HDL_vs_all_proteins.png", width = 9, height = 7)
if (!is.null(disc)) {
  dd <- tibble(hdl = S[[hdl_col]], disc = disc, sex = sexlab) |> filter(!is.na(hdl), !is.na(disc))
  d0 <- disc_tab |> filter(with == "HDL cholesterol", analysis == "all persons")
  p5 <- ggplot(dd, aes(hdl, disc)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", linewidth = 0.6) + geom_point(aes(colour = sex), size = 2.3) +
    scale_colour_manual(values = c(women = "#c0392b", men = "#2471a3", all = "grey30")) +
    labs(title = "Do the ELISA and Olink disagree more when HDL is high?",
         subtitle = sprintf("discordance = z(ELISA) - z(Olink GAL); rho with HDL-C %+.2f (p = %s)", d0$rho, fmt_p(d0$p)),
         x = "HDL cholesterol (lab)", y = "z(ELISA) - z(Olink GAL)", colour = NULL)
  fig(F, "discordance", p5, cfg, "3_galanin_hdl", "assay_discordance_vs_HDL.png", width = 8, height = 6)
}

# ---- tables --------------------------------------------------------------------------------------------------------------------------------------
save_csv(lip, cfg, "3_galanin_hdl", "galanin_vs_lipids.csv")
save_csv(hv, cfg, "3_galanin_hdl", "HDL_vs_all_proteins.csv")
save_csv(spec, cfg, "3_galanin_hdl", "HDL_specificity.csv")
save_csv(lm_tab, cfg, "3_galanin_hdl", "effect_sizes.csv")
if (nrow(hp_tab)) save_csv(hp_tab, cfg, "3_galanin_hdl", "lipoprotein_proteins.csv")
if (nrow(disc_tab)) save_csv(disc_tab, cfg, "3_galanin_hdl", "discordance_vs_HDL.csv")
figs_save(F, cfg, "3_galanin_hdl", "figures.rds")
writexl::write_xlsx(list(answers = answers, galanin_vs_lipids = lip, effect_sizes = lm_tab, HDL_specificity = spec,
                         lipoprotein_proteins = hp_tab, HDL_vs_all_proteins = hv, discordance_vs_HDL = disc_tab,
                         agreement_by_HDL = split_tab) |> keep(\(x) is.data.frame(x) && ncol(x) > 0),
                    out_path(cfg, "3_galanin_hdl", "galanin_hdl.xlsx"))
for (i in seq_len(nrow(answers))) msg("%s: %s", answers$item[i], answers$verdict[i])
