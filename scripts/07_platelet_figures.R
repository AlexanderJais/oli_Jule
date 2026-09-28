# 07 - Paper figures: is serum galanin released from platelets?
# Uses the Olink data of the LEIP sera only. Asks whether galanin (Olink GAL) rises and falls with the
# proteins that platelets release or shed when blood clots. The platelet set comes from the literature
# (config: platelets) - it is NOT chosen by correlation with galanin, so the test is not circular.
#   Figure 1  discovery: GAL in every sample; GAL vs all measurable proteins (platelet and neuroendocrine
#             proteins marked); platelet gene sets along that ranking (GSEA enrichment curves)
#   Figure 2  the platelet set: do its proteins co-vary (one release factor)? GAL vs the platelet release
#             score and vs each platelet protein; specificity against random protein sets of the same
#             size and against all other proteins
#   Figure 3  robustness (adjusted, within sex, per plate, leaving out a person or a protein), contrast with
#             the neuroendocrine proteins, and whether GAL's links to CRP and age depend on platelet release
#   S1-S3     the platelet proteins (detection, fit to the score); main axes of the serum proteome; lab
#             assays vs Olink in the same sera (data quality)
# Out: output/paper_platelets/  Figure*.pdf/.png/.tiff, figure_legends.md, source_data.xlsx, key_results.csv

source("R/utils.R")
source("R/leip.R")
source("R/figures.R")
cfg <- load_config()
set.seed(cfg$seed %||% 1)
d <- leip_load(cfg)
S <- d$S
out <- "paper_platelets"
clear_outputs(cfg, out)
if (!d$has_gal) stop("Galanin (GAL) is not in the Olink data - no figures.")
g <- d$npx
pc <- cfg$platelets %||% list()
has_sex <- "sex_male" %in% names(S) && n_distinct(na.omit(S$sex_male)) == 2
crp_col <- pc$crp_column %||% "c_CRP"
plate_short <- \(x) { s <- str_match(x, "(?i)plate[ _-]?0*(\\d+)")[, 2]; if_else(is.na(s), x, paste("Plate", s)) }
fmt_n <- \(x) format(x, big.mark = ",")
src <- list()                                                # source data, one sheet per panel
key <- tibble(item = character(), value = character())      # key numbers (key_results.csv)
kv <- \(item, value) key <<- bind_rows(key, tibble(item = item, value = as.character(value)))

# ---- protein sets ------------------------------------------------------------------------------------------------------------
pset <- bind_rows(tibble(gene = unlist(pc$released), group = "released"), tibble(gene = unlist(pc$membrane), group = "membrane")) |>
  distinct(gene, .keep_all = TRUE) |>
  mutate(OlinkID = map_chr(gene, \(x) find_assay(d$det, x)[1])) |>
  left_join(d$det |> select(OlinkID, frac_above_lod, measurable), by = "OlinkID") |>
  mutate(status = case_when(is.na(OlinkID) ~ "not on the panel", OlinkID %in% d$gal ~ "galanin itself",
                            !coalesce(measurable, FALSE) ~ "below LOD in most samples", TRUE ~ "used"))
plt <- pset |> filter(status == "used") |> distinct(OlinkID, .keep_all = TRUE)
if (nrow(plt) < 3) stop("Fewer than 3 proteins of the platelet set are measurable - check config: platelets.")
score_of <- \(ids) rowMeans(scale(d$Y[, ids, drop = FALSE]), na.rm = TRUE)       # mean z-score per sample
pscore <- score_of(plt$OlinkID)
nes <- tibble(gene = unlist(cfg$neuroendocrine_proteins)) |> mutate(OlinkID = map_chr(gene, \(x) find_assay(d$det, x)[1])) |>
  left_join(d$det |> select(OlinkID, measurable), by = "OlinkID") |>
  filter(!is.na(OlinkID), !OlinkID %in% d$gal, coalesce(measurable, FALSE)) |> distinct(OlinkID, .keep_all = TRUE)
nscore <- if (nrow(nes) >= 3) score_of(nes$OlinkID) else NULL
kv("platelet proteins used", sprintf("%d of %d (%s)", nrow(plt), nrow(pset), paste(plt$gene, collapse = ", ")))
kv("neuroendocrine proteins used", sprintf("%d (%s)", nrow(nes), paste(nes$gene, collapse = ", ")))

# ---- GAL vs every measurable protein ----------------------------------------------------------------------------------------
Mp <- d$M[, setdiff(colnames(d$M), d$gal), drop = FALSE]
u <- spearman_vs(g, Mp)
prot <- tibble(OlinkID = u$id, n = u$n, rho = u$rho, ci_low = u$ci_low, ci_high = u$ci_high, p = u$p) |>
  mutate(fdr = p.adjust(p, "BH")) |> left_join(d$det |> select(OlinkID, Assay, protein), by = "OlinkID") |>
  mutate(set = case_when(OlinkID %in% plt$OlinkID ~ "platelet set", OlinkID %in% nes$OlinkID ~ "neuroendocrine set", TRUE ~ "other")) |>
  arrange(p)
m_prot <- sum(!is.na(prot$p))

# gene-set enrichment along the ranking (GSEA): the pre-specified platelet gene sets
st <- prot |> filter(!is.na(rho)) |> transmute(gene = str_split(Assay, "_"), stat = rho) |> unnest(gene) |>
  group_by(gene) |> slice_max(abs(stat), n = 1, with_ties = FALSE) |> ungroup()
stats <- sort(setNames(st$stat, st$gene), decreasing = TRUE)
gs_want <- unlist(pc$gene_sets %||% c("GOBP_PLATELET_ACTIVATION", "REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2"))
gsets <- tryCatch(load_gene_sets(unique(c("C5:GO:BP", "C2:CP:REACTOME"))), error = \(e) { msg("Gene sets not available: %s", conditionMessage(e)); list() })
gsets <- keep(gsets[intersect(gs_want, names(gsets))], \(s) sum(names(stats) %in% s) >= 5)
pretty_set <- \(x) recode(x, GOBP_PLATELET_ACTIVATION = "GO: platelet activation",
                          REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2 = "Reactome: platelet degranulation",
                          .default = str_to_sentence(str_replace_all(str_remove(x, "^(GOBP|REACTOME|HALLMARK)_"), "_", " ")))
gsea <- if (length(gsets)) as_tibble(suppressWarnings(fgsea::fgseaMultilevel(gsets, stats, minSize = 5, maxSize = 5000, eps = 0))) else tibble()
curves <- imap(gsets, \(s, nm) { hit <- names(stats) %in% s; w <- abs(stats) * hit
  tibble(set = nm, rank = seq_along(stats), es = cumsum(w) / sum(w) - cumsum(!hit) / sum(!hit), hit = hit) }) |> bind_rows()

# ---- galanin vs the platelet release score ----------------------------------------------------------------------------------
r_all <- spearman1(pscore, g)
r_adj <- spearman1(pscore, g, S[d$cov])
ok <- !is.na(g) & !is.na(pscore)
rg <- rank(g[ok])
pool <- setdiff(colnames(d$M), c(d$gal, plt$OlinkID, nes$OlinkID))
Zp <- scale(d$M[ok, pool, drop = FALSE])
B <- pc$random_sets %||% 10000
null_of <- \(k) replicate(B, cor(rank(rowMeans(Zp[, sample.int(ncol(Zp), k), drop = FALSE], na.rm = TRUE)), rg))
null_p <- null_of(nrow(plt))
p_rand <- (1 + sum(null_p >= r_all$rho)) / (B + 1)                      # one-sided: the hypothesis is a positive link
r_ne <- if (!is.null(nscore)) spearman1(nscore, g) else tibble()
null_n <- if (!is.null(nscore)) null_of(nrow(nes)) else numeric()
p_rand_ne <- if (length(null_n)) (1 + sum(null_n >= r_ne$rho)) / (B + 1) else NA
Yg <- if (d$gal %in% colnames(d$M)) d$M else cbind(d$M, d$Y[, d$gal, drop = FALSE])
allr <- spearman_vs(pscore, Yg[, setdiff(colnames(Yg), plt$OlinkID), drop = FALSE]) |> filter(!is.na(rho)) |>
  left_join(d$det |> select(id = OlinkID, protein), by = "id") |> mutate(rank = rank(-rho, ties.method = "min")) |> arrange(rank)
gal_rank <- allr$rank[allr$id == d$gal]; n_allr <- nrow(allr)
kv("GAL vs platelet score", rho_txt(r_all)); kv("GAL vs platelet score, adjusted for age, sex, plate", rho_txt(r_adj))
kv("random protein sets of the same size: one-sided p", fmt_p(p_rand))
kv("GAL rank among all proteins by correlation with the platelet score", sprintf("%d of %d (top %.1f%%)", gal_rank, n_allr, 100 * gal_rank / n_allr))
if (nrow(r_ne)) kv("GAL vs neuroendocrine score", rho_txt(r_ne))

# robustness
loo <- \(x, y) { k <- which(!is.na(x) & !is.na(y)); map_dbl(k, \(i) { j <- setdiff(k, i); cor(x[j], y[j], method = "spearman") }) }
flagged <- coalesce(col_or(S, "SampleQC") != "PASS", FALSE) | coalesce(as.logical(col_or(S, "qc_outlier", FALSE)), FALSE) |
  coalesce(col_or(S, "low_serum", 0) == 1, FALSE) | coalesce(col_or(S, "lipamisch", 0) == 1, FALSE)
row_of <- \(what, r, type = "estimate") tibble(analysis = what, n = r$n, rho = r$rho, lo = r$ci_low, hi = r$ci_high, p = r$p, type = type)
sub_r <- \(keep) spearman1(pscore[keep], g[keep], min_n = 5)
lr <- loo(pscore, g)
jk <- map_dbl(plt$OlinkID, \(id) cor(score_of(setdiff(plt$OlinkID, id))[ok], g[ok], method = "spearman"))
plates <- sort(unique(na.omit(S$plate)))
rob <- bind_rows(
  row_of(sprintf("all sera (n = %d)", r_all$n), r_all),
  row_of(sprintf("adjusted for %s", covs_label(d$cov)), r_adj),
  if (crp_col %in% names(S)) row_of(sprintf("adjusted for %s, CRP", covs_label(d$cov)), spearman1(pscore, g, S[c(d$cov, crp_col)])),
  if (has_sex) row_of("women only", sub_r(S$sex_male %in% 0)),
  if (has_sex) row_of("men only", sub_r(S$sex_male %in% 1)),
  if (any(flagged)) row_of(sprintf("without flagged sera (n = %d)", sum(!flagged & ok)), sub_r(!flagged)),
  map(plates, \(pl) row_of(sprintf("%s only", plate_short(pl)), sub_r(S$plate %in% pl))),
  tibble(analysis = "leaving out one serum (range)", n = r_all$n - 1L, rho = r_all$rho, lo = min(lr), hi = max(lr), type = "range"),
  tibble(analysis = "leaving out one platelet protein (range)", n = r_all$n, rho = r_all$rho, lo = min(jk), hi = max(jk), type = "range"),
  map(c("released", "membrane"), \(grp) { ids <- plt$OlinkID[plt$group == grp]
    if (length(ids) >= 2) row_of(sprintf("%s proteins only (%d)", if (grp == "released") "released" else "membrane", length(ids)), spearman1(score_of(ids), g)) })) |>
  mutate(analysis = factor(analysis, rev(unique(analysis))))

# GAL's clinical links (CRP, age): do they depend on platelet release?
clin <- map(intersect(c(crp_col, "age"), names(S)), \(pm) {
  z0 <- setdiff(d$cov, pm); lab <- param_label(pm)
  bind_rows(row_of("galanin", spearman1(S[[pm]], g, min_n = 10)),
            row_of("galanin, adjusted", spearman1(S[[pm]], g, S[z0], min_n = 10)),
            row_of("galanin, adjusted + platelet score", spearman1(S[[pm]], g, bind_cols(S[z0], tibble(platelet_score = pscore)), min_n = 10)),
            row_of("platelet score", spearman1(S[[pm]], pscore, min_n = 10))) |> mutate(parameter = lab, adjusted_for = covs_label(z0))
}) |> bind_rows()

# ---- supplementary data: item-rest fit, proteome axes, lab vs Olink --------------------------------------------------------
fit <- map(plt$OlinkID, \(id) spearman1(d$Y[, id], score_of(setdiff(plt$OlinkID, id))) |> mutate(OlinkID = id)) |> bind_rows() |>
  left_join(plt |> select(OlinkID, gene, group), by = "OlinkID")
X <- apply(Mp, 2, \(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v }); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
pca <- prcomp(X, center = TRUE, scale. = TRUE); kp <- min(8, ncol(pca$x)); ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
axes <- map(seq_len(kp), \(i) bind_rows(spearman1(pca$x[, i], g) |> mutate(what = "galanin (GAL)"),
                                       spearman1(pca$x[, i], pscore) |> mutate(what = "platelet score")) |>
              mutate(component = sprintf("PC%d (%.0f%%)", i, ve[i]))) |> bind_rows()
bench <- imap(cfg$lab_vs_olink %||% list(), \(assay, lab) {
  if (lab == "galanin_elisa" || !lab %in% names(S) || sum(!is.na(S[[lab]])) < 10) return(NULL)
  oid <- find_assay(d$det, assay)[1]; if (is.na(oid)) return(NULL)
  ar <- spearman_vs(S[[lab]], if (oid %in% colnames(d$M)) d$M else cbind(d$M, d$Y[, oid, drop = FALSE]))
  spearman1(S[[lab]], d$Y[, oid]) |> mutate(lab = param_label(lab), olink = assay, rank = rank(-ar$rho, ties.method = "min", na.last = "keep")[ar$id == oid],
                                           of = sum(!is.na(ar$rho)))
}) |> compact() |> bind_rows()

# ---- Figure 1: discovery -----------------------------------------------------------------------------------------------------
gl <- d$lod |> filter(OlinkID == d$gal) |> select(SampleID, LOD, below_lod)
f1a <- tibble(SampleID = S$SampleID, plate = plate_short(S$plate), npx = g) |> left_join(gl, by = "SampleID") |> filter(!is.na(npx))
lodp <- f1a |> group_by(plate) |> summarise(LOD = median(LOD, na.rm = TRUE), .groups = "drop")
kw <- if (n_distinct(f1a$plate) > 1) kruskal.test(f1a$npx, factor(f1a$plate))$p.value else NA
n_above <- sum(!f1a$below_lod, na.rm = TRUE)
pA <- ggplot(f1a, aes(plate, npx)) +
  geom_errorbar(data = lodp, aes(x = plate, ymin = LOD, ymax = LOD), inherit.aes = FALSE, width = 0.6, linetype = 2, colour = fig_col[["muted"]], linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.2, colour = fig_col[["ink"]]) +
  labs(title = "Galanin in serum", subtitle = sprintf("%d of %d sera above LOD (dashed)", n_above, nrow(f1a)),
       x = NULL, y = "GAL (NPX)") + theme_pub()
vd <- prot |> filter(!is.na(p)) |> mutate(set = factor(set, c("other", "neuroendocrine set", "platelet set")))
lab1 <- bind_rows(vd |> filter(set == "platelet set", p < 0.05), vd |> slice_min(p, n = 8)) |> distinct(OlinkID, .keep_all = TRUE)
thr <- -log10(0.05 / m_prot)
pB <- ggplot(vd, aes(rho, -log10(p))) +
  geom_hline(yintercept = thr, linetype = 2, colour = fig_col[["muted"]], linewidth = 0.3) +
  annotate("text", x = min(vd$rho, na.rm = TRUE), y = thr, label = "Bonferroni 5%", hjust = 0, vjust = -0.5, size = 2, colour = fig_col[["muted"]]) +
  geom_point(data = vd |> filter(set == "other"), colour = fig_col[["light"]], size = 0.5) +
  geom_point(data = vd |> filter(set != "other"), aes(colour = set), size = 1.4) +
  ggrepel::geom_text_repel(data = lab1, aes(label = protein), size = 2, colour = fig_col[["ink"]], segment.size = 0.2,
                           segment.colour = fig_col[["muted"]], min.segment.length = 0, max.overlaps = 40, seed = 1) +
  scale_colour_manual(values = c(`platelet set` = fig_col[["platelet"]], `neuroendocrine set` = fig_col[["neuro"]]), name = NULL, drop = FALSE,
                      breaks = c("platelet set", "neuroendocrine set")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = sprintf("Galanin vs %s serum proteins", fmt_n(m_prot)), x = "Spearman rho with GAL", y = expression(-log[10]~italic(p))) +
  theme_pub() + theme(legend.position = "bottom")
pC <- if (nrow(curves)) {
  cv <- curves |> mutate(set = pretty_set(set))
  rows <- tibble(set = unique(cv$set), y0 = min(cv$es) - 0.07 * seq_along(unique(cv$set)))
  ticks <- cv |> filter(hit) |> left_join(rows, by = "set")
  ann <- gsea |> transmute(set = pretty_set(pathway), lab = sprintf("%s: NES %.2f, p = %s", set, NES, fmt_p(pval)))
  ggplot(cv, aes(rank, es, colour = set)) + geom_hline(yintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
    geom_line(linewidth = 0.6) +
    geom_segment(data = ticks, aes(x = rank, xend = rank, y = y0, yend = y0 + 0.05), linewidth = 0.3) +
    annotate("text", x = length(stats), y = max(cv$es), label = paste(ann$lab, collapse = "\n"), hjust = 1, vjust = 1, size = 2.1, colour = fig_col[["ink"]]) +
    scale_colour_manual(values = setNames(c(fig_col[["platelet"]], fig_col[["third"]])[seq_along(unique(cv$set))], unique(cv$set)), name = NULL) +
    labs(title = "Platelet gene sets along the galanin ranking (GSEA)", x = sprintf("%s proteins ranked by correlation with GAL (positive to negative)", fmt_n(length(stats))),
         y = "running enrichment score") + theme_pub() + theme(legend.position = "bottom")
} else ggplot() + labs(title = "Platelet gene sets: not available") + theme_void()
fig1 <- compose(list(compose(list(pA, pB), ncol = 2, rel_widths = c(0.7, 1.3), labels = c("A", "B")), compose(list(pC), ncol = 1, labels = "C")),
                ncol = 1, rel_heights = c(1.15, 0.85), labels = NULL)
save_figure(fig1, cfg, out, "Figure1", 180, 165)
src$Fig1A <- f1a; src$Fig1B <- prot; src$Fig1C_gsea <- gsea |> mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

# ---- Figure 2: the platelet set -----------------------------------------------------------------------------------------------
H <- cbind(g, d$Y[, plt$OlinkID, drop = FALSE]); colnames(H) <- c("GAL (galanin)", plt$gene)
R <- suppressWarnings(cor(H, method = "spearman", use = "pairwise.complete.obs"))
pair <- R[-1, -1][upper.tri(R[-1, -1])]
ord <- colnames(R)[hclust(as.dist(1 - R))$order]
hm <- as_tibble(R, rownames = "a") |> pivot_longer(-a, names_to = "b", values_to = "rho") |> mutate(a = factor(a, ord), b = factor(b, ord))
p2A <- ggplot(hm, aes(a, b, fill = rho)) + geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = fig_col[["low"]], mid = fig_col[["mid"]], high = fig_col[["high"]], limits = c(-1, 1), name = "rho") +
  coord_fixed() + labs(title = "Do the platelet proteins co-vary?", subtitle = sprintf("%d proteins; median pairwise rho %.2f", nrow(plt), median(pair)), x = NULL, y = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.line = element_blank(), axis.ticks = element_blank(),
                      panel.grid = element_blank(), legend.position = "right")
sc <- tibble(score = pscore, gal = g) |> filter(!is.na(score), !is.na(gal))
p2B <- ggplot(sc, aes(score, gal)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = fig_col[["ink2"]], fill = fig_col[["grid"]], linewidth = 0.5) +
  geom_point(colour = fig_col[["platelet"]], size = 1.4) +
  labs(title = "Galanin vs platelet release score", subtitle = sprintf("%s\nadjusted for %s: rho = %.2f, p = %s", rho_txt(r_all), covs_label(d$cov), r_adj$rho, fmt_p(r_adj$p)),
       x = sprintf("platelet release score (mean z of %d proteins)", nrow(plt)), y = "GAL (NPX)") + theme_pub()
fp <- spearman_vs(g, d$Y[, plt$OlinkID, drop = FALSE]) |> left_join(plt |> select(id = OlinkID, gene, group), by = "id") |>
  bind_rows(r_all |> mutate(gene = "platelet score", group = "score")) |>
  mutate(gene = factor(gene, c("platelet score", setdiff(gene[order(rho)], "platelet score"))),
         group = factor(group, c("released", "membrane", "score"), c("released", "membrane", "score")))
p2C <- ggplot(fp, aes(rho, gene, shape = group)) + geom_vline(xintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high, colour = group == "score"), linewidth = 0.4) +
  geom_point(aes(colour = group == "score"), size = 1.6) +
  scale_colour_manual(values = c(`FALSE` = fig_col[["platelet"]], `TRUE` = fig_col[["ink"]]), guide = "none") +
  scale_shape_manual(values = c(released = 16, membrane = 17, score = 18), name = NULL) +
  labs(title = "Galanin vs each platelet protein", x = "rho with GAL (95% CI)", y = NULL) + theme_pub() + row_guides() +
  theme(legend.position = "bottom")
nd <- bind_rows(tibble(rho = null_p, what = "random sets"))
p2D <- ggplot(nd, aes(rho)) + geom_histogram(bins = 50, fill = fig_col[["light"]], colour = NA) +
  geom_vline(xintercept = r_all$rho, colour = fig_col[["platelet"]], linewidth = 0.6) +
  annotate("text", x = r_all$rho, y = Inf, label = "platelet set", colour = fig_col[["ink"]], hjust = 1.08, vjust = 1.5, size = 2.1) +
  (if (nrow(r_ne)) list(geom_vline(xintercept = r_ne$rho, colour = fig_col[["neuro"]], linewidth = 0.6),
                        annotate("text", x = r_ne$rho, y = Inf, label = "neuroendocrine", colour = fig_col[["ink"]], hjust = 1.08, vjust = 3.2, size = 2.1))) +
  labs(title = "Specificity: random protein sets", subtitle = sprintf("%s random sets of %d proteins
one-sided p = %s", fmt_n(B), nrow(plt), fmt_p(p_rand)),
       x = "rho of the set score with GAL", y = "random sets") + theme_pub()
p2E <- ggplot(allr, aes(rho)) + geom_histogram(bins = 50, fill = fig_col[["light"]], colour = NA) +
  geom_vline(xintercept = allr$rho[allr$id == d$gal], colour = fig_col[["ink"]], linewidth = 0.6) +
  annotate("text", x = allr$rho[allr$id == d$gal], y = Inf, label = sprintf("GAL: rank %d of %s", gal_rank, fmt_n(n_allr)), hjust = 1.05, vjust = 1.5, size = 2.1) +
  labs(title = "Which proteins track platelet release?", subtitle = "proteins outside the platelet set",
       x = "rho with the platelet score", y = "proteins") + theme_pub()
fig2 <- compose(list(compose(list(p2A, p2B), ncol = 2, rel_widths = c(1.1, 0.9), labels = c("A", "B")),
                     compose(list(p2C, p2D, p2E), ncol = 3, labels = c("C", "D", "E"))), ncol = 1, rel_heights = c(1.05, 0.95), labels = NULL)
save_figure(fig2, cfg, out, "Figure2", 180, 190)
src$Fig2A <- as_tibble(R, rownames = "protein"); src$Fig2B <- sc; src$Fig2C <- fp
src$Fig2D <- tibble(random_set_rho = null_p); src$Fig2E <- allr

# ---- Figure 3: robustness, contrast, clinical links ------------------------------------------------------------------------
p3A <- ggplot(rob, aes(rho, analysis)) + geom_vline(xintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi, linetype = type), colour = fig_col[["platelet"]], linewidth = 0.45) +
  geom_point(colour = fig_col[["platelet"]], size = 1.5) +
  scale_linetype_manual(values = c(estimate = 1, range = 3), labels = c(estimate = "95% CI", range = "min to max"), name = NULL) +
  labs(title = "Galanin vs platelet score: robustness", x = "Spearman rho", y = NULL) + theme_pub() + row_guides() + theme(legend.position = "bottom")
p3B <- if (!is.null(nscore)) {
  nsd <- tibble(score = nscore, gal = g) |> filter(!is.na(score), !is.na(gal))
  ggplot(nsd, aes(score, gal)) + geom_smooth(method = "lm", formula = y ~ x, colour = fig_col[["ink2"]], fill = fig_col[["grid"]], linewidth = 0.5) +
    geom_point(colour = fig_col[["neuro"]], size = 1.4) +
    labs(title = "Galanin vs neuroendocrine proteins", subtitle = sprintf("%d proteins stored with galanin\n%s", nrow(nes), rho_txt(r_ne)),
         x = "neuroendocrine score (mean z)", y = "GAL (NPX)") + theme_pub()
} else ggplot() + labs(title = "Neuroendocrine proteins: fewer than 3 measurable") + theme_void()
p3C <- if (nrow(clin)) {
  cl <- clin |> mutate(analysis = factor(analysis, rev(unique(analysis))), is_score = analysis == "platelet score")
  ggplot(cl, aes(rho, analysis)) + geom_vline(xintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
    geom_linerange(aes(xmin = lo, xmax = hi, colour = is_score), linewidth = 0.45) + geom_point(aes(colour = is_score), size = 1.5) +
    scale_colour_manual(values = c(`FALSE` = fig_col[["ink"]], `TRUE` = fig_col[["platelet"]]), guide = "none") +
    facet_wrap(~parameter, ncol = 1, scales = "free_y") + labs(title = "Clinical links", subtitle = "adjusted: age, sex, Olink plate\n(for age: sex, Olink plate)",
                                                               x = "Spearman rho (95% CI)", y = NULL) + theme_pub() + row_guides()
} else ggplot() + labs(title = "CRP and age not available") + theme_void()
fig3 <- compose(list(p3A, compose(list(p3B, p3C), ncol = 1, rel_heights = c(0.9, 1.1), labels = c("B", "C"))), ncol = 2, rel_widths = c(1.1, 0.9), labels = c("A", ""))
save_figure(fig3, cfg, out, "Figure3", 180, 150)
src$Fig3A <- rob; if (!is.null(nscore)) src$Fig3B <- tibble(neuroendocrine_score = nscore, gal = g); src$Fig3C <- clin

# ---- Supplementary figures ---------------------------------------------------------------------------------------------------
s1a <- pset |> filter(status != "not on the panel") |> mutate(pct = 100 * frac_above_lod, gene = factor(gene, gene[order(pct)]),
                                                             used = status == "used")
pS1A <- ggplot(s1a, aes(pct, gene)) + geom_vline(xintercept = 100 * (cfg$min_detect_frac %||% 0.5), linetype = 2, colour = fig_col[["muted"]], linewidth = 0.3) +
  geom_point(aes(colour = used), size = 1.6) +
  scale_colour_manual(values = c(`TRUE` = fig_col[["platelet"]], `FALSE` = fig_col[["light"]]), labels = c(`TRUE` = "used", `FALSE` = "not measurable"), name = NULL) +
  labs(title = "Platelet proteins: detection",
       subtitle = if (any(pset$status == "not on the panel")) str_wrap(paste("not on the panel:", paste(pset$gene[pset$status == "not on the panel"], collapse = ", ")), 60) else NULL,
       x = "% of sera above LOD", y = NULL) + theme_pub() + row_guides() + theme(legend.position = "bottom")
pS1B <- ggplot(fit |> mutate(gene = factor(gene, gene[order(rho)])), aes(rho, gene)) + geom_vline(xintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = fig_col[["platelet"]], linewidth = 0.4) + geom_point(colour = fig_col[["platelet"]], size = 1.5) +
  labs(title = "Fit to the score", subtitle = "each protein vs the score of the others", x = "Spearman rho (95% CI)", y = NULL) + theme_pub() + row_guides()
save_figure(compose(list(pS1A, pS1B), ncol = 2), cfg, out, "FigureS1", 180, 110)
pS2 <- ggplot(axes |> mutate(component = factor(component, unique(component))), aes(rho, component, colour = what)) +
  geom_vline(xintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(width = 0.5), linewidth = 0.4) +
  geom_point(position = position_dodge(width = 0.5), size = 1.5) +
  scale_colour_manual(values = c(`galanin (GAL)` = fig_col[["ink"]], `platelet score` = fig_col[["platelet"]]), name = NULL) +
  scale_y_discrete(limits = rev) + labs(title = "Main axes of the serum proteome (principal components)", subtitle = sprintf("%s proteins; %% of variance in brackets", fmt_n(ncol(X))),
                                        x = "Spearman rho with the component (95% CI)", y = NULL) + theme_pub() + row_guides() + theme(legend.position = "bottom")
save_figure(pS2, cfg, out, "FigureS2", 110, 100)
if (nrow(bench)) {
  pS3 <- ggplot(bench |> mutate(name = factor(sprintf("%s vs Olink %s", lab, olink), sprintf("%s vs Olink %s", lab, olink)[order(rho)])), aes(rho, name)) +
    geom_vline(xintercept = 0, colour = fig_col[["ink2"]], linewidth = 0.3) +
    geom_linerange(aes(xmin = ci_low, xmax = ci_high), colour = fig_col[["ink"]], linewidth = 0.4) + geom_point(colour = fig_col[["ink"]], size = 1.5) +
    geom_text(aes(x = 1.02, label = sprintf("rank %d of %s", rank, fmt_n(of))), hjust = 0, size = 2, colour = fig_col[["ink2"]]) +
    coord_cartesian(xlim = c(min(0, min(bench$ci_low, na.rm = TRUE)), 1.35), clip = "off") +
    labs(title = "Olink vs clinical laboratory assays in the same sera", subtitle = "rank: place of the matching Olink protein among all proteins correlated with the lab value",
         x = "Spearman rho (95% CI)", y = NULL) + theme_pub() + row_guides()
  save_figure(pS3, cfg, out, "FigureS3", 150, 25 + 7 * nrow(bench))
}
src$FigS1 <- pset |> left_join(fit |> select(OlinkID, fit_rho = rho, fit_p = p), by = "OlinkID"); src$FigS2 <- axes; src$FigS3 <- bench

# ---- legends, key numbers, source data --------------------------------------------------------------------------------------
gs_txt <- if (nrow(gsea)) paste(sprintf("%s NES %.2f (p = %s)", pretty_set(gsea$pathway), gsea$NES, fmt_p(gsea$pval)), collapse = "; ") else "not available"
cl_txt <- if (nrow(clin)) paste(sprintf("%s: %s %.2f (p = %s)", clin$parameter, clin$analysis, clin$rho, fmt_p(clin$p)), collapse = "; ") else "not available"
legends <- c(
  "# Figure legends - galanin and platelets (draft; numbers filled in from the data)", "",
  sprintf("Data: %d LEIP sera (population controls), Olink Explore HT (%s); %s of %s proteins measurable (>= %.0f%% of sera above LOD). Spearman correlations; adjusted = partial Spearman for %s. The platelet set was defined from the literature before testing (config: platelets) - proteins made by megakaryocytes/platelets and released or shed when blood clots.",
          nrow(S), cfg$npx_column %||% "PCNormalizedNPX", fmt_n(ncol(d$M)), fmt_n(nrow(d$det)), 100 * (cfg$min_detect_frac %||% 0.5), covs_label(d$cov)), "",
  "**Figure 1. Serum galanin and its protein neighbourhood.**",
  sprintf("(A) Olink galanin (GAL) in each serum by Olink plate; dashed: median limit of detection (LOD) per plate. %d of %d sera above LOD; difference between plates p = %s (Kruskal-Wallis).", n_above, nrow(f1a), fmt_p(kw)),
  sprintf("(B) Spearman correlation of GAL with %s measurable proteins. Blue: platelet set (%d proteins); orange: proteins stored and released with galanin in neuroendocrine vesicles (%d); dashed: Bonferroni 5%%. %d proteins pass FDR < 5%%.",
          fmt_n(m_prot), nrow(plt), nrow(nes), sum(prot$fdr < 0.05, na.rm = TRUE)),
  sprintf("(C) Gene-set enrichment (GSEA) of pre-specified platelet gene sets along the ranking of proteins by correlation with GAL; ticks: members of each set. %s.", gs_txt), "",
  "**Figure 2. Galanin and a platelet release score.**",
  sprintf("(A) Spearman correlations between the platelet proteins and GAL (clustered); median pairwise rho among platelet proteins %.2f.", median(pair)),
  sprintf("(B) GAL vs the platelet release score (mean z-score of %d platelet proteins): %s; adjusted for %s: rho = %.2f, p = %s. Line: least squares with 95%% confidence band.",
          nrow(plt), rho_txt(r_all), covs_label(d$cov), r_adj$rho, fmt_p(r_adj$p)),
  "(C) GAL vs each platelet protein and the score (Spearman rho, 95% CI); circles: released granule cargo; triangles: membrane proteins shed on activation.",
  sprintf("(D) Specificity: correlation of GAL with the score of %s random sets of %d measurable proteins (grey) compared with the platelet set (blue line; one-sided p = %s)%s.",
          fmt_n(B), nrow(plt), fmt_p(p_rand), if (nrow(r_ne)) sprintf(" and the neuroendocrine set (orange; rho = %.2f, one-sided p = %s)", r_ne$rho, fmt_p(p_rand_ne)) else ""),
  sprintf("(E) Correlation of every measurable protein outside the platelet set with the platelet score; GAL ranks %d of %s (top %.1f%%).", gal_rank, fmt_n(n_allr), 100 * gal_rank / n_allr), "",
  "**Figure 3. Robustness of the galanin-platelet correlation, contrast with neuroendocrine proteins, and clinical links.**",
  sprintf("(A) GAL vs the platelet score across analyses: adjusted, within sex, without flagged sera, within each Olink plate, leaving out one serum or one platelet protein at a time (dotted: minimum to maximum), and with released or membrane proteins only. Leaving out one serum: rho %.2f to %.2f; one protein: %.2f to %.2f.",
          min(lr), max(lr), min(jk), max(jk)),
  sprintf("(B) GAL vs the score of %d neuroendocrine vesicle proteins stored with galanin: %s.", nrow(nes), if (nrow(r_ne)) rho_txt(r_ne) else "n/a"),
  sprintf("(C) GAL's associations with CRP and age: unadjusted, adjusted (age, sex, Olink plate; for age: sex, plate), additionally adjusted for the platelet score, and the platelet score itself. %s.", cl_txt), "",
  "**Figure S1.** Platelet proteins of the pre-defined set: (A) share of sera above LOD (dashed: threshold for use); (B) each protein vs the score of the other platelet proteins.",
  sprintf("**Figure S2.** GAL and the platelet score vs the first %d principal components of the %s measurable proteins.", kp, fmt_n(ncol(X))),
  "**Figure S3.** Agreement of Olink with clinical laboratory assays of the same proteins in the same sera (data quality); rank = place of the matching Olink protein among all proteins correlated with the lab value.", "",
  "## Notes for the text",
  "- The platelet hypothesis arose from the unbiased analysis of these data (Figure 1); Figures 2-3 test it with a pre-defined set, but in the same sera. Independent confirmation needs other data: galanin in matched plasma vs serum, and in the releasate of thrombin-activated washed platelets.",
  "- Serum was used: platelets release their contents during clotting, so a platelet link in serum can mean platelets carry galanin, or that clotting affects the Olink galanin signal.")
writeLines(legends, out_path(cfg, out, "figure_legends.md"))
save_csv(key, cfg, out, "key_results.csv")
writexl::write_xlsx(c(list(key_results = key, platelet_set = pset), src) |> keep(\(x) is.data.frame(x) && ncol(x) > 0) |>
                      map(\(x) mutate(x, across(where(is.list), \(v) map_chr(v, paste, collapse = ";")))),
                    out_path(cfg, out, "source_data.xlsx"))
for (i in seq_len(nrow(key))) msg("%s: %s", key$item[i], key$value[i])
msg("Figures: %s", file.path(cfg$paths$output, out))
