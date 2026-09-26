# 06 - ISF vs serum on matched MicroAD visits
# For every protein measured in both matrices:
#   within-subject correlation (repeated-measures correlation, rmcorr): do serum and ISF move
#     together over visits within a person?
#   between-subject correlation (Spearman of per-person means): do people with high ISF levels
#     also have high serum levels?
# Done separately for the lesional-site ISF (L: lesional/ex-lesional) and the non-lesional site
# (NL; healthy controls' skin is included here).
# Out: output/isf_serum/isf_serum_correlation.csv (+ plots of the top proteins)

source("R/utils.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
clear_outputs(cfg, "isf_serum")

shared <- intersect(rownames(wide$ISF), rownames(wide$Serum))
msg("%d proteins pass QC in both matrices", length(shared))

mk <- meta |> filter(cohort == "MicroAD") |> select(SampleID, SubjectID, visit, matrix, site, group)
serum <- mk |> filter(matrix == "Serum", SampleID %in% colnames(wide$Serum)) |> select(SubjectID, visit, serum_id = SampleID)
isf <- mk |> filter(matrix == "ISF", SampleID %in% colnames(wide$ISF)) |>
  mutate(site2 = if_else(site == "L", "L", "NL")) |>           # healthy skin (H) counts as NL
  select(SubjectID, visit, site2, group, isf_id = SampleID)
pairs <- inner_join(isf, serum, by = c("SubjectID", "visit"))
save_csv(count(pairs, site2, group), cfg, "isf_serum", "matched_pairs.csv")

within_r <- function(x, y, subj) {
  # repeated-measures correlation (Bakdash & Marusich 2017): ANCOVA with subject as factor
  keep <- subj %in% names(which(table(subj) >= 2)) & !is.na(x) & !is.na(y)
  if (sum(keep) < 5 || length(unique(subj[keep])) < 2) return(c(r = NA, p = NA, df = NA))
  r <- rmcorr::rmcorr(participant = subj, measure1 = x, measure2 = y,
                      dataset = data.frame(subj = subj[keep], x = x[keep], y = y[keep]) |>
                        setNames(c("subj", "x", "y")) |> transform(subj = factor(subj)))
  c(r = r$r, p = r$p, df = r$df)
}

res <- map(split(pairs, pairs$site2), \(pp) {
  map(shared, \(a) {
    x <- wide$ISF[a, pp$isf_id]; y <- wide$Serum[a, pp$serum_id]
    w <- suppressWarnings(within_r(x, y, pp$SubjectID))
    means <- tibble(s = pp$SubjectID, x, y) |> group_by(s) |> summarise(x = mean(x), y = mean(y))
    b <- suppressWarnings(cor.test(means$x, means$y, method = "spearman", exact = FALSE))
    tibble(OlinkID = a, site = pp$site2[1], n_pairs = nrow(pp), n_subjects = n_distinct(pp$SubjectID),
           r_within = w[["r"]], p_within = w[["p"]], r_between = unname(b$estimate), p_between = b$p.value,
           isf_minus_serum = median(x - y, na.rm = TRUE))
  }) |> bind_rows()
}) |> bind_rows() |>
  group_by(site) |>
  mutate(fdr_within = p.adjust(p_within, "BH"), fdr_between = p.adjust(p_between, "BH")) |>
  ungroup() |>
  mutate(significant = coalesce(fdr_within < cfg$stats$fdr, FALSE) | coalesce(fdr_between < cfg$stats$fdr, FALSE)) |>
  left_join(assay_map, by = "OlinkID") |> relocate(Assay, .after = OlinkID)

save_csv(res, cfg, "isf_serum", "isf_serum_correlation.csv")
sig <- res |> filter(significant) |> distinct(OlinkID, Assay)
save_csv(sig, cfg, "isf_serum", "significant_proteins.csv")
msg("%d proteins significantly correlated between ISF and serum (either site, within or between subjects)", nrow(sig))
print(res |> group_by(site) |> summarise(within = sum(fdr_within < cfg$stats$fdr, na.rm = TRUE),
                                         between = sum(fdr_between < cfg$stats$fdr, na.rm = TRUE)))

top <- res |> filter(significant) |> slice_min(pmin(fdr_within, fdr_between, na.rm = TRUE), n = 9, with_ties = FALSE)
if (nrow(top)) {
  pd <- map(seq_len(nrow(top)), \(i) {
    pp <- pairs |> filter(site2 == top$site[i])
    tibble(Assay = sprintf("%s (%s site)", top$Assay[i], top$site[i]), SubjectID = pp$SubjectID, group = pp$group,
           isf = wide$ISF[top$OlinkID[i], pp$isf_id], serum = wide$Serum[top$OlinkID[i], pp$serum_id])
  }) |> bind_rows()
  p <- ggplot(pd, aes(isf, serum, colour = SubjectID)) + geom_point(size = 1.2) +
    geom_line(stat = "smooth", method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, alpha = 0.6) +
    facet_wrap(~Assay, scales = "free") + guides(colour = "none") +
    labs(title = "Top ISF-serum correlations (lines: within-subject fits)", x = "ISF NPX", y = "serum NPX")
  save_plot(p, cfg, "isf_serum", "top_correlations.png", width = 10, height = 8)
}
