# Relapse and lesional signatures per visit (part of run_all.R; shortcut: run_visits_v1_v4.R), MicroAD dISF.
#   1 relapse vs non-relapse at the lesion site, at each visit (V1-V4): limma, exact relabelling check, volcano; proteins
#     that differ in the same direction at several visits (the relapse signature over time)
#   2 GSEA of each of these comparisons, and the gene sets over the visits
#   3 lesion site vs non-lesional skin of the same patients (paired), at each visit: limma, volcano
#   4 AKR1B10 (config track_proteins): lesion site vs non-lesional over the visits
#   5 mTORC1 (config signature_sets): GSEA lesion site vs non-lesional at each visit, and its proteins at each visit
# In:  output of steps 01-02 (data/npx_clean.rds); MSigDB / KEGG gene sets
# Out: output/visits_V1_V4/ - 1_ 2_ 3_ 5_ Excel files; figures as PNG (300 dpi, 16:9, large fonts: paste into PowerPoint as they
#      are) and PDF; Visits_V1-V4_summary.pdf

source("R/utils.R")
source("R/qc_overview.R")
source("R/visit3.R")
source("R/relapse_v1.R")
source("R/relapse_followup.R")
source("R/visits_v1_v4.R")
cfg <- load_config()
vv <- cfg$visits_v1_v4 %||% list()
visits <- as.integer(unlist(vv$visits %||% 1:4)); inc_rel <- isTRUE(vv$include_relapse_visit)
n_lab <- vv$volcano_labels %||% 20; n_gs <- vv$gsea_top %||% 20
track <- unlist(vv$track_proteins %||% "AKR1B10"); sig_sets <- unlist(vv$signature_sets %||% "HALLMARK_MTORC1_SIGNALING")
fdr <- cfg$stats$fdr %||% 0.05; min_f <- cfg$qc$min_detect_frac %||% 0.5; min_n <- cfg$stats$min_group_n %||% 3
gmin <- cfg$enrichment$min_size %||% 10; gmax <- cfg$enrichment$max_size %||% 500
colours <- cfg$qc_overview$colours %||% list(above = "#69005F", below = "#FF506E")
out <- "visits_V1_V4"
clear_outputs(cfg, out)
set.seed(1)

clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
info <- clean |> distinct(OlinkID, Assay, UniProt)
ad <- clean |> filter(cohort == "MicroAD", group == "AD", matrix == "ISF")
vn <- set_names(visits, paste0("V", visits))
fmt_p <- \(x) ifelse(is.na(x), "-", formatC(x, digits = 2, format = "g"))
sig3 <- \(x) x |> mutate(across(where(is.double), \(v) signif(v, 3)))
top_prot <- \(r, dir, k = 5) { x <- r |> filter(p < 0.05, sign(logFC) == dir) |> arrange(p); if (nrow(x)) paste(head(x$Assay, k), collapse = ", ") else "none" }
manifest_txt <- sprintf("sample manifest: %s", basename(cfg$paths$manifest))
sets <- v3_gene_sets(cfg)

# ---- samples ------------------------------------------------------------------------------------------------------------------
# relapse vs non-relapse at the lesion site: relapsers' samples at their relapse visit (lesional again) are left out unless
# include_relapse_visit is true - the comparison is about patients in remission who relapse later
rel_samples <- function(v) ad |> filter(site == "L", visit_num == v, relapse %in% c("relapse", "non-relapse")) |>
  distinct(SampleID, SubjectID, relapse, relapse_visit, state) |>
  mutate(visit = paste0("V", v), at_relapse_visit = relapse == "relapse" & coalesce(relapse_visit <= v, FALSE),
         included = inc_rel | !at_relapse_visit, grp = factor(relapse, c("non-relapse", "relapse")))
# lesion site vs non-lesional skin: all AD patients with both samples at the visit (paired)
ln_samples <- function(v) {
  x <- ad |> filter(site %in% c("L", "NL"), visit_num == v) |> distinct(SampleID, SubjectID, site, state, relapse)
  x |> filter(SubjectID %in% (x |> count(SubjectID) |> filter(n == 2) |> pull(SubjectID))) |>
    mutate(visit = paste0("V", v), grp = factor(if_else(site == "L", "lesion site", "non-lesional"), c("non-lesional", "lesion site")))
}

# =============================================================================================================================
# 1 / 2  Relapse vs non-relapse at the lesion site, per visit: limma, exact relabelling check, GSEA
# =============================================================================================================================
rel <- map(vn, \(v) {
  s <- rel_samples(v); si <- s |> filter(included)
  n <- c(relapse = sum(si$relapse == "relapse"), `non-relapse` = sum(si$relapse == "non-relapse"))
  base <- list(v = v, samples = s, n = n, n_left_out = sum(!s$included))
  left_txt <- if (base$n_left_out) sprintf("; %d relapser(s) at the relapse visit left out", base$n_left_out) else ""
  if (min(n) < min_n) return(c(base, list(skipped = sprintf("not tested: fewer than %d patients in a group (relapse %d, non-relapse %d%s)",
                                                              min_n, n[["relapse"]], n[["non-relapse"]], left_txt))))
  r <- de_limma(clean, si, min_f = min_f, fdr = fdr)
  M <- patient_matrix(clean, si, r$OlinkID); M <- M[r$OlinkID, sort(colnames(M)), drop = FALSE]
  y <- si$grp[match(colnames(M), si$SubjectID)] == "relapse"
  rc <- relabel_check(M, y)                                    # every way to choose the relapse group
  null <- colSums(rc$up_hits) + colSums(rc$down_hits); obs <- sum(r$p < 0.05, na.rm = TRUE)
  r <- r |> left_join(rc$per_protein |> select(OlinkID, p_exact = p_exact_two_sided), by = "OlinkID") |> relocate(p_exact, .after = FDR)
  msg("relapse vs non-relapse %s: %d vs %d, %d proteins, %d with p < 0.05 (random labels: median %.0f), %d with FDR < %.2f",
      paste0("V", v), n[["relapse"]], n[["non-relapse"]], nrow(r), obs, median(null), sum(r$FDR < fdr, na.rm = TRUE), fdr)
  c(base, list(left_txt = left_txt, res = r, gsea = gsea_keyed(r, sets, gmin, gmax),
               chance = list(obs = obs, med = median(null), q90 = unname(quantile(null, 0.9)), P = mean(null >= obs), n_lab = rc$n_labellings)))
})
rel_ok <- keep(rel, \(x) is.null(x$skipped))

# the relapse signature over time: proteins with p < 0.05 in the same direction at several tested visits
sig_long <- imap(rel_ok, \(x, nm) x$res |> transmute(OlinkID, Assay, visit = nm, log2FC = logFC, p, FDR)) |> bind_rows()
sig_tab <- if (nrow(sig_long)) sig_long |> group_by(OlinkID, Assay) |>
  summarise(visits_higher = sum(p < 0.05 & log2FC > 0, na.rm = TRUE), visits_lower = sum(p < 0.05 & log2FC < 0, na.rm = TRUE),
            visits_tested = sum(!is.na(p)), mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") |>
  mutate(direction = case_when(visits_higher >= 2 & visits_lower == 0 ~ "higher in relapsers at >= 2 visits",
                               visits_lower >= 2 & visits_higher == 0 ~ "lower in relapsers at >= 2 visits", TRUE ~ NA_character_)) |>
  left_join(sig_long |> pivot_wider(id_cols = OlinkID, names_from = visit, values_from = c(log2FC, p)), by = "OlinkID") |>
  arrange(is.na(direction), desc(pmax(visits_higher, visits_lower)), desc(abs(mean_log2FC))) else tibble()
sig_hits <- if (nrow(sig_tab)) sig_tab |> filter(!is.na(direction)) else tibble()

# gene sets over the visits: FDR < fdr at one visit or more
gs_long <- imap(rel_ok, \(x, nm) if (nrow(x$gsea)) x$gsea |> transmute(pathway, collection, gene_set, visit = nm, NES, p, FDR) else NULL) |> bind_rows()
gs_across <- if (nrow(gs_long)) gs_long |> group_by(pathway) |> filter(any(FDR < fdr)) |> ungroup() |>
  pivot_wider(id_cols = c(pathway, collection, gene_set), names_from = visit, values_from = c(NES, FDR)) |> select(-pathway) else tibble()

# =============================================================================================================================
# 3  Lesion site vs non-lesional skin (paired), per visit; GSEA of each (for section 5)
# =============================================================================================================================
ln <- map(vn, \(v) {
  s <- ln_samples(v); np <- n_distinct(s$SubjectID)
  states <- s |> filter(site == "L") |> count(state) |> mutate(t = paste(n, state)) |> pull(t) |> paste(collapse = ", ")
  base <- list(v = v, samples = s, n_pairs = np, states = states)
  if (np < min_n) return(c(base, list(skipped = sprintf("not tested: fewer than %d patients with both samples (%d)", min_n, np))))
  r <- de_limma(clean, s, paired = TRUE, min_f = min_f, fdr = fdr)
  msg("lesion site vs non-lesional %s: %d pairs (%s), %d proteins, %d with FDR < %.2f", paste0("V", v), np, states, nrow(r), sum(r$FDR < fdr, na.rm = TRUE), fdr)
  c(base, list(res = r, gsea = gsea_keyed(r, sets, gmin, gmax)))
})
ln_ok <- keep(ln, \(x) is.null(x$skipped))

# =============================================================================================================================
# 4  Tracked proteins (AKR1B10): lesion site vs non-lesional over the visits
# =============================================================================================================================
tp <- assay_ids(info, track)
track_missing <- setdiff(track, tp$gene)
pairs_all <- map(ln, "samples") |> bind_rows()
tp_vals <- clean |> filter(SampleID %in% pairs_all$SampleID, OlinkID %in% tp$OlinkID) |>
  transmute(Assay, OlinkID, SubjectID, visit_num, visit = paste0("V", visit_num), site = if_else(site == "L", "lesion site", "non-lesional"), state,
            NPX = value, below_LOD = below_lod)
tp_means <- tp_vals |> group_by(Assay, OlinkID, visit_num, visit, site) |>
  summarise(n = sum(!is.na(NPX)), mean_NPX = mean(NPX, na.rm = TRUE), se = sd(NPX, na.rm = TRUE) / sqrt(n), pct_above_LOD = 100 * mean(!below_LOD, na.rm = TRUE), .groups = "drop")
tp_hc <- clean |> filter(cohort == "MicroAD", group == "HC", matrix == "ISF", OlinkID %in% tp$OlinkID) |> group_by(Assay) |>
  summarise(mean_NPX = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))), .groups = "drop")
# paired test per visit: the limma result of section 3; a protein below LOD there (not tested) gets a paired t-test
tp_stats <- if (!nrow(tp)) tibble(Assay = character(), visit = character(), log2FC = numeric(), p = numeric(), FDR = numeric(), test = character()) else
  map(names(ln), \(nm) map(seq_len(nrow(tp)), \(i) {
  x <- ln[[nm]]; o <- tp$OlinkID[i]
  r <- if (!is.null(x$res)) x$res |> filter(OlinkID == o) else tibble()
  if (nrow(r)) return(tibble(Assay = tp$Assay[i], visit = nm, log2FC = r$logFC, p = r$p, FDR = r$FDR, test = "limma, paired (as in the DEG table)"))
  w <- tp_vals |> filter(OlinkID == o, visit == nm) |> select(SubjectID, site, NPX) |> pivot_wider(names_from = site, values_from = NPX) |> drop_na()
  if (nrow(w) < 3) return(tibble(Assay = tp$Assay[i], visit = nm, log2FC = NA_real_, p = NA_real_, FDR = NA_real_, test = "too few pairs"))
  tibble(Assay = tp$Assay[i], visit = nm, log2FC = mean(w$`lesion site` - w$`non-lesional`), p = t.test(w$`lesion site`, w$`non-lesional`, paired = TRUE)$p.value,
         FDR = NA_real_, test = "paired t-test (mostly below LOD, so not in the DEG table)")
}) |> bind_rows()) |> bind_rows()
tp_table <- tp_means |> select(Assay, visit, site, n, mean_NPX, se, pct_above_LOD) |>
  pivot_wider(names_from = site, values_from = c(n, mean_NPX, se, pct_above_LOD), names_glue = "{.value}_{site}") |>
  left_join(tp_stats, by = c("Assay", "visit")) |> left_join(tp_hc |> rename(mean_NPX_healthy = mean_NPX, se_healthy = se), by = "Assay") |>
  rename_with(\(x) str_replace_all(x, " ", "_")) |> arrange(Assay, visit)

# =============================================================================================================================
# 5  Signature gene sets (mTORC1): GSEA lesion site vs non-lesional at each visit, and their proteins at each visit
# =============================================================================================================================
# readable gene set name for figures and tables, e.g. HALLMARK_MTORC1_SIGNALING -> "Hallmark mTORC1 signaling"
case_fix <- c(mtorc1 = "mTORC1", mtor = "mTOR", nfkb = "NF-kB", tnfa = "TNF-alpha", pi3k = "PI3K", akt = "AKT", jak = "JAK", stat3 = "STAT3",
              stat5 = "STAT5", il2 = "IL2", il6 = "IL6", kras = "KRAS", myc = "MYC", e2f = "E2F", g2m = "G2M", p53 = "p53", tgf = "TGF", wnt = "WNT",
              uv = "UV", dna = "DNA", t = "T", b = "B", gobp = "GO:BP", kegg = "KEGG")
set_name <- \(k) map_chr(k, \(x) {
  x <- str_to_lower(str_replace_all(x, "_", " "))
  for (w in names(case_fix)) x <- str_replace_all(x, paste0("\\b", w, "\\b"), case_fix[[w]])
  paste0(toupper(substr(x, 1, 1)), substring(x, 2))
})
sig_sets_found <- intersect(sig_sets, names(sets))
sig_nes <- map(sig_sets_found, \(k) imap(ln_ok, \(x, nm) {
  g <- if (nrow(x$gsea)) x$gsea |> filter(pathway == k) else tibble()
  tibble(set = set_name(k), msigdb_name = k, visit = nm, visit_num = x$v, NES = if (nrow(g)) g$NES else NA_real_, p = if (nrow(g)) g$p else NA_real_, FDR = if (nrow(g)) g$FDR else NA_real_,
         proteins = if (nrow(g)) g$measured_proteins_in_set else NA_integer_, leading_edge = if (nrow(g)) g$leading_edge else "not tested (too few measured proteins)")
}) |> bind_rows()) |> bind_rows()
sig_prot <- map(sig_sets_found, \(k) imap(ln_ok, \(x, nm) x$res |> filter(map_lgl(str_split(Assay, "_"), \(g) any(g %in% sets[[k]]))) |>
                                          transmute(set = set_name(k), Assay, OlinkID, visit = nm, log2FC = logFC, p, FDR)) |> bind_rows()) |> bind_rows()

# =============================================================================================================================
# Excel files
# =============================================================================================================================
mark_lab <- \(r) r |> mutate(labelled_in_volcano = OlinkID %in% (r |> filter(!is.na(p)) |> slice_min(p, n = n_lab, with_ties = FALSE) |> pull(OlinkID)))
skip_sheet <- \(x) tibble(note = x$skipped)
rel_summary <- imap(rel, \(x, nm) tibble(visit = nm, relapse = x$n[["relapse"]], `non-relapse` = x$n[["non-relapse"]], left_out_at_relapse_visit = x$n_left_out,
  proteins_tested = if (is.null(x$res)) NA_integer_ else nrow(x$res),
  p_below_0.05 = if (is.null(x$res)) NA_integer_ else sum(x$res$p < 0.05, na.rm = TRUE),
  higher_in_relapse = if (is.null(x$res)) NA_integer_ else sum(x$res$p < 0.05 & x$res$logFC > 0, na.rm = TRUE),
  lower_in_relapse = if (is.null(x$res)) NA_integer_ else sum(x$res$p < 0.05 & x$res$logFC < 0, na.rm = TRUE),
  random_labels_median = if (is.null(x$chance)) NA_real_ else x$chance$med, random_labels_P = if (is.null(x$chance)) NA_real_ else x$chance$P,
  FDR_below = if (is.null(x$res)) NA_integer_ else sum(x$res$FDR < fdr, na.rm = TRUE),
  gene_sets_FDR_higher = if (is.null(x$gsea) || !nrow(x$gsea)) NA_integer_ else sum(x$gsea$FDR < fdr & x$gsea$NES > 0),
  gene_sets_FDR_lower = if (is.null(x$gsea) || !nrow(x$gsea)) NA_integer_ else sum(x$gsea$FDR < fdr & x$gsea$NES < 0),
  top_higher = if (is.null(x$res)) x$skipped else top_prot(x$res, 1), top_lower = if (is.null(x$res)) "" else top_prot(x$res, -1))) |> bind_rows()
ln_summary <- imap(ln, \(x, nm) tibble(visit = nm, pairs = x$n_pairs, lesion_site_state = x$states,
  proteins_tested = if (is.null(x$res)) NA_integer_ else nrow(x$res),
  FDR_higher_at_lesion_site = if (is.null(x$res)) NA_integer_ else sum(x$res$FDR < fdr & x$res$logFC > 0, na.rm = TRUE),
  FDR_lower_at_lesion_site = if (is.null(x$res)) NA_integer_ else sum(x$res$FDR < fdr & x$res$logFC < 0, na.rm = TRUE),
  top_higher = if (is.null(x$res)) x$skipped else top_prot(x$res, 1), top_lower = if (is.null(x$res)) "" else top_prot(x$res, -1))) |> bind_rows()
sample_sheet <- \(lst, extra) map(lst, "samples") |> bind_rows() |> select(any_of(c("visit", "SampleID", "SubjectID", "site", "relapse", "state", extra)))

writexl::write_xlsx(c(list(
  README = tibble(item = c("question", "samples", "test", "chance check", "signature", "columns", "data"),
    note = c("Relapse vs non-relapse in dISF at the lesion site, at each visit: which proteins differ already before the relapse (relapse signature)?",
             sprintf("MicroAD AD patients with a relapse / non-relapse label (dropouts excluded), one lesion-site sample each. %s",
                     if (inc_rel) "Relapsers' samples at their relapse visit are included." else "Relapsers' samples at their relapse visit (lesional again) are left out (config include_relapse_visit), so all are in remission from V2 on."),
             sprintf("limma per protein, log2FC = relapse - non-relapse; proteins above LOD in >= %.0f%% of one group; FDR (Benjamini-Hochberg) < %.2f; a visit is tested only with >= %d patients per group", 100 * min_f, fdr, min_n),
             "p_exact: the same test repeated for every way to choose the relapse group (exact p of |t|); summary: proteins with p < 0.05 with the real labels vs the median of all random labellings (random_labels_P = share of labellings with at least as many)",
             "sheet signature_over_visits: proteins with p < 0.05 in the same direction at >= 2 tested visits (the relapse signature that persists), with log2FC and p per visit",
             sprintf("labelled_in_volcano: the %d proteins with the lowest p (labelled in the volcano plots)", n_lab), manifest_txt)),
  summary = rel_summary |> sig3()),
  map(rel, \(x) if (is.null(x$res)) skip_sheet(x) else mark_lab(x$res) |> sig3()),
  list(signature_over_visits = if (nrow(sig_tab)) sig_tab |> sig3() else tibble(note = "no visit tested"),
       samples = sample_sheet(rel, c("relapse_visit", "at_relapse_visit", "included")))),
  out_path(cfg, out, "1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx"))
writexl::write_xlsx(c(list(
  README = tibble(item = c("question", "method", "gene_sets_over_visits", "data"),
    note = c("Which gene sets differ between relapsers and non-relapsers at the lesion site, at each visit?",
             sprintf("GSEA (fgsea) on the limma t of 1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx (all tested proteins ranked); NES > 0 = higher in relapsers. MSigDB %s and KEGG; gene sets with %d-%d measured proteins",
                     paste(cfg$enrichment$collections %||% c("H", "C2:CP:REACTOME", "C5:GO:BP"), collapse = ", "), gmin, gmax),
             sprintf("gene sets with FDR < %.2f at one visit or more: NES and FDR at every visit", fdr), manifest_txt)),
  gene_sets_over_visits = gs_across |> sig3()),
  map(rel, \(x) if (is.null(x$res)) skip_sheet(x) else if (nrow(x$gsea)) x$gsea |> select(-pathway) else tibble(note = "no gene set tested"))),
  out_path(cfg, out, "2_relapse_vs_nonrelapse_GSEA_V1-V4.xlsx"))
tp_sheet <- if (length(track) == 1) track else "tracked_proteins"
writexl::write_xlsx(c(list(
  README = tibble(item = c("question", "samples", "test", tp_sheet, "data"),
    note = c("Lesion site vs non-lesional skin of the same patients, at each visit",
             "MicroAD AD patients with both a lesion-site and a non-lesional dISF sample at the visit (all relapse labels). The lesion site is lesional at V1 and cleared later, except at a relapse visit (column lesion_site_state of the summary)",
             sprintf("limma, paired (patient as factor), log2FC = lesion site - non-lesional; proteins above LOD in >= %.0f%% of one group; FDR < %.2f", 100 * min_f, fdr),
             sprintf("%s at each visit: mean NPX +/- SE at the lesion site and in non-lesional skin, the paired test (from the DEG table), healthy skin for reference%s",
                     paste(tp$Assay, collapse = ", "), if (length(track_missing)) sprintf("; not in the data: %s", paste(track_missing, collapse = ", ")) else ""), manifest_txt)),
  summary = ln_summary |> sig3()),
  map(ln, \(x) if (is.null(x$res)) skip_sheet(x) else mark_lab(x$res) |> sig3()),
  set_names(list(tp_table |> sig3(), tp_vals |> select(-OlinkID) |> sig3()), c(tp_sheet, paste0(tp_sheet, "_per_patient")) |> str_trunc(31, ellipsis = "")),
  list(samples = sample_sheet(ln, character()))),
  out_path(cfg, out, "3_lesional_vs_nonlesional_DEG_V1-V4.xlsx"))
writexl::write_xlsx(c(list(
  README = tibble(item = c("question", "method", "signature", "signature_proteins", "data"),
    note = c(sprintf("Is %s higher at the lesion site than in non-lesional skin at each visit?", paste(map_chr(sig_sets, set_name), collapse = ", ")),
             "GSEA (fgsea) on the paired limma t of 3_lesional_vs_nonlesional_DEG_V1-V4.xlsx; NES > 0 = higher at the lesion site",
             "NES, p, FDR and leading edge of the signature gene set(s) at each visit (config visits_v1_v4$signature_sets)",
             "every measured protein of the signature gene set(s): log2FC, p, FDR at each visit", manifest_txt)),
  signature = if (nrow(sig_nes)) sig_nes |> sig3() else tibble(note = sprintf("gene set(s) not found: %s", paste(sig_sets, collapse = ", "))),
  signature_proteins = if (nrow(sig_prot)) sig_prot |> pivot_wider(id_cols = c(set, Assay, OlinkID), names_from = visit, values_from = c(log2FC, p, FDR)) |> sig3() else tibble()),
  map(ln, \(x) if (is.null(x$res)) skip_sheet(x) else if (nrow(x$gsea)) x$gsea |> select(-pathway) else tibble(note = "no gene set tested"))),
  out_path(cfg, out, "5_lesional_vs_nonlesional_GSEA_V1-V4.xlsx"))

# =============================================================================================================================
# Figures: PNG (300 dpi, 16:9, large fonts - safe to paste into PowerPoint) and PDF, collected in order for the PDF summary
# =============================================================================================================================
fi <- qc_font_setup(cfg$qc_overview$font %||% "Nimbus Sans")
figs <- list()
fig <- function(p, file, title, note = NULL, width = 12.5, height = 6.4) {
  if (is.null(p)) return(invisible(NULL))
  if (!is.null(note)) note <- str_wrap(note, 115)     # long subtitles wrap instead of running off the slide
  pt <- p + labs(title = title, subtitle = note)
  for (type in c("pdf", "png")) {
    open_device(out_path(cfg, out, paste0(file, ".", type)), type, fi, width, height)
    print(set_family(pt, fi$family))
    close_device(fi)
  }
  figs[[length(figs) + 1]] <<- list(p = p, title = title, note = note)
  invisible(p)
}
for (nm in names(rel)) {
  x <- rel[[nm]]
  if (!is.null(x$skipped)) next                      # noted in the summary and the Excel file
  fig(volcano_slide(x$res, colours, n_lab, fdr, up = "higher in relapse", down = "lower in relapse", xlab = "log2 fold change (relapse - non-relapse)"),
      paste0("1_relapse_volcano_", nm), sprintf("Relapse vs non-relapse, dISF lesion site %s (%d vs %d patients)", nm, x$n[["relapse"]], x$n[["non-relapse"]]),
      sprintf("%d proteins; %d with p < 0.05 (random labels: median %.0f), %d with FDR < %.2f (black ring); labelled: the %d lowest p%s",
              nrow(x$res), x$chance$obs, x$chance$med, sum(x$res$FDR < fdr, na.rm = TRUE), fdr, n_lab, x$left_txt))
}
if (nrow(sig_hits)) {
  show <- sig_hits |> slice_head(n = 40)
  fig(heat_slide(sig_long |> filter(OlinkID %in% show$OlinkID) |> transmute(row = factor(Assay, rev(unique(show$Assay))), col = visit, value = log2FC, star = p < 0.05),
                 colours, "log2FC\nrelapse -\nnon-relapse"),
      "1b_relapse_signature_over_visits", "Relapse signature over the visits: proteins differing in the same direction at >= 2 visits",
      sprintf("%d proteins (top %d shown); * = p < 0.05 at that visit; full list: sheet signature_over_visits", nrow(sig_hits), nrow(show)), height = 6.4)
}
for (nm in names(rel_ok)) {
  x <- rel_ok[[nm]]
  fig(gsea_slide(x$gsea, colours, fdr, n_gs, up = "higher in relapse", down = "lower in relapse"),
      paste0("2_relapse_GSEA_", nm), sprintf("Gene sets, relapse vs non-relapse, dISF lesion site %s", nm),
      sprintf("GSEA on the limma t; %d gene sets with FDR < %.2f (%d higher, %d lower in relapsers); shown: the %d lowest p",
              sum(x$gsea$FDR < fdr), fdr, sum(x$gsea$FDR < fdr & x$gsea$NES > 0), sum(x$gsea$FDR < fdr & x$gsea$NES < 0), min(n_gs, nrow(x$gsea))))
}
if (nrow(gs_long)) {
  top_gs <- gs_long |> group_by(pathway) |> filter(any(FDR < fdr)) |> summarise(best = min(FDR), .groups = "drop") |> slice_min(best, n = 30, with_ties = FALSE)
  if (nrow(top_gs)) {
    d <- gs_long |> filter(pathway %in% top_gs$pathway) |>
      mutate(label = str_trunc(paste0(if_else(gene_set == toupper(gene_set), str_to_sentence(gene_set), gene_set), " (", collection, ")"), 65))
    ord <- d |> group_by(label) |> summarise(m = mean(NES), .groups = "drop") |> arrange(m) |> pull(label)
    fig(heat_slide(d |> transmute(row = factor(label, ord), col = visit, value = NES, star = FDR < fdr), colours, "NES\n(> 0 higher\nin relapse)", base = 12),
        "2b_relapse_GSEA_over_visits", "Gene sets over the visits, relapse vs non-relapse (lesion site)",
        sprintf("gene sets with FDR < %.2f at one visit or more (the 30 with the lowest FDR); * = FDR < %.2f at that visit", fdr, fdr))
  }
}
for (nm in names(ln)) {
  x <- ln[[nm]]
  if (!is.null(x$skipped)) next
  fig(volcano_slide(x$res, colours, n_lab, fdr, up = "higher at lesion site", down = "lower at lesion site", xlab = "log2 fold change (lesion site - non-lesional)"),
      paste0("3_lesional_volcano_", nm), sprintf("Lesion site vs non-lesional skin, dISF %s (%d patients, paired)", nm, x$n_pairs),
      sprintf("lesion site: %s. %d proteins; %d with FDR < %.2f (black ring); labelled: the %d lowest p", x$states, nrow(x$res), sum(x$res$FDR < fdr, na.rm = TRUE), fdr, n_lab))
}
if (nrow(tp_means)) {
  pl <- tp_stats |> mutate(visit_num = as.integer(str_remove(visit, "V")),
         label = if_else(is.na(p), "", sprintf("p = %s%s", fmt_p(p), if_else(coalesce(FDR < fdr, FALSE), "*", ""))))
  fig(levels_slide(tp_means, tp_hc, pl), "4_AKR1B10_lesional_vs_nonlesional_over_visits",
      sprintf("%s: lesion site vs non-lesional skin over the visits", paste(unique(tp_means$Assay), collapse = ", ")),
      sprintf("mean +/- SE of the patients with both samples at the visit; grey band: healthy skin (mean +/- SE); p: paired limma (* FDR < %.2f)%s", fdr,
              if (length(track_missing)) sprintf("; not in the data: %s", paste(track_missing, collapse = ", ")) else ""))
}
if (nrow(sig_nes) && any(!is.na(sig_nes$NES))) {
  fig(nes_slide(sig_nes |> filter(!is.na(NES)), colours), "5_mTORC1_lesional_vs_nonlesional_by_visit",
      sprintf("%s: lesion site vs non-lesional skin at each visit", paste(unique(sig_nes$set), collapse = ", ")),
      "GSEA on the paired limma t; NES > 0 = higher at the lesion site")
  prot_show <- sig_prot |> group_by(set, Assay) |> summarise(m = mean(abs(log2FC), na.rm = TRUE), any_sig = any(p < 0.05), .groups = "drop") |>
    arrange(desc(any_sig), desc(m)) |> group_by(set) |> slice_head(n = 35) |> ungroup()
  d <- sig_prot |> semi_join(prot_show, by = c("set", "Assay"))
  ord <- d |> group_by(Assay) |> summarise(m = mean(log2FC, na.rm = TRUE), .groups = "drop") |> arrange(m) |> pull(Assay)
  fig(heat_slide(d |> transmute(set, row = factor(Assay, unique(ord)), col = visit, value = log2FC, star = p < 0.05), colours, "log2FC\nlesion site -\nnon-lesional", base = 12) +
        facet_wrap(~set, scales = "free_y"),
      "5b_mTORC1_proteins_by_visit", sprintf("%s proteins: lesion site vs non-lesional at each visit", paste(unique(sig_nes$set), collapse = ", ")),
      sprintf("measured proteins of the gene set (up to 35: those with p < 0.05 at a visit first); * = p < 0.05 (paired limma); all in 5_lesional_vs_nonlesional_GSEA_V1-V4.xlsx"))
}

# =============================================================================================================================
# Summary text and PDF summary
# =============================================================================================================================
rel_lines <- imap_chr(rel, \(x, nm) if (!is.null(x$skipped)) sprintf("%s: %s", nm, x$skipped) else
  sprintf("%s (%d vs %d%s): %d proteins with p < 0.05 (random labels: median %.0f, P = %s), %d with FDR < %.2f; gene sets FDR < %.2f: %d higher, %d lower in relapsers. Higher: %s. Lower: %s.",
          nm, x$n[["relapse"]], x$n[["non-relapse"]], x$left_txt, x$chance$obs, x$chance$med, fmt_p(x$chance$P), sum(x$res$FDR < fdr, na.rm = TRUE), fdr, fdr,
          sum(x$gsea$FDR < fdr & x$gsea$NES > 0), sum(x$gsea$FDR < fdr & x$gsea$NES < 0), top_prot(x$res, 1), top_prot(x$res, -1)))
sig_line <- if (nrow(sig_hits)) sprintf("Relapse signature over the visits (p < 0.05, same direction, >= 2 visits): %d higher (%s), %d lower (%s).",
                                        sum(sig_hits$visits_higher >= 2), paste(head(sig_hits$Assay[sig_hits$visits_higher >= 2], 8), collapse = ", "),
                                        sum(sig_hits$visits_lower >= 2), paste(head(sig_hits$Assay[sig_hits$visits_lower >= 2], 8), collapse = ", ")) else
  "No protein differs in the same direction (p < 0.05) at two or more tested visits."
ln_lines <- imap_chr(ln, \(x, nm) if (!is.null(x$skipped)) sprintf("%s: %s", nm, x$skipped) else
  sprintf("%s (%d pairs; lesion site %s): %d proteins with FDR < %.2f (%d higher, %d lower at the lesion site).", nm, x$n_pairs, x$states,
          sum(x$res$FDR < fdr, na.rm = TRUE), fdr, sum(x$res$FDR < fdr & x$res$logFC > 0, na.rm = TRUE), sum(x$res$FDR < fdr & x$res$logFC < 0, na.rm = TRUE)))
tp_lines <- if (nrow(tp_table)) sprintf("%s %s: lesion site %.2f vs non-lesional %.2f NPX, log2FC %+.2f (p = %s%s)", tp_table$Assay, tp_table$visit,
                                        tp_table$mean_NPX_lesion_site, tp_table$`mean_NPX_non-lesional`, tp_table$log2FC, fmt_p(tp_table$p),
                                        if_else(coalesce(tp_table$FDR < fdr, FALSE), sprintf(", FDR < %.2f", fdr), "")) else sprintf("%s: not in the data", paste(track, collapse = ", "))
sig_lines <- if (nrow(sig_nes)) sprintf("%s %s: NES %s (p = %s, FDR = %s)", sig_nes$set, sig_nes$visit, ifelse(is.na(sig_nes$NES), "-", sprintf("%+.2f", sig_nes$NES)), fmt_p(sig_nes$p), fmt_p(sig_nes$FDR)) else
  sprintf("%s: gene set not found", paste(sig_sets, collapse = ", "))
summary_lines <- c(
  sprintf("MicroAD dISF, %s. Visits %s. log2FC: relapse - non-relapse, and lesion site - non-lesional.", manifest_txt, paste(names(vn), collapse = ", ")),
  "## 1 / 2  Relapse vs non-relapse at the lesion site", rel_lines, sig_line,
  "## 3  Lesion site vs non-lesional skin (paired)", ln_lines,
  sprintf("## 4  %s, lesion site vs non-lesional", paste(track, collapse = ", ")), tp_lines,
  sprintf("## 5  %s, lesion site vs non-lesional (GSEA)", paste(map_chr(sig_sets, set_name), collapse = ", ")), sig_lines,
  "## Keep in mind",
  sprintf("Small groups (relapse vs non-relapse %s): the random-label medians show how many p < 0.05 proteins chance alone gives. Proteins that differ at several visits (same patients) make the stronger case.",
          paste(sprintf("%s %d vs %d", names(rel), map_int(rel, \(x) x$n[["relapse"]]), map_int(rel, \(x) x$n[["non-relapse"]])), collapse = ", ")))

source("R/report.R")
pdf_open(out_path(cfg, out, "Visits_V1-V4_summary.pdf"))
section("Summary", page_text("Relapse and lesional signatures per visit (V1-V4) - summary", c(summary_lines,
  "Figures: PNG files (300 dpi, 16:9, large fonts) can be pasted into PowerPoint as they are; the same as PDF."),
                             subtitle = sprintf("O-MicroAD, Olink Explore HT, dISF. Generated %s. Tables: 1_ 2_ 3_ 5_ Excel files in this folder.", format(Sys.Date())), size = 9.5, width = 145))
section("Relapse table", page_table("Relapse vs non-relapse at the lesion site, per visit",
  rel_summary |> transmute(visit, n = sprintf("%s/%s", relapse, `non-relapse`), tested = proteins_tested, `p<0.05` = p_below_0.05, higher = higher_in_relapse,
                           lower = lower_in_relapse, random_median = random_labels_median, random_P = fmt_p(random_labels_P), FDR = FDR_below,
                           sets_up = gene_sets_FDR_higher, sets_down = gene_sets_FDR_lower),
  subtitle = "n = relapse / non-relapse; random_median / random_P: p < 0.05 proteins under all random labellings (exact relabelling)",
  note = sprintf("FDR < %.2f; sets_up / sets_down: gene sets with FDR < %.2f higher / lower in relapsers", fdr, fdr)))
section("Lesional table", page_table("Lesion site vs non-lesional skin, per visit (paired)",
  ln_summary |> transmute(visit, pairs, lesion_site_state = str_trunc(lesion_site_state, 38), tested = proteins_tested, FDR_higher = FDR_higher_at_lesion_site, FDR_lower = FDR_lower_at_lesion_site),
  note = sprintf("FDR < %.2f; higher / lower at the lesion site", fdr)))
for (f in figs) section(f$title, page_plot(f$p + labs(title = f$title, subtitle = f$note)))
invisible(grDevices::dev.off())
msg("Per-visit analysis: %s", file.path(cfg$paths$output, out))
