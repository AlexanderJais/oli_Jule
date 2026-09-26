# 09 - Aim 1: descriptive proteome profile of dISF
# Which proteins are reliably measurable in dermal ISF, how this differs by skin state,
# what biology the detectable ISF proteome covers, and which proteins vary most.
# NPX is a relative, protein-specific scale, so levels are NOT compared between proteins;
# the profile uses detection rates, variability and pathway composition instead.
# Out: output/isf_profile/*  (incl. isf_profile.xlsx)

source("R/utils.R")
cfg   <- load_config()
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
wide  <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
min_f <- cfg$qc$min_detect_frac

isf <- clean |> filter(matrix == "ISF") |>
  mutate(state_group = case_when(group == "HC" ~ "healthy (HC)",
                                 group == "CPUO" ~ paste("CPUO", state),
                                 TRUE ~ paste("AD", state)))

# ---- detection per protein and skin state --------------------------------------------------------
det_state <- isf |>
  group_by(OlinkID, Assay, UniProt, state_group) |>
  summarise(frac = mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = state_group, values_from = frac, names_prefix = "detected: ")
det_all <- isf |>
  group_by(OlinkID, Assay) |>
  summarise(n_samples = n_distinct(SampleID), frac_detected = mean(!below_lod, na.rm = TRUE),
            sd_npx = sd(value, na.rm = TRUE), .groups = "drop")
ad_cols <- c("detected: AD lesional", "detected: AD non-lesional", "detected: healthy (HC)")
profile <- det_all |> left_join(det_state, by = c("OlinkID", "Assay")) |>
  mutate(
    detection_class = case_when(frac_detected >= 0.9 ~ "robust (>= 90%)",
                                frac_detected >= min_f ~ sprintf("detected (%.0f-90%%)", 100 * min_f),
                                frac_detected > 0.1 ~ "sporadic (10-50%)",
                                TRUE ~ "not detected (<= 10%)")) |>
  arrange(desc(frac_detected))
# detectable in lesional AD skin but not in non-lesional or healthy skin
profile$lesion_restricted <- if (all(ad_cols %in% names(profile)))
  profile[[ad_cols[1]]] >= min_f & profile[[ad_cols[2]]] < min_f & profile[[ad_cols[3]]] < min_f else NA
save_csv(profile, cfg, "isf_profile", "isf_detection_profile.csv")
msg("dISF: %s", paste(names(table(profile$detection_class)), table(profile$detection_class), collapse = ", "))
msg("%d proteins detectable only in lesional AD skin", sum(profile$lesion_restricted, na.rm = TRUE))

by_state <- isf |> group_by(state_group, SampleID) |>
  summarise(n_detected = sum(!below_lod, na.rm = TRUE), .groups = "drop")
p <- ggplot(by_state, aes(reorder(state_group, n_detected, median), n_detected)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") + geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
  coord_flip() + labs(title = "Proteins above LOD per dISF sample", x = NULL, y = "proteins detected")
save_plot(p, cfg, "isf_profile", "detected_per_sample.png", width = 8, height = 4)

# ---- pathway composition of the detectable ISF proteome (over-representation vs all assays) ----
detected_genes <- profile |> filter(frac_detected >= min_f) |> pull(Assay) |> str_split("_") |> unlist() |> unique()
universe <- clean |> distinct(Assay) |> pull(Assay) |> str_split("_") |> unlist() |> unique()
sets <- map(cfg$enrichment$collections, \(cl) {
  parts <- str_split_fixed(cl, ":", 2)
  g <- if (parts[2] == "") msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1])
       else msigdbr::msigdbr(species = "Homo sapiens", collection = parts[1], subcollection = parts[2])
  split(g$gene_symbol, g$gs_name)
}) |> unlist(recursive = FALSE)
ora <- fgsea::fora(sets, genes = detected_genes, universe = universe,
                   minSize = cfg$enrichment$min_size, maxSize = cfg$enrichment$max_size) |>
  as_tibble() |> mutate(overlapGenes = map_chr(overlapGenes, paste, collapse = ";")) |>
  arrange(pval)
save_csv(ora, cfg, "isf_profile", "isf_detected_pathways.csv")
msg("Detectable dISF proteome: %d of %d panel genes; %d gene sets over-represented (padj < %.2f)",
    length(detected_genes), length(universe), sum(ora$padj < cfg$stats$fdr), cfg$stats$fdr)

# ---- most variable proteins across dISF samples (heatmap of z-scores) --------------------------
m <- wide$ISF
if (!is.null(m) && nrow(m) > 1) {
  top <- names(sort(apply(m, 1, sd, na.rm = TRUE), decreasing = TRUE))[seq_len(min(50, nrow(m)))]
  info <- isf |> distinct(SampleID, state_group, SubjectID, visit_num)
  z <- t(scale(t(m[top, , drop = FALSE])))
  hm <- as_tibble(z, rownames = "OlinkID") |>
    pivot_longer(-OlinkID, names_to = "SampleID", values_to = "z") |>
    left_join(info, by = "SampleID") |>
    left_join(distinct(clean, OlinkID, Assay), by = "OlinkID") |>
    mutate(z = pmax(pmin(z, 3), -3),
           SampleID = factor(SampleID, levels = info |> arrange(state_group, SubjectID, visit_num) |> pull(SampleID)),
           Assay = factor(Assay, levels = rev(unique(Assay[order(match(OlinkID, top))]))))
  p <- ggplot(hm, aes(SampleID, Assay, fill = z)) + geom_tile() +
    scale_fill_gradient2(low = "steelblue", high = "firebrick") +
    facet_grid(~state_group, scales = "free_x", space = "free_x") +
    labs(title = "50 most variable dISF proteins (z-score per protein)", x = "samples", y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size = 7))
  save_plot(p, cfg, "isf_profile", "top_variable_heatmap.png", width = 12, height = 9)
}

writexl::write_xlsx(list(detection_profile = profile, detected_pathways = ora),
                    out_path(cfg, "isf_profile", "isf_profile.xlsx"))
