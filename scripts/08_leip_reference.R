# 08 - LEIP population reference for proteins correlated between ISF and serum
# For each protein that is significantly correlated between ISF and serum (script 06):
#   1. normal range in the population (LEIP) and where AD patients fall relative to it
#   2. associations with clinical parameters in healthy people (possible confounders)
#   3. detectability in the population
#   4. source check: LEIP vs in-study healthy controls (pre-analytical differences)
# Out: output/leip_reference/*.csv and leip_reference.xlsx

source("R/utils.R")
cfg   <- load_config()
clear_outputs(cfg, "leip_reference")
meta  <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
sig   <- read_csv(file.path(cfg$paths$output, "isf_serum", "significant_proteins.csv"), show_col_types = FALSE)
if (!nrow(sig)) { msg("No ISF-serum correlated proteins - nothing to check in LEIP."); quit(save = "no") }

serum <- clean |> filter(matrix == "Serum", OlinkID %in% sig$OlinkID)
leip  <- serum |> filter(cohort == "LEIP")
if (!nrow(leip)) stop("No LEIP samples in the cleaned serum data.")

# 1 + 3: reference ranges and detectability -----------------------------------------------------
ref <- leip |>
  group_by(OlinkID, Assay) |>
  summarise(n_leip = sum(!is.na(value)), leip_mean = mean(value, na.rm = TRUE), leip_sd = sd(value, na.rm = TRUE),
            leip_median = median(value, na.rm = TRUE),
            leip_p05 = quantile(value, 0.05, na.rm = TRUE), leip_p95 = quantile(value, 0.95, na.rm = TRUE),
            leip_frac_detected = mean(!below_lod, na.rm = TRUE), .groups = "drop")

pos <- serum |>
  filter(group %in% c("AD", "HC"), cohort != "LEIP") |>
  left_join(ref, by = c("OlinkID", "Assay")) |>
  mutate(z = (value - leip_mean) / leip_sd,
         above_p95 = value > leip_p95, below_p05 = value < leip_p05)
save_csv(pos |> select(SampleID, SubjectID, cohort, group, visit, lesion_state, OlinkID, Assay, value, z,
                       above_p95, below_p05), cfg, "leip_reference", "samples_vs_leip.csv")
pos_sum <- pos |>
  group_by(OlinkID, Assay, group) |>
  summarise(n = sum(!is.na(z)), median_z = median(z, na.rm = TRUE), pct_above_p95 = 100 * mean(above_p95, na.rm = TRUE),
            pct_below_p05 = 100 * mean(below_p05, na.rm = TRUE),
            .groups = "drop") |>
  pivot_wider(names_from = group, values_from = c(n, median_z, pct_above_p95, pct_below_p05))

# 4: LEIP vs in-study healthy controls ---------------------------------------------------------------
hc <- serum |> filter(group == "HC")
src <- map(unique(serum$OlinkID), \(a) {
  x <- na.omit(hc$value[hc$OlinkID == a]); y <- na.omit(leip$value[leip$OlinkID == a])
  if (length(x) < 3 || length(y) < 3) return(NULL)
  tibble(OlinkID = a, hc_minus_leip = median(x) - median(y),
         p_source = suppressWarnings(wilcox.test(x, y, exact = FALSE))$p.value)
}) |> bind_rows() |> mutate(fdr_source = p.adjust(p_source, "BH"), source_shift = fdr_source < cfg$stats$fdr)

# 2: clinical associations in LEIP ------------------------------------------------------------------
params <- intersect(c("age", "BMI", "WHR", "c_fett", "HOMA_IR", "c_CRP", "C_CHOL", "C_HDL", "C_LDL",
                      "C_TRIGLY", "c_apo", "MDRD_kurz", "Gluc0_mg_dl"), names(leip))
assoc <- NULL
if (length(params)) {
  assoc_num <- map(params, \(pm) leip |> group_by(OlinkID, Assay) |>
                     filter(sum(!is.na(.data[[pm]])) >= 10) |>
                     summarise(parameter = pm, n = sum(!is.na(.data[[pm]])),
                               rho = suppressWarnings(cor(value, .data[[pm]], method = "spearman", use = "complete.obs")),
                               p = suppressWarnings(cor.test(value, .data[[pm]], method = "spearman", exact = FALSE))$p.value,
                               .groups = "drop")) |> bind_rows()
  assoc_sex <- if ("sex" %in% names(leip)) leip |> filter(!is.na(sex)) |> group_by(OlinkID, Assay) |>
    summarise(parameter = "sex (M - F)", n = n(),
              rho = median(value[sex == "M"], na.rm = TRUE) - median(value[sex == "F"], na.rm = TRUE),
              p = suppressWarnings(wilcox.test(value[sex == "M"], value[sex == "F"], exact = FALSE))$p.value,
              .groups = "drop") else NULL
  assoc <- bind_rows(assoc_num, assoc_sex) |>
    mutate(fdr = p.adjust(p, "BH"), significant = fdr < cfg$stats$fdr) |>
    arrange(p)
  save_csv(assoc, cfg, "leip_reference", "leip_clinical_associations.csv")
  if (nrow(assoc)) {
    p <- assoc |> filter(parameter != "sex (M - F)") |>
      ggplot(aes(parameter, Assay, fill = rho)) + geom_tile() +
      geom_text(aes(label = if_else(significant, "*", "")), size = 5) +
      scale_fill_gradient2(limits = c(-1, 1)) +
      labs(title = "LEIP: Spearman correlation with clinical parameters (* FDR < 0.05)", x = NULL, y = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_plot(p, cfg, "leip_reference", "leip_clinical_heatmap.png", width = 9, height = 2 + 0.25 * n_distinct(assoc$Assay))
  }
}

# summary table: one row per protein --------------------------------------------------------------
top_assoc <- if (!is.null(assoc)) assoc |> filter(significant) |> group_by(OlinkID) |>
  summarise(clinical_associations = paste(sprintf("%s (%.2f)", parameter, rho), collapse = "; ")) else
  tibble(OlinkID = character(), clinical_associations = character())
corr <- read_csv(file.path(cfg$paths$output, "isf_serum", "isf_serum_correlation.csv"), show_col_types = FALSE) |>
  filter(OlinkID %in% sig$OlinkID) |>
  select(OlinkID, site, r_within, fdr_within, r_between, fdr_between) |>
  pivot_wider(names_from = site, values_from = c(r_within, fdr_within, r_between, fdr_between))
summary_tbl <- ref |> left_join(corr, by = "OlinkID") |> left_join(pos_sum, by = c("OlinkID", "Assay")) |>
  left_join(src, by = "OlinkID") |> left_join(top_assoc, by = "OlinkID") |>
  arrange(desc(abs(coalesce(median_z_AD, 0))))
save_csv(summary_tbl, cfg, "leip_reference", "leip_reference_summary.csv")
writexl::write_xlsx(list(summary = summary_tbl, clinical_associations = assoc %||% tibble(), source_check = src),
                    out_path(cfg, "leip_reference", "leip_reference.xlsx"))
msg("%d proteins checked in LEIP (n = %d); %d show a LEIP vs in-study control shift; %d have clinical associations",
    nrow(ref), n_distinct(leip$SampleID), sum(src$source_shift, na.rm = TRUE), nrow(top_assoc))
