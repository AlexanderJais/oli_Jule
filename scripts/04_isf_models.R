# 04 - ISF (dermal interstitial fluid) models, MicroAD
# Per protein; subject as random effect (dream), plate as fixed effect in every model.
# Out: output/models/ISF_results.csv, ISF_summary.csv, volcano plots

source("R/utils.R")
source("R/models.R")
source("R/design.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")
clean <- read_step(cfg, "data", "npx_clean.rds", step = "scripts/02_import_qc.R")
assay_map <- clean |> distinct(OlinkID, Assay)
expr <- wide$ISF
clear_outputs(cfg, "models", "^ISF_")

info <- isf_design(meta) |> filter(SampleID %in% colnames(expr))
specs <- isf_specs(info)

res <- run_model_specs(specs, expr, info, cfg, "ISF", assay_map)

# Relapse on the within-visit difference ex-lesional minus non-lesional (same subject, visit, plate:
# removes plate and systemic day-to-day variation).
# Pivot on subject + visit only, so a pair is kept even if its two samples differ in plate or date;
# covariates are taken from the ex-lesional sample.
pairs <- isf_delta_pairs(info)
if (nrow(pairs) >= 4) {
  delta <- expr[, pairs$AD_xL, drop = FALSE] - expr[, pairs$AD_NL, drop = FALSE]
  colnames(delta) <- pairs$SampleID
  rd <- fit_contrasts(delta, pairs, ~ 0 + relapse2 + weeks + (1 | SubjectID),
                      c(relapse_vs_non = "relapse2relapse - relapse2non_relapse"),
                      "relapse_delta_xL_minus_NL", cfg$stats$min_group_n)
  if (!is.null(rd)) {
    rd <- annotate_results(rd, assay_map, cfg$stats$fdr)
    save_csv(rd, cfg, "models", "ISF_relapse_delta_results.csv")
    msg("relapse_delta_xL_minus_NL: %d pairs, %d significant", nrow(pairs), sum(rd$significant))
  }
}
