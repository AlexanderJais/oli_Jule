# End-to-end test on simulated data: simulate -> run all steps -> check that the built-in
# effects are found and that false positives stay rare.
#   Rscript tests/test_pipeline.R
source("R/utils.R")

stopifnot(system2("Rscript", "tests/simulate_explore_ht.R") == 0)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
stopifnot(system2("Rscript", "run_all.R") == 0)

out   <- "output_sim"
truth <- read_csv("data_sim/truth.csv", show_col_types = FALSE)
rd <- \(f) read_csv(file.path(out, f), show_col_types = FALSE) |> left_join(truth |> select(OlinkID, role), by = "OlinkID")
hits <- \(r, mdl, ct, rl) with(r |> filter(model == mdl, contrast == ct), c(found = sum(significant & role == rl), of = sum(role == rl)))
fp   <- \(r, mdl, ct) with(r |> filter(model == mdl, contrast == ct), mean(significant[role == "null"]))

isf <- rd("models/ISF_results.csv"); ser <- rd("models/Serum_results.csv")
delta <- rd("models/ISF_relapse_delta_results.csv") |> mutate(model = "delta")
cor <- rd("isf_serum/isf_serum_correlation.csv")
checks <- list(
  "ISF lesional vs non-lesional"   = hits(isf, "states_all_visits", "AD_L_vs_NL", "ISF_lesional"),
  "ISF AD non-lesional vs healthy" = hits(isf, "states_all_visits", "AD_NL_vs_HC", "ISF_AD_vs_HC"),
  "ISF relapse (xL - NL delta)"    = hits(delta, "delta", "relapse_vs_non", "ISF_relapse"),
  "Serum AD vs in-study HC"        = hits(ser, "AD_vs_HC_in_study", "AD_vs_HC", "Serum_AD_vs_HC"),
  "Serum RELAD relapse"            = hits(ser, "RELAD_relapse", "relapse_vs_non", "Serum_relapse"),
  "ISF-serum coupling (L site)"    = with(cor |> filter(site == "L"), c(found = sum(significant & role == "ISF_serum_coupled"), of = sum(role == "ISF_serum_coupled")))
)
res <- tibble(check = names(checks), found = map_dbl(checks, "found"), of = map_dbl(checks, "of")) |>
  mutate(power = found / of)
print(res)
fpr <- c(ISF = fp(isf, "states_all_visits", "AD_L_vs_NL"), Serum = fp(ser, "AD_vs_HC_in_study", "AD_vs_HC"))
print(round(fpr, 3))

biobank <- ser |> filter(model == "AD_vs_HC_in_study", role == "Biobank_shift")
leip <- read_csv(file.path(out, "leip_reference/leip_clinical_associations.csv"), show_col_types = FALSE)
th2 <- read_csv(file.path(out, "enrichment/gsea_results.csv"), show_col_types = FALSE) |>
  filter(pathway == "CUSTOM_AD_TH2_AXIS", model == "states_all_visits", contrast == "AD_L_vs_NL")

stopifnot(
  "an effect was not recovered (power < 0.7)" = all(res$power >= 0.7),
  "too many false positives" = all(fpr <= 0.05),
  "biobank-only shift leaked into the in-study comparison" = sum(biobank$significant) <= 1,
  "BMI association in LEIP not found" = any(leip$significant & leip$parameter == "BMI"),
  "Th2 set not enriched in lesional skin" = nrow(th2) == 1 && th2$padj < 0.05 && th2$NES > 0
)
msg("All pipeline checks passed.")
