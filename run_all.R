# Run the whole pipeline in order.
#   Rscript run_all.R                 # real data (config.yml)
#   OLINK_CONFIG=data_sim/config_sim.yml Rscript run_all.R   # simulated data
steps <- c("scripts/01_metadata.R", "scripts/02_import_qc.R", "scripts/03_explore.R",
           "scripts/04_isf_models.R", "scripts/05_serum_models.R", "scripts/06_isf_vs_serum.R",
           "scripts/07_enrichment.R", "scripts/08_leip_reference.R")
for (s in steps) {
  message("\n==== ", s, " ====")
  status <- system2("Rscript", s)
  if (status != 0) stop(s, " failed (exit ", status, ")")
}
message("\nAll steps finished. Results are in the output folder set in the config.")
