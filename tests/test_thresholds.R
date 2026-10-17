# Test of the threshold analysis (run_thresholds.R) on simulated data: every threshold is listed with a verdict, fewer
# proteins are kept and tested the stricter the detection filter, the simulated relapse markers pass the V1 relapse
# criteria, group sizes and the smallest exact p follow the simulated design, and all files are written.
#   Rscript tests/test_thresholds.R
source("R/utils.R")
rscript <- file.path(R.home("bin"), "Rscript")
stopifnot(system2(rscript, "tests/make_synthetic_manifest.R") == 0)
stopifnot(system2(rscript, "tests/simulate_explore_ht.R") == 0)
old_cfg <- Sys.getenv("OLINK_CONFIG", NA)
Sys.setenv(OLINK_CONFIG = "data_sim/config_sim.yml")
run_ok <- system2(rscript, "run_thresholds.R") == 0
if (is.na(old_cfg)) Sys.unsetenv("OLINK_CONFIG") else Sys.setenv(OLINK_CONFIG = old_cfg)
stopifnot("threshold analysis failed on simulated data" = run_ok)

out <- "output_sim/thresholds"
rd <- \(sh) readxl::read_excel(file.path(out, "thresholds.xlsx"), sh)
verdict <- rd("verdict"); inv <- rd("inventory"); det <- rd("detection"); cmp <- rd("comparisons"); grp <- rd("group_sizes")
funnel <- rd("relapse_funnel"); rem <- rd("removed_by_filter")
cfg_sim <- yaml::read_yaml("data_sim/config_sim.yml")
mono <- \(x) all(diff(x) <= 0)
stopifnot(
  "verdict incomplete" = all(c("Detection (LOD) filter", "False discovery rate", "Minimum group size", "Olink sample QC: FAIL", "Samples not in the manifest",
                               "Missing values in the main models") %in% verdict$threshold) && !any(is.na(verdict$verdict)) &&
    all(verdict$status %in% c("may limit these data", "by design", "not limiting")),
  "inconsistencies not listed" = nrow(rd("inconsistencies")) >= 10,
  "removed samples without the share of failed values" = "pct_values_FAIL" %in% names(rd("removed_samples")),
  "inventory missing the config thresholds" = all(c("qc$min_detect_frac", "stats$fdr", "stats$min_group_n", "stats$visit_min_subjects") %in% inv$config_key),
  "current detection threshold not marked" = sum(det$current) == 2 && all(det$detection_threshold[det$current] == cfg_sim$qc$min_detect_frac),
  "assays kept must fall with a stricter filter" = all(det |> group_by(matrix) |> summarise(ok = mono(kept_in_one_group) && mono(kept_of_all_samples)) |> pull(ok)),
  "proteins tested must fall with a stricter filter" = all(cmp |> filter(!is.na(proteins_tested)) |> group_by(comparison) |> summarise(ok = mono(proteins_tested)) |> pull(ok)),
  "five comparisons expected" = n_distinct(cmp$comparison) == 5 && nrow(rem) == 5,
  "relapse chance level missing" = all(!is.na(cmp$chance_median_p_below_0.05[cmp$short == "relapse V1"])),
  "simulated relapse markers lost by the V1 criteria" = funnel$proteins[funnel$direction == "higher in relapsers" & grepl("above LOD", funnel$criterion) & funnel$kind == "screen"] >= 4,
  "group sizes wrong (simulation: 4 relapsers, 6 non-relapsers at V1)" = grp$relapse[grp$visit == "V1"] == 4 && grp$non_relapse[grp$visit == "V1"] == 6 &&
    abs(grp$smallest_exact_p[grp$visit == "V1"] - 1 / choose(10, 4)) < 1e-9,
  "V4 relapse comparison should be below the minimum group size (2 relapsers)" = !grp$relapse_comparison_tested[grp$visit == "V4"],
  "files missing" = all(file.exists(file.path(out, c("Thresholds_summary.pdf", "1_assays_kept_by_detection_threshold.png", "2_comparisons_by_detection_threshold.png",
                                                     "3_removed_proteins_vs_chance.png", "5_group_sizes_per_visit.png", "6_relapse_criteria_V1.png"))))
)
msg("All threshold checks passed.")
