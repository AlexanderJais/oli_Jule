# O-MicroAD – Olink Explore HT analysis

R pipeline for the Olink Explore HT run of O-MicroAD: dermal interstitial fluid (dISF) and serum
from the MicroAD study (atopic dermatitis patients followed over up to 6 visits, healthy volunteers,
CPUO), plus serum from RELAD / RELAD2 and LEIP biobank controls.

## Study aims and where they are answered

| aim | steps | key outputs |
|---|---|---|
| 1. Profile the dISF proteome of mild-to-moderate AD | 09, 04, 07 | `isf_profile/isf_profile.xlsx` (what is detectable, by skin state, pathways), `models/ISF_*` (lesional / ex-lesional / non-lesional / healthy) |
| 2. Compare the dISF and blood proteome | 10, 06, 08 | `matrix_comparison/matrix_comparison.xlsx` (detected where, relatively enriched in dISF, skin vs blood disease signals), `isf_serum/*` (correlation), `leip_reference/*` |
| 3. Track changes across the disease course | 11, 04 (time models) | `trajectories/trajectory_results.csv`, per-patient plots in `trajectories/plots/` |

| Per visit / time course (lesion site vs healthy skin at V1–V6; proteins regulated at all visits) | 13 | `visit_course/visit_course.xlsx`, `visit_course/time_course_*.png`, `models/volcano/ISF_by_visit_V*.png` |
| **Executive summary of everything** | 14 | `Executive_summary.pdf` |
| Focus: CD137 (4-1BB, assay TNFRSF9) | 12 | `focus/TNFRSF9/TNFRSF9_report.xlsx` + figures |

Secondary: relapse (04, 05, 11; exploratory) and RELAD/RELAD2 serum relapse (05).

ISF and serum are analysed separately. They are only combined for the ISF–serum correlation
(step 06) and the LEIP population check (step 08).

## Quick start

Open `oli_Jule.Rproj` in RStudio (or `setwd()` to this folder): all paths are relative to the project folder.

```bash
Rscript install_packages.R        # once
# put the data files in data/ (see below), then
Rscript run_all.R                 # runs steps 01-14; results in output/, summary in output/Executive_summary.pdf
```

To run a single step, use `Rscript scripts/0X_....R`. Each step reads what the previous one saved.

## Data (never committed – `data/` and `output/` are git-ignored)

| file | what | where it is set |
|---|---|---|
| `data/manifest.xlsx` | Olink sample submission sheet; the `manifest` sheet is the master | `paths$manifest` |
| `data/npx/*.parquet` | Olink NPX files, here `O-MicroAD_ISF_NPX_2026-09-24.parquet` and `O-MicroAD_Serum_NPX_2026-09-24.parquet`; all files in the folder are read, and the file name must contain ISF or Serum | `paths$npx_dir` |
| `data/LEIP_clinical_parameters_n35.xlsx` | LEIP clinical data, sheet `Key_parameters` (optional) | `paths$leip_clinical` |
| `data/severity.xlsx` | optional: `SubjectID`, `Visit` (V1–V6), plus numeric scores (e.g. SCORAD, EASI, itch NRS). Step 11 then models each protein against each score | `paths$severity` |
| `data/Explore_HT_Fixed_LOD.csv` | Olink fixed LOD file for Explore HT, version ≥ 6.0.0, from olink.com (recommended) | `paths$fixed_lod` |

All settings (thresholds, FDR, number of cores) are in `config.yml`.

## Steps

| script | does | main output |
|---|---|---|
| `01_metadata.R` | Cleans the manifest: harmonised groups (AD / HC / CPUO / Biobank), skin site and state, visit and days since V1, relapse visit and visits before relapse, sample volume, and LEIP clinical data. Checks IDs against the plate layout, and flags inconsistencies without dropping anything. | `metadata/sample_metadata.csv`, `metadata/data_flags.csv` |
| `02_import_qc.R` | Reads the parquet, checks sample IDs and normalisation, and works out LOD. Sample QC uses Olink flags plus a median/IQR outlier check. Keeps an assay in a matrix if it is detected in ≥ 50 % of at least one group, and reports control CVs and low-volume samples. | `qc/*`, `data/npx_clean.rds`, `data/npx_wide.rds` |
| `03_explore.R` | Per matrix: PCA coloured by state, group, visit, plate and cohort; variance partitioning. | `explore/*` |
| `04_isf_models.R` | ISF models per protein (see below). | `models/ISF_*` |
| `05_serum_models.R` | Serum models per protein (see below). | `models/Serum_*` |
| `06_isf_vs_serum.R` | On matched MicroAD visits: within-subject (repeated-measures) and between-subject correlation of ISF and serum, separately for the lesional and non-lesional site. | `isf_serum/*` |
| `07_enrichment.R` | GSEA (fgsea) for every contrast: MSigDB Hallmark, Reactome, GO:BP, plus a custom AD/Th2 set. | `enrichment/*` |
| `09_isf_profile.R` | **Aim 1.** Descriptive dISF profile: detection class per protein (overall and by skin state), proteins detectable only in lesional skin, pathway over-representation of the detectable proteome, heatmap of the most variable proteins. NPX is protein-specific, so there is no ranking of levels between proteins. | `isf_profile/*` |
| `10_matrix_comparison.R` | **Aim 2.** Detected in dISF only / serum only / both. Relative dISF/serum enrichment: paired, centred log2 ratio, FDR plus ≥ 2-fold (`stats$min_rel_log2`). Concordance of disease effects in skin vs blood. Whether serum tracks the lesional-minus-non-lesional skin difference visit by visit. | `matrix_comparison/*` |
| `11_trajectories.R` | **Aim 3.** Residual lesional signal (ex-lesional minus non-lesional) vs weeks since clearance, and vs weeks to relapse (relapsers; exploratory). Serum vs weeks to relapse. Optional severity models, and per-patient trajectory plots of the top proteins (dISF lesional site, non-lesional site, serum; relapse marked). | `trajectories/*` |
| `12_focus_proteins.R` | **Dedicated analysis of pre-specified proteins** (`focus_proteins` in `config.yml`; default CD137 = TNFRSF9). Runs even if the protein fails the detection filter. Reports detection per matrix and group, and the same models as steps 04/05 plus the xL − NL relapse model for this protein alone (lmerTest / lm). The unadjusted p-value is the primary test, with the proteome-wide FDR shown alongside. Also ISF–serum correlation, severity (if available), LEIP clinical associations, every pipeline result for the protein, and figures. | `focus/<protein>/*` |
| `13_visit_course.R` | dISF per visit (V1–V6, visits with ≥ `stats$visit_min_subjects` patients): tracked lesion site vs healthy skin, vs non-lesional skin, and non-lesional vs healthy. Same model as step 04, with a volcano plot per visit. Proteins regulated at **all** visits in the same direction are listed at FDR < 0.05 at every visit (strict) and at p < 0.05 at every visit (nominal), with time-course plots of effects and NPX levels. | `visit_course/*`, `models/ISF_by_visit_*` |
| `14_summary_report.R` | Executive summary PDF: data and QC, automatically extracted key findings per aim, all comparisons, visit course, volcano plots, pathways, dISF vs serum, serum, focus proteins, methods and caveats. Each section is skipped with a note if its step did not run. | `Executive_summary.pdf` |
| `08_leip_reference.R` | For the proteins significant in step 06: LEIP normal range, where AD patients fall in it, clinical associations in LEIP, detectability, and LEIP vs in-study controls. | `leip_reference/*` (incl. `.xlsx`) |

### Models

The models are limma/dream (variancePartition) with empirical Bayes moderation. Subject is a random
effect wherever a person contributes several samples. NPX differences are on the log2 scale, and
BH FDR is applied within each contrast.

**ISF (04)**, all models adjusted for plate:
- `states_all_visits`: lesional, ex-lesional and non-lesional AD skin vs each other and vs healthy skin; CPUO lesional vs non-lesional.
- `baseline_V1`: the same comparisons at V1 only.
- `time_ex_lesional`, `time_non_lesional`: change per week in AD skin.
- `relapse_ex_lesional`, exploratory: cleared skin of relapsers vs non-relapsers.
- `relapse_delta_xL_minus_NL`, exploratory: the same question on the within-visit ex-lesional minus non-lesional difference. Plate and systemic day-to-day variation cancel out.

**Serum (05):**
- AD vs controls, fitted twice: against the in-study healthy controls (adjusted for cohort and plate) and against LEIP biobank controls (adjusted for plate).
  - `Serum_AD_vs_controls_agreement.csv` marks proteins that agree in both comparisons.
  - `HC_vs_Biobank` shows proteins affected by the biobank source.
- `MicroAD_active_vs_cleared`: serum when the tracked lesion is active vs cleared.
- `MicroAD_relapse`: relapsers vs non-relapsers at cleared visits.
- `RELAD_relapse`: RELAD and RELAD2, adjusted for cohort and plate, plus a sensitivity analysis without the samples with conflicting relapse labels.

## Design decisions

- **`PCNormalizedNPX` is analysed, not intensity-normalised NPX.** Plate 1 is ISF only, plate 2 is mixed, and plates 3–4 are serum. Intensity normalisation assumes randomised samples of one matrix and would distort plate 2. Step 02 reports the `Normalization` column of the delivered file.
- **LOD** is computed per row with OlinkAnalyze's own routine (`olink_lod`), so count-based assays (`LODMethod = lod_count`, about 18 % in the fixed LOD file v10.2.0) get their sample-specific LOD.
  - Preferred: Olink's fixed LOD file (`paths$fixed_lod`), matched on `DataAnalysisRefID`.
  - Fallback for rows without a match: the Olink negative-control method on the 8 negative controls. OlinkAnalyze normally requires ≥ 10, so this LOD is less precise; in simulation it was off by up to ±0.6 NPX.
  - `qc/lod.csv` gives the LOD source and the per-sample LOD range for each assay.
- **Relapse in ISF is exploratory.** All four ISF relapsers are on plate 1, and plate 2 holds only non-relapsers.
- **Metadata issues are flagged, not fixed**, in `metadata/data_flags.csv`. As of manifest v3 this covers the RELAD / RELAD2 label conflicts, 5 low-volume ISF samples, and LEIP_35 without clinical data.
- **Age and sex** are currently only available for LEIP, so the serum models are not adjusted for them.

## Testing without real data

```bash
Rscript tests/test_pipeline.R
```

This needs no study data. It builds a synthetic manifest with the same layout and design, simulates the two Olink NPX files (ISF, serum) with known effects, runs all
steps with `data_sim/config_sim.yml` (output in `output_sim/`), and checks that the effects are
recovered and false positives stay rare.
