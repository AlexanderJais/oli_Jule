# O-MicroAD – technical notes

R pipeline for the Olink Explore HT run of O-MicroAD: dermal interstitial fluid (dISF) and serum
from the MicroAD study (atopic dermatitis patients followed over up to 6 visits, healthy volunteers,
CPUO), plus serum from RELAD / RELAD2 and LEIP biobank controls.

## Study aims and where they are answered

| aim / topic | steps | key outputs |
|---|---|---|
| 1. Profile the dISF proteome of mild-to-moderate AD | 09, 04, 07 | `isf_profile/isf_profile.xlsx` (what is detectable, by skin state, pathways), `models/ISF_*` (lesional / ex-lesional / non-lesional / healthy) |
| 2. Compare the dISF and blood proteome | 10, 06, 08, 14 | `matrix_comparison/matrix_comparison.xlsx` (detected where, relatively enriched in dISF per skin type, skin vs blood disease signals), `isf_serum/*` (correlation), `leip_reference/*`, `serum_vs_disf/*` (overlap per visit) |
| 3. Track changes across the disease course | 11, 13, 04 (time models) | `trajectories/trajectory_results.csv`, per-patient plots in `trajectories/plots/`, `visit_course/*` (per visit, proteins regulated at all visits) |
| Focus proteins: CD137 (TNFRSF9), TNFSF9, KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, PNOC, POSTN | 12 | `focus/focus_overview.xlsx`, `focus/focus_overview_heatmap.png`, one folder per protein |
| Key questions: mast cells, CD137 / CD137L, relapse, dISF vs serum | 15 | `key_questions/answers.csv`, `key_questions/key_questions.xlsx`, `key_questions/evidence/*.png` |
| **Executive summary of everything** | 16 | `Executive_summary.pdf`, `Executive_summary_tables.xlsx` (full lists) |
| **Data export: all proteins, all samples (CSV); RELAD2 separately** | 17 | `export/` |

Relapse analyses (04, 05, 11, 14, 15) are exploratory in MicroAD; RELAD/RELAD2 serum gives larger groups (05).

ISF and serum are modelled separately. They are combined only where the two matrices are compared:
correlation (06), LEIP check (08), relative enrichment and concordance (10), overlap (14) and Q5 (15).

## Quick start

Open `oli_Jule.Rproj` in RStudio (or `setwd()` to this folder): all paths are relative to the project folder.

```bash
Rscript install_packages.R        # once
# put the data files in data/ (see below), then
Rscript run_all.R                 # runs all steps; results in output/, summary in output/Executive_summary.pdf
```

To run a single step, use `Rscript scripts/NN_....R` (or `start_at <- NN; source("run_all.R")` in RStudio to resume from there). Each step reads what the earlier ones saved.

## Data (never committed – `data/` and `output/` are git-ignored)

| file | what | where it is set |
|---|---|---|
| `data/manifest.xlsx` | Olink sample submission sheet; the `manifest` sheet is the master (v4: incl. `Sex`, `NoRELAD2`) | `paths$manifest` |
| `data/npx/*.parquet` | Olink NPX files, here `O-MicroAD_ISF_NPX_2026-09-24.parquet` and `O-MicroAD_Serum_NPX_2026-09-24.parquet`; all files in the folder are read, and the file name must contain ISF or Serum | `paths$npx_dir` |
| `data/LEIP_clinical_parameters_n35.xlsx` | LEIP clinical data, sheet `Key_parameters` (optional) | `paths$leip_clinical` |
| `data/severity.xlsx` | optional: `SubjectID`, `Visit` (V1–V6 or 1–6), plus numeric scores (e.g. SCORAD, EASI, itch NRS). Steps 11 and 12 then model the proteins against each score | `paths$severity` |
| `data/Explore_HT_Fixed_LOD.csv` | Olink fixed LOD file for Explore HT, version ≥ 6.0.0, from olink.com (recommended) | `paths$fixed_lod` |

If a file is not at the configured path, a single file with the usual name in the same folder is used instead
(e.g. `Explore HT_Fixed LOD.csv`, `*Sample Submission Sheet*.xlsx`). `run_all.R` first lists what it found.
All settings (thresholds, FDR, focus proteins, mast cell markers, number of cores) are in `config.yml`.

## Steps

`run_all.R` runs the steps in this order; `start_at <- N; source("run_all.R")` resumes at step N.

| script | does | main output |
|---|---|---|
| `01_metadata.R` | Cleans the manifest: harmonised groups (AD / HC / CPUO / Biobank), skin site and state, clinical state, sex, visit and days since V1, relapse visit and visits before relapse, sample volume, and LEIP clinical data. Checks IDs against the plate layout, and flags inconsistencies without dropping anything. | `metadata/sample_metadata.csv`, `metadata/data_flags.csv` |
| `02_import_qc.R` | Reads the parquet files, checks sample IDs, matrix and normalisation, and works out the LOD. Sample QC uses Olink flags plus a median/IQR outlier check. Removes assays Olink excluded (no values) and keeps an assay in a matrix if it is detected in ≥ 50 % of at least one group. Reports control CVs and low-volume samples. | `qc/*`, `data/npx_clean.rds`, `data/npx_wide.rds`, `data/npx_all_samples.rds` |
| `03_explore.R` | Per matrix: PCA coloured by state, group, visit, plate and cohort; variance partitioning. | `explore/*` |
| `04_isf_models.R` | ISF models per protein (see below). | `models/ISF_*` |
| `05_serum_models.R` | Serum models per protein (see below). | `models/Serum_*` |
| `06_isf_vs_serum.R` | On matched MicroAD visits: within-subject (repeated-measures) and between-subject correlation of ISF and serum, separately for the lesional and non-lesional site. | `isf_serum/*` |
| `07_enrichment.R` | GSEA (fgsea) for every contrast of steps 04 and 05: MSigDB Hallmark, Reactome, GO:BP, plus a custom AD/Th2 set. | `enrichment/*` |
| `08_leip_reference.R` | For the proteins significant in step 06: LEIP normal range, where AD patients fall in it, clinical associations in LEIP, detectability, and LEIP vs in-study controls. | `leip_reference/*` (incl. `.xlsx`) |
| `09_isf_profile.R` | **Aim 1.** Descriptive dISF profile: detection class per protein (overall and by skin state), proteins detectable only in lesional skin, pathway over-representation of the detectable proteome, heatmap of the most variable proteins. NPX is protein-specific, so there is no ranking of levels between proteins. | `isf_profile/*` |
| `10_matrix_comparison.R` | **Aim 2.** Detected in dISF only / serum only / both. Relative dISF/serum enrichment, separately for AD lesional, ex-lesional and non-lesional skin and healthy skin: paired, centred log2 ratio, FDR plus ≥ 2-fold (`stats$min_rel_log2`). Concordance of disease effects in skin vs blood. Whether serum tracks the lesional-minus-non-lesional skin difference visit by visit. | `matrix_comparison/*` |
| `11_trajectories.R` | **Aim 3.** Residual lesional signal (ex-lesional minus non-lesional) vs weeks since clearance, and vs weeks to relapse (relapsers; exploratory). Serum vs weeks to relapse. Optional severity models, and per-patient trajectory plots of the top proteins (dISF lesional site, non-lesional site, serum; relapse marked). | `trajectories/*` |
| `12_focus_proteins.R` | **Dedicated analysis of pre-specified proteins** (`focus_proteins` in `config.yml`). Runs even if the protein fails the detection filter. Reports detection per matrix and group, and the same models as steps 04/05 plus the xL − NL relapse model for this protein alone (lmerTest / lm). The unadjusted p-value is the primary test, with the proteome-wide FDR shown alongside. Also ISF–serum correlation, severity (if available), LEIP clinical associations, every pipeline result for the protein, and figures. | `focus/<protein>/*`, `focus/focus_overview.*` |
| `13_visit_course.R` | dISF per visit (V1–V6, visits with ≥ `stats$visit_min_subjects` patients): tracked lesion site vs healthy skin, vs non-lesional skin, and non-lesional vs healthy. Same model as step 04, with a volcano plot per visit. "Regulated at all visits" = modelled and significant at every analysed visit in the same direction, at FDR < `stats$fdr` (strict) or p < 0.05 (nominal); with time-course plots of effects and NPX levels. | `visit_course/*`, `models/ISF_by_visit_*` |
| `14_serum_vs_disf.R` | The same question in dISF and serum, per visit and pooled: AD vs healthy (dISF lesion site or non-lesional skin vs healthy skin; MicroAD serum AD vs healthy) and relapse vs non-relapse (only visits before the relapse). Significant sets (FDR, and p < 0.05 as exploratory) are split into both (same / opposite direction), dISF only, dISF only because the protein isn't measurable in serum, and serum only. Output: Venn diagrams, stacked bars per visit, dISF volcano plots coloured by what serum shows, and dISF vs serum effect plots. | `serum_vs_disf/*` |
| `15_key_questions.R` | Answers Q1–Q5 with pre-specified single tests (see *Key questions* below), with one evidence figure per protein / score. | `key_questions/*` |
| `16_summary_report.R` | Executive summary PDF: key-question answers with evidence, data and QC, key findings per aim, all comparisons, visit course, volcano plots, pathways, dISF vs serum, serum, focus proteins, methods and caveats. Each section is replaced by a note if its step did not run. Plus a workbook with the full list behind every shortened list in the PDF. | `Executive_summary.pdf`, `Executive_summary_tables.xlsx` |
| `17_export_data.R` | CSV export of the Olink data: `samples.csv` (metadata + QC), `proteins.csv` (annotation, LOD, detection), wide tables per matrix (samples × proteins) as delivered (`NPX`) and as analysed (`PCNormalizedNPX`), a long table (`NPX_long.csv.gz`) with LOD and QC flags, and `export/RELAD2/` with the RELAD2 samples only. Values below LOD are kept as measured. For German Excel set `export: sep: ";"` in `config.yml`. | `export/*` |

### Models

The models are limma/dream (variancePartition) with empirical Bayes moderation. Subject is a random
effect wherever a person contributes several samples; models with a random effect but only one fixed
term (e.g. step 10) use limma with `duplicateCorrelation` instead, because dream rejects them in newer
variancePartition versions. NPX differences are on the log2 scale, and BH FDR is applied within each
contrast. Assays with > 20 % missing values in a model's samples are left out of that model; remaining
missing values are set to the protein median.

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

### Key questions (step 15)

- **Mast cell score:** each marker in `key_questions$mast_cell_markers` (KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, CPA3, CMA1, KIT, HDC; those on the panel) is z-standardised within its matrix; a sample's score is the mean z of its markers (at least 2 measured).
- **Q1** tests each marker and the score with the step 04/05 models: AD vs healthy (dISF lesional, non-lesional; serum), dISF lesional vs non-lesional, and relapse vs non-relapse (dISF xL and xL − NL, serum MicroAD, RELAD/RELAD2). Verdict from p < 0.05.
- **Q2:** Spearman correlation (across samples) and repeated-measures correlation (within patients) of TNFRSF9 / TNFSF9 with the score and each marker in AD dISF.
- **Q3:** all relapse tests plus the trend in the weeks before relapse (relapsers).
- **Q4:** *marker* = changes with lesion activity (L vs xL, L vs NL, serum active vs cleared); *predictor* = AUC (bootstrap 95 % CI) of values **before** the relapse, per patient. With < 10 patients per group a hit is reported only as "possible predictor".
- **Q5:** detectable proteins, step 14 overlap, key-protein effects and relapse AUCs in dISF vs serum.
- All are single pre-specified tests (unadjusted p-values); the answers are generated automatically and should be read with the evidence figures.

## Design decisions

- **`PCNormalizedNPX` is analysed, not intensity-normalised NPX.** Plate 1 is ISF only, plate 2 is mixed, and plates 3–4 are serum. Intensity normalisation assumes randomised samples of one matrix and would distort plate 2. Step 02 reports the `Normalization` column of the delivered file.
- **LOD** is computed per row with OlinkAnalyze's own routine (`olink_lod`), so count-based assays (`LODMethod = lod_count`, about 18 % in the fixed LOD file v10.2.0) get their sample-specific LOD.
  - Preferred: Olink's fixed LOD file (`paths$fixed_lod`), matched on `DataAnalysisRefID`.
  - Fallback for rows without a match: the Olink negative-control method on the 8 negative controls. OlinkAnalyze normally requires ≥ 10, so this LOD is less precise; in simulation it was off by up to ±0.6 NPX.
  - `qc/lod.csv` gives the LOD source and the per-sample LOD range for each assay.
- **Relapse in ISF is exploratory.** All four ISF relapsers are on plate 1, and plate 2 holds only non-relapsers, so relapse and plate cannot be fully separated.
- **Metadata issues are flagged, not fixed**, in `metadata/data_flags.csv`. As of manifest v4 only the 5 low-volume ISF samples remain flagged.
- **Sex** is available for MicroAD (manifest column `Sex`) and LEIP, **age** only for LEIP; the models are not adjusted for them.

## Testing without real data

```bash
Rscript tests/test_pipeline.R
```

This needs no study data. It builds a synthetic manifest with the same layout and design, simulates the two Olink NPX files (ISF, serum) with known effects, runs all
steps with `data_sim/config_sim.yml` (output in `output_sim/`), and checks that the effects are
recovered and false positives stay rare.
