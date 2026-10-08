# O-MicroAD – technical notes

R pipeline for the Olink Explore HT run of O-MicroAD: dermal interstitial fluid (dISF) and serum
from the MicroAD study (atopic dermatitis patients followed over up to 6 visits, healthy volunteers,
CPUO), plus serum from RELAD / RELAD2 and LEIP biobank controls. Many file and folder names (e.g.
`models/ISF_*`, `isf_profile/`, `isf_serum/`), the `matrix` column and `config.yml` use `ISF` for dISF;
some later outputs (`serum_vs_disf/*`, sheets of `Executive_summary_tables.xlsx`) write `dISF`. Output paths in this document are inside `output/`.

## Study aims and where they are answered

| Aim / topic | Steps | Key outputs |
|---|---|---|
| 1. Profile the dISF proteome of mild-to-moderate AD | 09, 04, 07 | `isf_profile/isf_profile.xlsx` (what is detectable, by skin state, pathways), `models/ISF_*` (lesional / ex-lesional / non-lesional / healthy) |
| 2. Compare the dISF and blood proteome | 10, 06, 08, 14 | `matrix_comparison/matrix_comparison.xlsx` (detected where, relatively enriched in dISF per skin type, skin vs blood disease signals), `isf_serum/*` (correlation), `leip_reference/*`, `serum_vs_disf/*` (overlap per visit) |
| 3. Track changes across the disease course | 11, 13, 04 (time models) | `trajectories/trajectory_results.csv`, per-patient plots in `trajectories/plots/`, `visit_course/*` (per visit, proteins regulated at all visits) |
| Focus proteins: CD137 (TNFRSF9), CD137L (TNFSF9), KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, PNOC, POSTN | 12 | `focus/focus_overview.xlsx`, `focus/focus_overview_heatmap.png`, one folder per protein |
| Key questions: mast cells, CD137 / CD137L, relapse, dISF vs serum | 15 | `key_questions/answers.csv`, `key_questions/key_questions.xlsx`, `key_questions/evidence/*.png` |
| Which proteins correlate with TNFRSF9 (CD137) in dISF (IL-33, IL-4 …) | 16 | `tnfrsf9_correlation/TNFRSF9_correlation.xlsx`, scatter plots |
| Serum vs dISF: same or different signatures (MicroAD only) | 17 | `signatures/signatures.xlsx` |
| All RELAD / RELAD2 serum results | 05 | `relad/RELAD_RELAD2_serum_results.xlsx` |
| **Executive summary of everything** | 18 | `Executive_summary.pdf`, `Executive_summary_tables.xlsx` (full lists) |
| **Data export: all proteins, all samples (CSV); RELAD2 separately** | 19 | `export/` |

Relapse analyses (04, 05, 11, 12, 14, 15, 16, 17) are exploratory in MicroAD; RELAD/RELAD2 serum gives larger groups (05).

dISF and serum are modelled separately. They are combined only where the two matrices are compared:
correlation (06), relative enrichment and concordance (10), focus-protein dISF–serum correlation (12),
overlap (14), Q5 (15), TNFRSF9 in dISF vs serum (16) and signatures (17, including models on paired
dISF − serum differences). Step 08 uses serum only, for the proteins found in step 06.

## Quick start

Open `oli_Jule.Rproj` in RStudio (or `setwd()` to this folder): all paths are relative to the project folder.

```bash
Rscript install_packages.R        # once
# put the data files in data/ (see below), then
Rscript run_all.R                 # runs all steps; results in output/, summary in output/Executive_summary.pdf
```

To run a single step, use e.g. `Rscript scripts/12_focus_proteins.R`. In RStudio, `start_at <- 12; source("run_all.R")`
resumes at step 12 and runs all later steps. Each step reads what the earlier ones saved.

## Data (never committed – `data/` and `output/` are git-ignored)

| File | What | Where it is set |
|---|---|---|
| `data/manifest.xlsx` | Olink sample submission sheet; the `manifest` sheet is the master (v4: incl. `Sex`, `NoRELAD2`) | `paths$manifest` |
| `data/npx/*.parquet` | Olink NPX files, here `O-MicroAD_ISF_NPX_2026-09-24.parquet` and `O-MicroAD_Serum_NPX_2026-09-24.parquet`; all files in the folder are read, and the file name must contain ISF or Serum | `paths$npx_dir` |
| `data/LEIP_clinical_parameters_n35.xlsx` | LEIP clinical data, sheet `Key_parameters` (optional) | `paths$leip_clinical` |
| `data/reference/proteinatlas.tsv.zip` | optional: Human Protein Atlas table (public) for the protein classes in step 02b and the tissue origin in step 17; download with `source("tools/download_hpa.R")` | `paths$hpa` |
| `data/severity.xlsx` | optional: `SubjectID`, `Visit` (V1–V6 or 1–6), plus numeric scores (e.g. SCORAD, EASI, itch NRS). Steps 11 and 12 then model the proteins against each score | `paths$severity` |
| `data/Explore_HT_Fixed_LOD.csv` | Olink fixed LOD file for Explore HT, version ≥ 6.0.0, from olink.com (recommended) | `paths$fixed_lod` |

If a file is not at the configured path, a single file with the usual name in the same folder is used instead
(e.g. `Explore HT_Fixed LOD.csv`, `*Sample Submission Sheet*.xlsx`). `run_all.R` first lists what it found.
All settings (paths, thresholds, FDR, focus proteins, mast cell markers, TNFRSF9 targets, enrichment and
signature collections, CSV separator, number of cores) are in `config.yml`.

## Steps

A plain-language explanation of every step (question, what it does, how to read it) is in [`STEPS.md`](STEPS.md).

`run_all.R` runs the steps in this order (see *Quick start* for resuming at a later step).

| Script | What it does | Main output |
|---|---|---|
| `01_metadata.R` | Cleans the manifest: harmonised groups (AD / HC / CPUO / Biobank), skin site and state, clinical state, sex, visit and days since V1, relapse visit and visits before relapse, sample volume, and LEIP clinical data. Checks IDs against the plate layout, and flags inconsistencies without dropping anything. | `metadata/sample_metadata.csv`, `metadata/data_flags.csv` |
| `02_import_qc.R` | Reads the parquet files, checks sample IDs, matrix and normalisation, and works out the LOD. Sample QC uses Olink flags plus a median/IQR outlier check. Removes assays Olink excluded (no values) and keeps an assay in a matrix if it is detected in ≥ 50 % of at least one group. Reports control CVs and low-volume samples. | `qc/*`; `data/npx_clean.rds`, `data/npx_wide.rds`, `data/npx_all_samples.rds` (the `data/` folder inside `output/`, not the input folder) |
| `02b_qc_overview.R` | Proteins above / below LOD per matrix: above LOD = ≥ `qc$min_detect_frac` of all samples of the matrix (or the step 02 filter, `qc_overview$rule`); serum restricted to `qc_overview$serum_cohorts` (default MicroAD). Protein classes from the HPA column "Protein class" (14 classes, each a union of HPA values, see `R/qc_overview.R`); assays matched per component by UniProt, then gene symbol, then synonym; classes overlap. Ring charts and stacked class bars (ggplot2 + grid), font from `qc_overview$font`: Nimbus Sans is shipped in `fonts/` (URW base35, AGPL-3 with font exception) and embedded via showtext (text drawn as outlines); step 02b stops if showtext is missing; other fonts: an installed font via cairo, else Helvetica. With showtext the PDF lists no font and the figure text is not selectable. | `qc_overview/*` |
| `03_explore.R` | Per matrix: PCA coloured by state, group, visit and plate (dISF) or by group, cohort, plate and relapse (serum); variance partitioning. | `explore/*` |
| `04_isf_models.R` | dISF models per protein (see *Models* below). | `models/ISF_*` |
| `05_serum_models.R` | Serum models per protein (see *Models* below). | `models/Serum_*`, `relad/RELAD_RELAD2_serum_results.xlsx` |
| `06_isf_vs_serum.R` | On matched MicroAD visits: within-subject (repeated-measures correlation, rmcorr) and between-subject (Spearman of per-person means) correlation of dISF and serum, separately for the lesion site (AD lesion site plus CPUO lesional skin) and the non-lesional site (AD and CPUO non-lesional skin plus healthy volunteers' skin). Only AD patients have several visits, so only they contribute to the within-subject correlation. BH FDR per site; `significant_proteins.csv` = FDR < `stats$fdr` within or between subjects at either site. | `isf_serum/*` |
| `07_enrichment.R` | GSEA (fgsea) for every contrast of steps 04 and 05: MSigDB Hallmark, Reactome, GO:BP, plus a custom AD/Th2 set. | `enrichment/*` |
| `08_leip_reference.R` | For the proteins significant in step 06: LEIP normal range, where AD patients fall in it, clinical associations in LEIP, detectability, and LEIP vs in-study controls. | `leip_reference/*` (incl. `.xlsx`) |
| `09_isf_profile.R` | **Aim 1.** Descriptive dISF profile: detection class per protein (overall and by skin state), proteins detectable only in lesional skin, pathway over-representation of the detectable proteome, clustered heatmaps (pheatmap, Ward on Euclidean distance of z-scores) of the 50 most variable proteins and of the top 25 up / 25 down lesional vs non-lesional proteins of step 04. NPX is protein-specific, so there is no ranking of levels between proteins. | `isf_profile/*` |
| `10_matrix_comparison.R` | **Aim 2.** MicroAD only. Detected (above LOD in ≥ `qc$min_detect_frac` of all MicroAD samples of the matrix) in dISF only / serum only / both. Relative dISF/serum enrichment for proteins detected in both, separately for AD lesional, ex-lesional and non-lesional skin and healthy skin: paired, centred log2 ratio, FDR plus ≥ 2-fold (`stats$min_rel_log2`). Concordance of disease effects: dISF L vs NL ↔ serum `MicroAD_active_vs_cleared`; dISF L vs HC and NL vs HC ↔ serum `MicroAD_AD_vs_HC`. Whether serum tracks the L − NL skin difference visit by visit. | `matrix_comparison/*` |
| `11_trajectories.R` | **Aim 3.** Residual lesional signal (ex-lesional minus non-lesional) vs weeks since clearance, and vs weeks to relapse (relapsers; exploratory). Serum vs weeks to relapse. Optional severity models, and per-patient trajectory plots of the top proteins (dISF lesion site, non-lesional site, serum; relapse marked). | `trajectories/*` |
| `12_focus_proteins.R` | **Dedicated analysis of pre-specified proteins** (`focus_proteins` in `config.yml`). Runs even if the protein fails the detection filter. Reports detection per matrix and group, and the same models as steps 04/05 plus the xL − NL relapse model for this protein alone (lmerTest / lm). The unadjusted p-value is the primary test, with the proteome-wide FDR shown alongside. Also dISF–serum correlation, severity (if available), LEIP clinical associations, every pipeline result for the protein, and figures. | `focus/<protein>/*`, `focus/focus_overview.*` |
| `13_visit_course.R` | dISF per visit (V1–V6, visits with ≥ `stats$visit_min_subjects` patients): tracked lesion site vs healthy skin, vs non-lesional skin, and non-lesional vs healthy. Same model as step 04, with a volcano plot per visit. "Regulated at all visits" = modelled and significant at every analysed visit in the same direction, at FDR < `stats$fdr` (strict) or p < 0.05 (nominal); with time-course plots of effects and NPX levels. | `visit_course/*`, `models/ISF_by_visit_*` |
| `14_serum_vs_disf.R` | The same question in dISF and serum, per visit and pooled: AD vs healthy (dISF lesion site or non-lesional skin vs healthy skin; MicroAD serum AD vs healthy) and relapse vs non-relapse (only visits before the relapse). **Under review:** the pooled (all visits) relapse comparison also includes the lesional V1 samples. Significant sets (FDR, and p < 0.05 as exploratory) are split into both (same / opposite direction), dISF only, dISF only because the protein isn't measurable in serum, and serum only. Output: Venn diagrams, stacked bars per visit, dISF volcano plots coloured by what serum shows, and dISF vs serum effect plots. | `serum_vs_disf/*` |
| `15_key_questions.R` | Answers Q1–Q5 with pre-specified tests and one evidence figure per protein / score (see *Key questions* below). The single-primary-test logic of Q2–Q5 is under review. | `key_questions/*` |
| `16_tnfrsf9_correlation.R` | TNFRSF9 vs the target proteins (config `tnfrsf9_correlation`) in dISF; see *TNFRSF9 correlation (step 16)* below. | `tnfrsf9_correlation/*` |
| `17_disf_serum_signatures.R` | MicroAD only: do serum and dISF carry the same or different signatures? See *Serum vs dISF signatures (step 17)* below. | `signatures/*` |
| `18_summary_report.R` | Executive summary PDF: overview (samples, data and QC), proteins above / below LOD by protein class (02b), key-question answers with evidence, TNFRSF9 correlations (16), serum vs dISF signatures (17), key findings per aim, all comparisons, visit course, dISF volcano plots, pathways, dISF vs serum, serum, serum vs dISF per visit (14), focus proteins, methods and caveats. Each section is replaced by a note if its step did not run. Plus a workbook with the full list behind every shortened list in the PDF (including the 02b class counts and the TNFRSF9 and signature lists). | `Executive_summary.pdf`, `Executive_summary_tables.xlsx` |
| `19_export_data.R` | CSV export of the Olink data: `samples.csv` (metadata + QC), `proteins.csv` (annotation, LOD, detection), wide tables per matrix (samples × proteins) as delivered (`NPX`) and as analysed (`PCNormalizedNPX`), a long table (`NPX_long.csv.gz`) with LOD and QC flags, and `export/RELAD2/` with the RELAD2 samples only. Values below LOD are kept as measured. For German Excel set `sep: ";"` under `export:` in `config.yml`. | `export/*` |

### Models

The models are limma/dream (variancePartition) with empirical Bayes moderation. Subject is a random
effect wherever a person contributes several samples; models with a random effect but only one fixed
term (e.g. step 10) use limma with `duplicateCorrelation` instead, because dream rejects them in newer
variancePartition versions. NPX differences are on the log2 scale, and BH FDR is applied within each
contrast. Assays with > 20 % missing values in a model's samples are left out of that model; remaining
missing values are set to the protein median. The group × visit F-tests of step 17 use limma +
`duplicateCorrelation`, because they test several coefficients jointly; plate is dropped from such a
model when it is confounded with visit.

**dISF (04)**, adjusted for plate (except the delta model):
- `states_all_visits`: AD lesional vs non-lesional, ex-lesional vs non-lesional, lesional vs ex-lesional; AD lesional and non-lesional vs healthy skin; AD lesional vs CPUO lesional; CPUO lesional vs non-lesional.
- `baseline_V1`: V1 only: AD lesional vs non-lesional, lesional vs healthy, non-lesional vs healthy; CPUO lesional vs non-lesional.
- `time_ex_lesional`, `time_non_lesional`: change per week since V1 in AD ex-lesional and non-lesional skin (drop-outs excluded from the non-lesional model).
- `relapse_ex_lesional`, exploratory: cleared skin of relapsers vs non-relapsers, adjusted for weeks since V1.
- `relapse_delta_xL_minus_NL`, exploratory: the same question on the within-visit ex-lesional minus non-lesional difference, adjusted for weeks since V1. It has no plate term: plate and systemic day-to-day variation cancel out within a pair.

**Serum (05):**
- AD vs controls, fitted twice: against the in-study healthy controls (adjusted for cohort and plate) and against LEIP biobank controls (adjusted for plate).
  - `Serum_AD_vs_controls_agreement.csv` marks proteins that agree in both comparisons.
  - `HC_vs_Biobank` shows proteins affected by the biobank source.
- `MicroAD_active_vs_cleared`: serum when the tracked lesion is active vs cleared.
- `MicroAD_relapse`: relapsers vs non-relapsers at cleared visits.
- `RELAD_relapse`: RELAD and RELAD2, adjusted for cohort and plate, plus a sensitivity analysis without the samples with conflicting relapse labels; also each cohort alone and `RELAD_AD_vs_HC`. All in `relad/RELAD_RELAD2_serum_results.xlsx`.
- `MicroAD_AD_vs_HC`: MicroAD AD (V1) vs healthy volunteers. Used for every dISF vs serum comparison (steps 10, 15), so that RELAD/RELAD2 serum never enters them; steps 06, 14, 16 and 17 use MicroAD serum only as well.

### Key questions (step 15)

- **Mast cell score:** each marker in `key_questions$mast_cell_markers` (KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, CPA3, CMA1, KIT, HDC; those on the panel) is z-standardised within its matrix; a sample's score is the mean z of its markers (at least 2 measured).
- **Q1** tests each marker and the score with the step 04/05 models: AD vs healthy (dISF lesional, non-lesional; serum), dISF lesional vs non-lesional, and relapse vs non-relapse (dISF xL and xL − NL, serum MicroAD, RELAD/RELAD2). Verdict from p < 0.05.
- **Q2:** Spearman correlation (across samples) and repeated-measures correlation (within patients) of TNFRSF9 / TNFSF9 with the score and each marker in AD dISF.
- **Q3:** all relapse tests plus the trend in the weeks before relapse (relapsers).
- **Q4:** *marker* = changes with lesion activity (L vs xL, L vs NL, serum active vs cleared); *predictor* = AUC (bootstrap 95 % CI) of values **before** the relapse, per patient. With < 10 patients per group a hit is reported only as "possible predictor".
- **Q5:** detectable proteins, step 14 overlap, key-protein effects and relapse AUCs in dISF vs serum.
- Verdicts are generated automatically and should be read with the evidence figures. Q1–Q4 use unadjusted p < 0.05 (a firm Q4 "predicts relapse" also needs ≥ 10 patients per group and a CI that excludes 0.5). Several answers combine more than one test: Q3 counts nominal hits across all relapse tests, and Q4 across up to five predictor sets. Q5 compares dISF with serum on detectable-protein counts, step 14 FDR overlap counts, nominal key-protein hits and |AUC − 0.5|.
- **Under review:** the wording of Q2–Q5 and the choice of one primary test per question.

### TNFRSF9 correlation (step 16)

dISF of AD patients and healthy volunteers; CPUO is not used.
- **Detectability first:** % of samples above LOD for TNFRSF9 and each target, per site and group. Values below LOD are used as measured; each targeted correlation is repeated on the pairs where both proteins are above LOD.
- **Strata:** lesion site (AD), non-lesional / healthy skin, and all dISF together.
- **Targeted:** Spearman, pooled, within groups (AD, healthy, relapse, non-relapse, lesional, ex-lesional) and per visit. *Pooled-only* flag: pooled p < 0.05, but no within-group subset smaller than the pooled set has p < 0.05 in the same direction (per-visit subsets are not counted). Such a correlation probably reflects group or skin-state differences, not a link between the two proteins.
- **Proteome-wide:** every protein that passes the dISF detection filter. Partial Spearman (rank-based residuals on skin state, visit and plate; used for the ranking) plus dream (`protein ~ TNFRSF9 + state + visit + plate + (1|subject)`). Ranked and BH-adjusted per stratum.
- **Within patients** (AD, per site): repeated-measures correlation (rmcorr) and delta-delta Spearman (visit-to-visit changes).
- **dISF vs serum** (MicroAD only): dISF vs serum TNFRSF9 at the same visit, and TNFRSF9 vs the targets in serum.

### Serum vs dISF signatures (step 17)

MicroAD samples only; RELAD, RELAD2 and LEIP serum are not used.
- **Measured in both** = passes both detection filters and is above LOD in ≥ `qc$min_detect_frac` of at least one MicroAD serum group (AD or healthy). dISF proteins not detected in MicroAD serum are reported separately.
- **Effect concordance** per comparison, per visit and pooled: Spearman, OLS slope of serum on dISF, % same sign. Comparisons: AD vs healthy (dISF lesion site or non-lesional skin vs healthy skin; serum AD vs healthy) and relapse vs non-relapse (dISF lesion site; serum).
- **Relapse vs non-relapse** uses only visits before the relapse: per visit (V1 is lesional for everyone; later visits only while the lesion is cleared), and pooled on cleared visits only (V1 excluded), adjusted for weeks since V1.
- **Group × visit F-tests** (limma + `duplicateCorrelation`). Healthy controls have one visit, so the test asks whether the AD − healthy difference changes over the visits; relapse × visit likewise within AD. dISF-significant proteins get a rule-based temporal profile class (resolving / persistent / late-rising / reversing / fluctuating), compared with the same proteins in serum (attenuation slope < 1 = weaker in serum).
- **Compartment × group:** AD vs healthy on pair-centred dISF − serum differences (as in step 10), i.e. proteins whose disease effect differs between dISF and serum.
- **Paired correlation** of dISF and serum levels, taken from step 06.
- **Signature sets** (dISF only, serum only, shared concordant, shared discordant; FDR and nominal) from the all-visits effects, with over-representation by `fgsea::fora` (background = assayed panel) and tissue origin from the Human Protein Atlas (`paths$hpa`).
- **Relapse, predictive:** the visit before the relapse vs non-relapsers' cleared visits, one mean per patient; for the dISF lesion site, non-lesional skin, ex-lesional minus non-lesional, and serum.

## Design decisions

- **`PCNormalizedNPX` is analysed, not intensity-normalised NPX.** Plate 1 is dISF only, plate 2 is mixed, and plates 3–4 are serum. Intensity normalisation assumes randomised samples of one matrix and would distort plate 2. Step 02 reports the `Normalization` column of the delivered file.
- **LOD** is computed per row with OlinkAnalyze's own routine (`olink_lod`), so count-based assays (`LODMethod = lod_count`, about 18 % in the fixed LOD file v10.2.0) get their sample-specific LOD.
  - Preferred: Olink's fixed LOD file (`paths$fixed_lod`), matched on `DataAnalysisRefID`.
  - Fallback for rows without a match: the Olink negative-control method on the 8 negative controls. OlinkAnalyze normally requires ≥ 10, so this LOD is less precise; in simulation it was off by up to ±0.6 NPX.
  - `qc/lod.csv` gives the LOD source and the per-sample LOD range for each assay.
- **Relapse in dISF is exploratory.** All four dISF relapsers are on plate 1, and plate 2 holds only non-relapsers, so relapse and plate cannot be fully separated.
- **Metadata issues are flagged, not fixed**, in `metadata/data_flags.csv`; the overview page of the executive summary counts them per issue.
- **Sex** is available for MicroAD (manifest column `Sex`) and LEIP, **age** only for LEIP; the models are not adjusted for them.

## Testing without real data

```bash
Rscript tests/test_pipeline.R
```

This needs no study data. It builds a synthetic manifest with the same layout and design, simulates
the two Olink NPX files (dISF, serum) with known effects, runs all steps with
`data_sim/config_sim.yml` (output in `output_sim/`), and checks that the effects are recovered and
false positives stay rare.
