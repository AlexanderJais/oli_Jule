# O-MicroAD – Olink proteomics of skin fluid and blood

This repository contains the R scripts that analyse the **Olink Explore HT** data of the O-MicroAD
study: about 5,400 proteins measured in **dermal interstitial fluid (dISF)** and in **serum** of
patients with atopic dermatitis (AD) and control persons.

This page tells you what the analysis does, how to run it, and **where to find which result**.
Technical details (models, design decisions) are in [`docs/TECHNICAL.md`](docs/TECHNICAL.md).

> **Start here:** after a run, open **`output/Executive_summary.pdf`**. It summarises all
> results on about 20 pages and tells you which folder holds the details.

---

## 1. The study in short

| | dISF (skin fluid) | Serum (blood) |
|---|---|---|
| **MicroAD** – AD patients followed over up to 6 visits (V1–V6) | 2 skin sites per visit: the tracked **lesion site** and a **non-lesional** site | 1 sample per visit |
| MicroAD – healthy volunteers (1 visit) | healthy skin | ✔ |
| MicroAD – CPUO patients (chronic pruritus, 1 visit) | lesional + non-lesional | ✔ |
| **RELAD / RELAD2** – AD patients with / without later relapse, healthy controls (1 sample each) | – | ✔ |
| **LEIP** – biobank serum, population reference (1 sample each) | – | ✔ |

**What happens at the lesion site over time:** at V1 the lesion is active (*lesional*). After
treatment it clears (*ex-lesional*). In some patients it comes back (**relapse**, lesional again);
in others it doesn't (*non-relapse*).

**Study aims**
1. Profile the dISF proteome of mild-to-moderate AD.
2. Compare the dISF proteome with blood (serum).
3. Follow the changes over the disease course (visits, clearing, relapse).

---

## 2. Words you will meet

| Word | Meaning |
|---|---|
| **NPX** | Olink's protein value, on a log2 scale: +1 NPX = twice as much protein. Only compare NPX **within one protein**; NPX values of different proteins are not comparable. |
| **PCNormalizedNPX** | NPX adjusted to the plate controls on each plate. This is the value the analysis uses. |
| **LOD** | Limit of detection. Values below it are mostly noise. |
| **Detected** | The protein is above LOD in ≥ 50 % of the samples of at least one group. Only detected proteins go into the statistics. |
| **L / lesion site** | The skin site followed over time: lesional at V1 and at relapse, ex-lesional in between |
| **NL** | Non-lesional (clinically normal-looking) skin of the same AD patient |
| **xL** | Ex-lesional: the lesion site after it has cleared |
| **HC** | Healthy controls |
| **logFC / estimate** | The difference between two groups in NPX (log2). +1 = twice as high in the first group. |
| **p-value** | How surprising the difference would be by chance, for **one** protein |
| **FDR (adj.P.Val)** | The p-value corrected for testing ~5,400 proteins at once. **Significant = FDR < 0.05.** With few samples (e.g. one visit) the FDR is strict, so few proteins pass. |
| **Volcano plot** | x axis = difference (logFC), y axis = −log10 p. Top left and top right = strongest changes. The dashed line is the FDR cutoff of that panel. |
| **Nominal / p < 0.05** | Significant without correction for many proteins. Only for exploring and generating hypotheses. |

---

## 3. How to run it

1. Open **`oli_Jule.Rproj`** in RStudio. This sets the working folder correctly.
2. Put the data files into the `data/` folder, as described in [`data/README.md`](data/README.md): manifest, the two Olink NPX files, the LEIP clinical data and the Olink LOD file. **Data files are never uploaded to GitHub.**
3. The first time only, run `source("install_packages.R")`.
4. Run `source("run_all.R")`. The full run takes about 15–45 minutes, and results appear in `output/`.
5. To restart from a later step, e.g. after changing something in step 12, run: `start_at <- 12; source("run_all.R")`.

The first lines of the run show whether all input files were found (`[ok]` / `[MISSING]`).

To check that R and the packages work, run `source("tests/test_pipeline.R")`. It uses invented
data only, takes about 5 minutes, and must end with "All pipeline checks passed".

---

## 4. What is where – the `output/` folder

Each script writes into its own folder. Most folders contain an **`.xlsx` file** that collects
all their tables. Open that one first.

### `Executive_summary.pdf` – the overview (step 16)
Key numbers and findings for each aim, all comparisons in one table, the most important figures,
and the methods and caveats. **Read this first.**

### `key_questions/` – answers to the key questions (step 15)
| Question | How it is answered | Where |
|---|---|---|
| Q1 Are mast cell markers elevated in AD, or only in relapse vs non-relapse? | Each marker (KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, CPA3, CMA1, KIT, HDC – those measured on the panel) and a combined **mast cell score**: AD vs healthy (dISF lesional, non-lesional, serum) and relapse vs non-relapse (dISF, serum MicroAD, RELAD/RELAD2) | `Q1_mast_cell_markers.png`, sheet `Q1_*` |
| Q2 Is CD137 (TNFRSF9) or CD137L (TNFSF9) a marker for mast cells in AD? | Correlation with the mast cell score and each marker in AD dISF, within patients over visits and across samples | `Q2_cd137_vs_mast_score.png`, sheet `Q2_correlations` |
| Q3 Do CD137 / CD137L correlate with relapse? | All relapse tests, plus the trend in the weeks before relapse | sheet `Q3_relapse` |
| Q4 Marker or predictor of relapse? | *Marker* = changes with lesion activity (lesional vs cleared). *Predictor* = values **before** the relapse separate relapsers from non-relapsers (AUC with 95 % CI) | `Q4_relapse_prediction_auc.png`, sheets `Q4_*` |
| Q5 Is dISF superior to serum? | Measurable proteins; significant proteins in dISF vs serum for the same question (step 14); key-protein effects; relapse AUC in dISF vs serum for the same patients | sheets `Q5_*` |

`answers.csv` holds one line per answer: the verdict and the numbers behind it. The same answers
are on the first pages of the executive summary. With 4 relapsing patients a relapse "hit" is
called **possible (exploratory)**; only groups of ≥ 10 (RELAD/RELAD2) can give a firm "predicts
relapse". Marker lists are set in `config.yml` under `key_questions`.

### `metadata/` – the samples (step 01)
| File | Content |
|---|---|
| `sample_metadata.csv` | One row per sample: patient, visit, skin site, skin state, group, relapse, plate, volume … |
| `data_flags.csv` | Inconsistencies found in the manifest, e.g. contradicting relapse labels or low sample volume. They are **flagged, not corrected**. |

### `qc/` – quality control (step 02)
| File | Content |
|---|---|
| `sample_qc.csv`, `sample_qc_median_iqr.png` | Olink QC per sample and outliers (flagged, not removed) |
| `assay_detection.csv`, `assay_detection.png` | For each protein and matrix: share of samples above LOD, and whether it is analysed (`keep`) |
| `lod.csv` | LOD per protein and where it comes from (Olink fixed LOD file) |
| `pca_by_plate.png` | Do the plates differ? The points should not separate by plate. |
| `sample_control_cv.csv` | Reproducibility of the sample controls (CV) |
| `assays_without_values.csv` | Proteins Olink excluded (no values) |
| `sample_id_mismatches.csv`, `matrix_mismatches.csv` | Should be empty. Otherwise sample IDs in the Olink files and the manifest don't match. |

### `explore/` – first look at the data (step 03)
PCA plots of each matrix, coloured by skin state, group, visit, plate or cohort. The
`variance_partition_*` files show how much of the variation comes from the person, the skin state,
the plate or the cohort.

### `models/` – the main group comparisons (steps 04, 05, 13)
| File | Content |
|---|---|
| `ISF_results.csv` | dISF: every protein × every comparison (see the list of comparison names below) |
| `ISF_by_visit_results.csv` | dISF comparisons **separately for each visit** V1–V6 |
| `ISF_relapse_delta_results.csv` | Relapse test on the difference ex-lesional minus non-lesional skin |
| `Serum_results.csv` | Serum: every protein × every comparison |
| `Serum_AD_vs_controls_agreement.csv` | Serum proteins that differ in AD **against both** control groups (in-study healthy and LEIP) |
| `*_summary.csv` | Number of significant proteins per comparison (quick overview) |
| `volcano/*.png` | One volcano plot per model, one panel per comparison |

**How to read a results table:** `Assay` = protein, `contrast` = comparison, `logFC` = difference
(log2), `P.Value` = p-value, `adj.P.Val` = FDR, `significant` = TRUE if FDR < 0.05. Filter for
`significant == TRUE` and sort by `P.Value`.

### `visit_course/` – visit by visit and the time course (step 13)
| File | Content |
|---|---|
| `significant_per_visit.png` | Number of significant proteins at each visit |
| `lesion_site_state_per_visit.csv` | How many lesion sites are lesional vs ex-lesional at each visit |
| `consistency_across_visits.csv` | For each protein: at how many visits it is significant, and in which direction |
| `time_course_effects_*.png` | Proteins regulated **at all visits**: difference vs healthy skin at each visit (with 95 % CI) |
| `time_course_heatmap_*.png` | The same as a heatmap |
| `time_course_levels.png` | NPX levels over the visits: lesion site, non-lesional skin, and healthy skin (grey band) |
| `visit_course.xlsx` | All of the above as tables |

Volcano plots for each visit: `models/volcano/ISF_by_visit_V1.png` … `V6.png`.

### `enrichment/` – pathways (step 07)
`gsea_results.csv`: for every comparison, which biological pathways (Hallmark, Reactome, GO,
Th2 set) are shifted. `NES > 0` = pathway higher in the first group of the comparison;
significant if `padj < 0.05`.

### `isf_profile/` – what is in skin fluid? (step 09, aim 1)
| File | Content |
|---|---|
| `isf_detection_profile.csv` | For each protein: how often it is detected in dISF, per skin state. `lesion_restricted` = detectable only in lesional AD skin. |
| `isf_detected_pathways.csv` | Which pathways the detectable dISF proteome covers |
| `top_variable_heatmap.png` | The 50 dISF proteins that vary most between samples |

### `isf_serum/` – do skin fluid and blood go together? (step 06)
`isf_serum_correlation.csv`: for each protein, the correlation of dISF with serum. `r_within` =
within the same patient over the visits; `r_between` = between patients.
`top_correlations.png` shows the strongest ones.

### `matrix_comparison/` – dISF vs serum proteome (step 10, aim 2)
| File | Content |
|---|---|
| `detection_by_matrix.csv` | Measurable in dISF only, serum only, both, or neither |
| `relative_enrichment.csv` | Proteins relatively **enriched in skin fluid** compared with blood, i.e. candidates for local production in the skin |
| `disease_signal_concordance.*` | Do disease differences seen in dISF also appear in serum? |

### `serum_vs_disf/` – overlap and what dISF adds, per visit (step 14)
The same question is asked in dISF and in serum: **AD vs healthy** and **relapse vs
non-relapse**, the latter using only visits before the relapse. It is done at each visit and for
all visits together, for the lesion site and the non-lesional site.

| File | Content |
|---|---|
| `venn_*_nominal.png`, `venn_*_FDR.png` | Venn diagrams: significant in dISF / in serum / in both |
| `overlap_bars_*.png` | The same as bars per visit. Red = information only dISF provides. |
| `volcano_dISF_coloured_*.png` | dISF volcano plots, coloured by whether serum shows the same |
| `effects_dISF_vs_serum_*.png` | Effect in dISF vs effect in serum for each protein |
| `overlap_summary.csv`, `overlap_protein_lists.csv`, `serum_vs_disf.xlsx` | Numbers and the protein names in each category |

`FDR` files use strict significance. `nominal` files use p < 0.05, which is exploratory but shows
more, since the groups per visit are small.

### `trajectories/` – disease course and relapse (step 11, aim 3)
Does the remaining lesional signal fade after clearing? Does it rise before a relapse? Results
are in `trajectory_results.csv`, and `plots/` holds one figure per protein showing each patient
over time (red = lesion site, blue = non-lesional, yellow = serum, dashed line = relapse). If a
severity file (SCORAD/EASI) is present, the proteins that follow the severity are listed too.

### `leip_reference/` – the population reference (step 08)
For the proteins where skin fluid and blood are correlated: the normal range in the LEIP
population, where AD patients fall in it, and whether the protein depends on age, sex, BMI, CRP,
lipids … in healthy people. `leip_reference.xlsx` collects all of it.

### `focus/` – the pre-specified proteins (step 12)
CD137 (TNFRSF9), CD137L (TNFSF9), KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, PNOC and POSTN.
The list is set in `config.yml` under `focus_proteins`.

| File | Content |
|---|---|
| `focus_overview_heatmap.png` | All focus proteins × main comparisons at a glance (stars = p-value) |
| `focus_skin_states.png` | All focus proteins in dISF by skin state |
| `focus_overview.xlsx` | The numbers behind it |
| `<protein>/<protein>_report.xlsx` | Everything about one protein: detection, all tests, correlation with serum, LEIP, sample values |
| `<protein>/*.png` | Skin states, V1 lesional vs non-lesional per patient, serum groups, course per patient, dISF vs serum |

For these proteins the **p-value** of the single test is the main result, because they were
chosen in advance. The FDR is shown alongside.

### `export/` – the data as CSV (step 17)
| File | Content |
|---|---|
| `ISF_NPX_wide.csv`, `Serum_NPX_wide.csv` | The Olink result: one row per sample, one column per protein (NPX as delivered) |
| `*_PCNormalizedNPX_wide.csv` | The same with the values used in the analysis |
| `samples.csv`, `proteins.csv` | Sample information and protein information (UniProt, LOD, detection) |
| `NPX_long.csv.gz` | Everything in one long table, including LOD and QC flags |
| `RELAD2/` | The RELAD2 samples only |

Tip: if Excel shows everything in one column, set `sep: ";"` under `export:` in `config.yml`
and rerun step 17 (`start_at <- 17; source("run_all.R")`).

### `leip_biobank/` – LEIP biobank only: proteins vs clinical data, and galanin (step 18)
Only the LEIP biobank sera (population controls) are used. Every protein is correlated with every
clinical parameter of the LEIP clinical file (sheet `Key_parameters` plus all further SORB
variables in sheet `All_SORB_parameters`), and galanin gets its own analysis.
**Start with `LEIP_biobank_summary.pdf`**: the answers are on its first pages.

| Question | How it is answered | Where |
|---|---|---|
| G1 Is galanin (Olink assay GAL) measurable in LEIP serum? | Share of LEIP samples above LOD | `answers.csv` |
| G2 Do the Olink GAL levels correspond to the galanin ELISA? | Spearman correlation (also adjusted for Olink plate, within each plate, above LOD only, without flagged samples); same tertile in both methods; GAL's rank among all proteins correlated with the ELISA. **Benchmark:** the same comparison for every other protein measured both by the lab and by Olink (e.g. A-FABP/FABP4, progranulin/GRN, chemerin/RARRES2) | `galanin/GAL_Olink_vs_ELISA.png`, `galanin/lab_vs_Olink_benchmark.png` |
| G3 Which clinical parameters does galanin correlate with? | Olink GAL and the ELISA vs each parameter, unadjusted and adjusted for age and sex (Olink GAL also for plate). Both methods side by side | `galanin/galanin_vs_clinical.png`, `galanin/GAL_top_clinical_scatter.png` |
| G4 Which proteins does galanin correlate with? | Olink GAL (and the ELISA) vs every other protein; gene-set enrichment of the ranking | `galanin/GAL_vs_proteins_volcano.png`, `galanin/GAL_top_proteins_scatter.png`, `galanin/GAL_gene_sets.png` |
| P1 Which clinical parameters show up in the serum proteome? | Every protein × every parameter | `n_significant_per_parameter.png`, `heatmap_key_parameters.png`, `volcano_key_parameters.png`, `top_associations.png` |
| P2 Sanity check: do known associations show up? | E.g. leptin – BMI, cystatin C – eGFR, GDF15 – age (list in `config.yml`) | `sanity_checks.csv` |

| File | Content |
|---|---|
| `LEIP_biobank_summary.pdf` | Answers, the key figures and tables, methods and caveats |
| `answers.csv` | One line per answer: the verdict and the numbers behind it |
| `leip_biobank.xlsx` | Which clinical variables are analysed (and why the others are not), samples, detection, all significant associations, top 10 proteins per parameter, sanity checks |
| `associations_all.csv.gz` | Every protein × every parameter: `rho`, `p`, `fdr`; `rho_adj`, `p_adj`, `fdr_adj` = adjusted for age, sex and Olink plate |
| `n_significant_per_parameter.csv` | Per parameter: number of significant proteins, and how many would be expected by chance at p < 0.05 |
| `sample_check.csv` | Manifest vs clinical file: same person and Olink plate for every sample? |
| `galanin/galanin.xlsx` | All galanin results: agreement with the ELISA, tertiles, plates, benchmark, clinical parameters, proteins, gene sets, sample values |

**How to read it:** `rho` runs from −1 to +1; +0.5 means the protein tends to be high when the
parameter is high. With 34 persons, only |rho| ≳ 0.35 reaches p < 0.05, and only |rho| ≳ 0.6
passes the proteome-wide FDR, so weaker true correlations are missed. Olink and ELISA can only
agree in *ranking* (NPX is relative), not in absolute values. Clinical column names are used in
plain ASCII (e.g. `Ins0_µU_ml` becomes `Ins0_uU_ml`). Settings: `leip_biobank` in `config.yml`.
Step 18 only needs steps 01–02: `start_at <- 18; source("run_all.R")`.

### `data/` (inside `output/`)
Intermediate files used by the scripts (`.rds`). You don't need to open them.

---

## 5. Comparison names used in the tables

| Name | Meaning (first group vs second group) |
|---|---|
| `AD_L_vs_NL` | AD lesional skin vs non-lesional skin of the same patients |
| `AD_xL_vs_NL` | Ex-lesional (cleared) skin vs non-lesional skin: what remains after clearing |
| `AD_L_vs_xL` | Lesional vs ex-lesional: what goes away with clearing |
| `AD_L_vs_HC`, `AD_NL_vs_HC` | AD lesional / non-lesional skin vs healthy skin |
| `CPUO_L_vs_NL`, `AD_L_vs_CPUO_L` | Pruritus (CPUO) lesional vs non-lesional; AD lesional vs CPUO lesional |
| `per_week_xL`, `per_week_NL` | Change per week in ex-lesional / non-lesional skin |
| `relapse_vs_non` | Patients who relapse vs who don't |
| `Lsite_vs_HC`, `Lsite_vs_NL`, `NL_vs_HC` (per visit) | At one visit: lesion site vs healthy / vs non-lesional; non-lesional vs healthy |
| `AD_vs_HC` (serum) | AD patients vs healthy controls of the same studies |
| `AD_vs_Biobank`, `HC_vs_Biobank` | vs LEIP biobank serum (the latter checks for biobank handling effects) |
| `active_vs_cleared` (serum) | Serum at visits with an active lesion vs after clearing |

Model names: `states_all_visits` = all visits together; `baseline_V1` = V1 only;
`V1`…`V6` = one visit; `*_relapse*` = relapse analyses.

---

## 6. Please keep in mind

- **Relapse results are exploratory.** There are only 4 relapsing patients in MicroAD.
- **Few samples per visit** (6–11 patients) make the FDR strict. "Not significant" does **not** mean "no difference".
- **NPX is relative:** compare a protein between groups, never protein A with protein B.
- **LEIP biobank serum** was collected and stored differently. Trust AD-vs-LEIP differences only if they also appear against the in-study healthy controls.
- **Metadata issues** (see `metadata/data_flags.csv`) are flagged, not fixed.

---

## 7. What is in the repository

| Folder / file | Content |
|---|---|
| `run_all.R` | Runs everything in order |
| `config.yml` | All settings: file paths, thresholds, focus proteins, CSV separator |
| `scripts/01_…` to `scripts/18_…` | One script per analysis step. The number is also the step in `start_at`. |
| `R/` | Shared functions used by the scripts |
| `data/` | Your input files (not uploaded to GitHub) – see `data/README.md` |
| `output/` | All results (not uploaded to GitHub) |
| `tests/` | Test with invented data |
| `docs/TECHNICAL.md` | Statistical methods and design decisions |

| Step | Script | What it does |
|---|---|---|
| 01 | `01_metadata.R` | Read and check the manifest |
| 02 | `02_import_qc.R` | Read the Olink files, LOD, quality control |
| 03 | `03_explore.R` | PCA, sources of variation |
| 04 | `04_isf_models.R` | dISF comparisons |
| 05 | `05_serum_models.R` | Serum comparisons |
| 06 | `06_isf_vs_serum.R` | dISF–serum correlation |
| 07 | `07_enrichment.R` | Pathways |
| 08 | `08_leip_reference.R` | LEIP population reference |
| 09 | `09_isf_profile.R` | dISF proteome profile (aim 1) |
| 10 | `10_matrix_comparison.R` | dISF vs serum proteome (aim 2) |
| 11 | `11_trajectories.R` | Disease course (aim 3) |
| 12 | `12_focus_proteins.R` | Focus proteins (CD137 …) |
| 13 | `13_visit_course.R` | Per visit and time course |
| 14 | `14_serum_vs_disf.R` | Serum vs dISF overlap per visit |
| 15 | `15_key_questions.R` | Answers to the key questions (mast cells, CD137, relapse, dISF vs serum) |
| 16 | `16_summary_report.R` | Executive summary PDF |
| 17 | `17_export_data.R` | CSV export of the data |
| 18 | `18_leip_biobank.R` | LEIP biobank only: proteins vs clinical parameters, galanin (Olink vs ELISA) |

---

## 8. Something went wrong?

| Message | Fix |
|---|---|
| `cannot open file 'R/utils.R'` / "Working directory must be the project folder" | Open `oli_Jule.Rproj`, or use `setwd()` to go to the repository folder |
| `[MISSING] manifest` / `NPX parquet files` | Check the file names and places in `data/README.md` |
| `LOD source: negative controls` | The Olink fixed LOD file was not found. Put it in `data/`; its name must contain "Fixed LOD". |
| A step fails | Read the last lines of the message, fix the problem, and restart from that step with `start_at <- N; source("run_all.R")` |
| Package error | Run `source("install_packages.R")` again |

When asking for help, send the error message and the names of files, columns or samples,
**not the data itself**.
