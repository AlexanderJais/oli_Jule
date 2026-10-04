# O-MicroAD – Olink proteomics of skin fluid and blood

This repository contains the R scripts that analyse the **Olink Explore HT** data of the O-MicroAD
study: about 5,400 proteins measured in **dermal interstitial fluid (dISF)** and in **serum** of
patients with atopic dermatitis (AD) and control persons.

This page tells you how to run the analysis and **where to find which result**. Two more pages:
- [`docs/STEPS.md`](docs/STEPS.md) – what each of the 19 steps does, in plain words;
- [`docs/TECHNICAL.md`](docs/TECHNICAL.md) – technical details (models, design decisions).

> **Start here:** after a run, open **`output/Executive_summary.pdf`**. It summarises all
> results (the answers to the key questions come right after the overview page) and tells you
> which folder holds the details.

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

**Study aims** (→ where the answers are in `output/`)
1. Profile the dISF proteome of mild-to-moderate AD → `isf_profile/`, `models/ISF_*`.
2. Compare the dISF proteome with blood (serum) → `matrix_comparison/`, `isf_serum/`, `serum_vs_disf/`, `signatures/`, `leip_reference/`.
3. Follow the changes over the disease course (visits, clearing, relapse) → `trajectories/`, `visit_course/`.

---

## 2. Words you will meet

| Word | Meaning |
|---|---|
| **dISF / ISF** | Dermal interstitial fluid (skin fluid). Some file names (e.g. `models/ISF_results.csv`, `explore/pca_ISF_*.png`), the `matrix` column and `config.yml` use the short form `ISF`; other files and sheets say `dISF`. |
| **NPX** | Olink's protein value, on a log2 scale: +1 NPX = twice as much protein. Only compare NPX **within one protein**; NPX values of different proteins are not comparable. |
| **PCNormalizedNPX** | NPX adjusted to the plate controls on each plate. This is the value the analysis uses. |
| **LOD** | Limit of detection. Values below it are mostly noise. |
| **Detected** | The protein is above LOD in ≥ 50 % of the samples of at least one group. Only detected proteins go into the statistics. |
| **L** | Lesional skin (active eczema), e.g. in `AD_L_vs_NL` |
| **Lesion site** (`Lsite`) | The skin site followed over time: lesional at V1 and at relapse, ex-lesional in between. In `sample_metadata.csv` it is `site = L`. |
| **NL** | Non-lesional (clinically normal-looking) skin of the same AD patient |
| **xL** | Ex-lesional: the lesion site after it has cleared |
| **HC** | Healthy controls |
| **logFC / estimate** | The difference between two groups in NPX (log2). +1 = twice as high in the first group. |
| **p-value** | How surprising the difference would be by chance, for **one** protein |
| **FDR (adj.P.Val)** | The p-value corrected for testing all proteins of that comparison at once (all detected proteins: up to ~5,400, fewer in dISF). **Significant = FDR < 0.05.** With few samples (e.g. one visit) the FDR is strict, so few proteins pass. |
| **Volcano plot** | x axis = difference (logFC), y axis = −log10 p. Top left and top right = strongest changes. The dashed line is the FDR cutoff of that panel. |
| **Nominal / p < 0.05** | Significant without correction for many proteins. Only for exploring and generating hypotheses. |
| **PCA** | A plot that places similar samples close together. |
| **Correlation (Spearman ρ)** | How closely two values rise and fall together: −1 = opposite, 0 = unrelated, +1 = together. |
| **Partial correlation** | The correlation left after removing differences due to skin state / group, visit and plate, so two proteins that are merely both higher in lesional skin don't count as linked. |
| **Within-patient correlation** | Do two values go up and down together from visit to visit in the same patient? Differences between patients don't count. Also called repeated-measures correlation (rmcorr). |
| **AUC** | How well a value separates two groups: 0.5 = no better than chance; the further from 0.5 (towards 1 or towards 0), the better. In `key_questions/`, above 0.5 = higher in the patients who later relapse, below 0.5 = lower. |

---

## 3. How to run it

1. Open **`oli_Jule.Rproj`** in RStudio. This sets the working folder correctly.
2. Put the data files into the `data/` folder, as described in [`data/README.md`](data/README.md).
   **Required:** the manifest, and the two Olink NPX files (`.parquet`) in the subfolder `data/npx/`.
   **Recommended:** the Olink fixed LOD file. **Optional:** the LEIP clinical data and the severity scores.
   **Data files are never uploaded to GitHub.**
3. The first time only, run `source("install_packages.R")` (needs R 4.3 or newer).
4. Optional, once (needs internet): `source("tools/download_hpa.R")` downloads the public Human Protein
   Atlas table that step 17 uses for the tissue origin of proteins. Without it, step 17 still runs but
   leaves the tissue origin out.
5. Run `source("run_all.R")`. The full run can take an hour or more, and the results appear in `output/`.
6. **To resume** from a later step, e.g. after changing a setting for step 12: `start_at <- 12; source("run_all.R")`.
   Earlier results are reused, and all later steps run again (including the summary, step 18). Start
   from the **earliest** step you changed (step numbers are in [section 7](#7-what-is-in-the-repository)).
   `start_at` stays set until you restart R: before the next full run, type `rm(start_at)`.

The first lines of the run list the input files: `[ok]` = found, `[MISSING]` = required file
missing (the run stops), `[--]` = optional or recommended file not found (the run continues).

To check that R and the packages work, run `source("tests/test_pipeline.R")`. It uses invented data
only, takes about 15 minutes, and must end with "All pipeline checks passed". Its results go to
`output_sim/` and are **invented**: never use them as study results. **Afterwards restart R**
(RStudio: Session → Restart R) before you run the real analysis; otherwise `run_all.R` keeps using
the test settings.

---

## 4. What is where – the `output/` folder

Each script writes into its own folder. Most folders contain an **`.xlsx` file** that collects
all their tables. Open that one first.

### `Executive_summary.pdf` – the overview (step 18)
First an overview of the samples and the quality control, then the answers to the key questions
with the data behind them, the TNFRSF9 correlations (step 16) and the serum vs dISF signatures
(step 17), key numbers and findings for each aim, all comparisons in one table, the most important
figures, and the methods and caveats. **Read this first.**

### `Executive_summary_tables.xlsx` – the full lists behind the summary (step 18)
Wherever the PDF shows a shortened list ("up: TNC, LAIR2 …"), the complete list is in this
workbook. The first sheet, `index`, says what each sheet contains:
- proteins detectable only in lesional skin;
- all significant proteins for lesional / ex-lesional / non-lesional / healthy comparisons;
- proteins regulated at every visit, with per-visit values;
- all proteins per visit for lesion site vs non-lesional skin and vs healthy skin;
- dISF-vs-serum enrichment for each skin type;
- serum lists, key-question answers and focus proteins;
- TNFRSF9 (step 16): all targeted correlations, and all dISF proteins ranked by their correlation with TNFRSF9;
- serum vs dISF (step 17): the signature sets (FDR < 0.05) and the proteins with p < 0.05 at the visit before relapse.

### `key_questions/` – answers to the key questions (step 15)
All tables are in `key_questions.xlsx`: start with the sheet `answers`. The sheets named below are in this workbook.

| Question | How it is answered | Where |
|---|---|---|
| Q1 Are mast cell markers elevated in AD, or only in relapse vs non-relapse? | Each mast cell marker and a combined **mast cell score**: AD vs healthy (dISF lesional, non-lesional, serum), dISF lesional vs non-lesional, and relapse vs non-relapse (dISF, serum MicroAD, RELAD/RELAD2) | `Q1_mast_cell_markers.png`, sheets `Q1_*` |
| Q2 Is CD137 (TNFRSF9) or CD137L (TNFSF9) a marker for mast cells in AD? | Correlation with the mast cell score and each marker in AD dISF, within patients over visits and across samples | `Q2_cd137_vs_mast_score.png`, sheet `Q2_correlations` |
| Q3 Do CD137 / CD137L correlate with relapse? | All relapse tests, plus the trend in the weeks before relapse | sheet `Q3_relapse` |
| Q4 Marker or predictor of relapse? | *Marker* = changes with lesion activity (lesional vs cleared). *Predictor* = values **before** the relapse separate relapsers from non-relapsers (AUC with 95 % CI) | `Q4_relapse_prediction_auc.png`, sheets `Q4_*` |
| Q5 Is dISF superior to serum? | Measurable proteins; significant proteins in dISF vs serum for the same question (step 14); key-protein effects; relapse AUC in dISF vs serum for the same patients | sheets `Q5_*` |

`evidence/<protein>.png` shows, for each mast cell marker, the mast cell score, CD137 and CD137L,
the data behind the statement: dISF by skin state, serum AD vs controls, and values before relapse
in relapsers vs non-relapsers. The test results are in the subtitle. The **mast cell score** is the
mean of the z-standardised markers measured in a sample (KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1,
CPA3, CMA1, KIT, HDC – those on the panel). A single marker can be significant while the score is
not, if the other markers don't move with it.

`answers.csv` holds one line per answer: the verdict and the numbers behind it. The same answers
open the executive summary, right after the overview page. With only 4 relapsing patients in MicroAD a
relapse "hit" is called **possible predictor (exploratory)**; only groups of ≥ 10 (RELAD/RELAD2)
can give a firm "predicts relapse". Marker lists are set in `config.yml` under `key_questions`.

The verdicts are generated automatically from p < 0.05: read them together with the evidence
figures. **Under review:** the wording of Q2–Q5 and the choice of one primary test per question.

### `tnfrsf9_correlation/` – which proteins correlate with TNFRSF9 / CD137 in dISF? (step 16)
TNFRSF9 is compared with IL33 (IL-33), IL4 (IL-4), CSF2, IL6, IL18, CXCL8, IL1RL1, KIT, KITLG, TPSAB1,
TPSB2 and FCER1A (the list is in `config.yml` under `tnfrsf9_correlation`; the tables use these
Olink names). TPSB2 is not on the Explore HT panel.

| File | Content |
|---|---|
| `TNFRSF9_correlation.xlsx` | Everything below in one workbook. Start with `README`, then `LOD_summary`. |
| `LOD_summary` (sheet) | % of dISF samples above LOD per site and group. **Read first:** proteins mostly below LOD give unreliable correlations. |
| `targeted_correlations.csv` | Spearman ρ, p and n for TNFRSF9 vs each target: pooled, within each group (AD, healthy, relapse, non-relapse, lesional, ex-lesional) and per visit; for the lesion site, for non-lesional / healthy skin and for all dISF together (column `site`). The `*_both_above_LOD` columns repeat each correlation with only the samples where both proteins are above LOD. `interpretation` = "pooled only" when the pooled correlation has p < 0.05 but no within-group one does: it then probably reflects group differences, not a link between the two proteins. |
| `proteome_wide_correlation.csv` | Every dISF protein that passes the detection filter, ranked by its partial correlation with TNFRSF9 (adjusted for skin state / group, visit and plate), with FDR; a mixed model is in the `*_mixed` columns. Ranked separately for all dISF, the lesion site and non-lesional / healthy skin (column `site`). `target = TRUE` marks IL33, IL4 … |
| `longitudinal_within_patient.csv` | Do TNFRSF9 and the targets change together from visit to visit in the same patient? |
| `cross_compartment.csv`, `serum_correlations.csv` | dISF TNFRSF9 vs serum TNFRSF9 at the same visit; TNFRSF9 vs targets in serum (MicroAD only) |
| `scatter/*.png`, `targeted_heatmap.png`, `proteome_wide_correlation.png` | Scatter plots per target (coloured by group, per site; open symbols = below LOD) and overviews |

### `signatures/` – do serum and dISF carry the same or different signatures? (step 17)
**MicroAD only** – RELAD, RELAD2 and LEIP serum are not used. Open `signatures.xlsx` (sheet `README` explains each sheet).

| Analysis | Where |
|---|---|
| 1 Effect-size concordance: serum vs dISF log2FC per comparison and visit (ρ, slope, % same sign); dISF proteins not detected in serum listed separately | `effect_concordance*.png`, sheet `concordance`; sheets `not_detected_in_serum*` only if some dISF proteins are not detected in serum |
| 2 Group × visit models per compartment; temporal profiles (resolving / persistent / late-rising …) and whether serum shows a weaker version | `temporal_profiles.png`, sheets `time_models_*`, `temporal_profiles*` |
| 3 Compartment × group: proteins whose disease effect differs between dISF and serum | sheet `compartment_x_group` |
| 4 Within-patient paired correlation: systemic spill-over vs local production candidates | sheet `paired_correlation` |
| 5 Signature sets (dISF-only, serum-only, shared-concordant, shared-discordant) with Reactome/GO enrichment and tissue origin | sheets `signature_sets`, `signature_set_counts`, `enrichment`; `signature_origin_counts` only with the Human Protein Atlas file |
| 6 Relapse, predictive: the visit before the relapse vs non-relapsers | `relapse_predictive_volcano.png`, sheets `relapse_predictive*` |

**How to read:** start with `concordance`. ρ near 1 and a high % same sign = serum shows the same
biology as dISF; a slope between 0 and 1 = the same signal, but weaker in serum. Sheets that would
be empty are left out of the workbook.

Tissue origin uses the Human Protein Atlas table. Download it once with `source("tools/download_hpa.R")`
(public data, saved in `data/reference/`), then run `start_at <- 17; source("run_all.R")`, which
also updates the summary. Without the file, the tissue origin is left out.

### `relad/` – all RELAD and RELAD2 serum results (step 05)
`RELAD_RELAD2_serum_results.xlsx`: relapse vs non-relapse (both cohorts together, RELAD alone,
RELAD2 alone, without conflicting labels), AD vs healthy within RELAD/RELAD2, group means with % above
LOD per protein, all sample labels and the NPX values. Start with the sheets `README` and `summary`;
"without conflicting labels" is the sheet `RELAD_relapse_unflagged`.

RELAD and RELAD2 are **not** used when dISF is compared with serum: steps 06, 10, 14, 16 and 17 use
MicroAD serum only, and so do the Q5 comparisons of step 15 (the Q5 count of measurable serum proteins
uses the step 02 detection filter over all serum samples). (Q1, Q3 and Q4 of step 15 also report
RELAD/RELAD2 serum, as separate tests.)

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
| `lod.csv` | LOD per protein and where it comes from (Olink fixed LOD file, or the negative controls if that file is missing) |
| `pca_by_plate.png` | Do the plates differ? The points should not separate by plate. |
| `sample_control_cv.csv` | Reproducibility of the sample controls (CV) |
| `assays_without_values.csv` | Proteins Olink excluded (no values) |
| `sample_id_mismatches.csv`, `matrix_mismatches.csv` | Should be empty. Otherwise sample IDs in the Olink files and the manifest don't match. |

### `explore/` – first look at the data (step 03)
PCA plots of each matrix (`pca_<matrix>_<colour>.png`): dISF coloured by skin state, group, visit
and plate; serum by group, cohort, plate and relapse. The `variance_partition_*` files show how much
of the variation comes from the person, the skin state and the plate (dISF), or from the group, the
cohort and the plate (serum).

### `models/` – the main group comparisons (steps 04, 05, 13)
| File | Content |
|---|---|
| `ISF_results.csv` | dISF: every protein × every comparison (see the list of comparison names below) |
| `ISF_by_visit_results.csv` | dISF comparisons **separately for each visit** V1–V6 |
| `ISF_relapse_delta_results.csv` | Relapse test on the difference ex-lesional minus non-lesional skin |
| `Serum_results.csv` | Serum: every protein × every comparison |
| `Serum_AD_vs_controls_agreement.csv` | All serum proteins with AD vs in-study healthy and AD vs LEIP side by side. `agree_both_controls = TRUE` marks proteins that differ in AD **against both** control groups, in the same direction |
| `*_summary.csv` | Number of significant proteins per comparison (quick overview) |
| `volcano/*.png` | One volcano plot per model, one panel per comparison |

**How to read a results table:** `Assay` = protein, `model` = which analysis, `contrast` =
comparison, `logFC` = difference (log2), `P.Value` = p-value, `adj.P.Val` = FDR, `significant` =
TRUE if FDR < 0.05. The same contrast name can occur in several models (e.g. `AD_vs_HC`,
`relapse_vs_non`), so filter on `model` **and** `contrast` first, then on `significant == TRUE`,
and sort by `P.Value`. Model and comparison names are explained in [section 5](#5-comparison-names-used-in-the-tables).

### `visit_course/` – visit by visit and the time course (step 13)
| File | Content |
|---|---|
| `significant_per_visit.png` | Number of significant proteins at each visit |
| `lesion_site_state_per_visit.csv` | How many lesion sites are lesional vs ex-lesional at each visit |
| `consistency_across_visits.csv` | For each protein: at how many visits it is significant, and in which direction. "Regulated at all visits" = significant at **every** analysed visit, same direction |
| `time_course_effects_*.png` | Proteins regulated **at all visits**: the difference (vs healthy or vs non-lesional skin, see file name) at each visit, with 95 % CI |
| `time_course_heatmap_*.png` | The same as a heatmap |
| `time_course_levels.png` | NPX levels over the visits: lesion site, non-lesional skin, and healthy skin (grey band) |
| `visit_course.xlsx` | All of the above as tables |

Volcano plots for each visit: `models/volcano/ISF_by_visit_V1.png` … `V6.png`.

### `enrichment/` – pathways (step 07)
`gsea_results.csv`: for every comparison of steps 04 and 05, which biological pathways (Hallmark, Reactome, GO,
Th2 set) are shifted. `NES > 0` = pathway higher in the first group of the comparison;
significant if `padj < 0.05`.

### `isf_profile/` – what is in skin fluid? (step 09, aim 1)
| File | Content |
|---|---|
| `isf_profile.xlsx` | Both tables below in one workbook (sheets `detection_profile`, `detected_pathways`) |
| `isf_detection_profile.csv` | For each protein: how often it is detected in dISF, per skin state. `lesion_restricted` = detectable only in lesional AD skin. |
| `isf_detected_pathways.csv` | Which pathways the detectable dISF proteome covers |
| `detected_per_sample.png` | Number of proteins above LOD in each dISF sample, by skin state |
| `top_variable_heatmap.png` | The 50 dISF proteins that vary most between samples, chosen **without** using any group information. Rows and samples are clustered (similar ones side by side); colour bars show skin state, patient, visit and plate. Use it to see which proteins move together and whether samples group by skin state or by patient. |
| `lesional_vs_nonlesional_heatmap.png` | The disease signal: the 25 proteins most clearly higher and the 25 most clearly lower in lesional than in non-lesional skin (step 04, all visits), shown in lesional, ex-lesional, non-lesional and healthy skin. Rows are clustered; the black bar marks FDR-significant proteins. |
| `isf_profile.xlsx` (heatmap sheets) | The numbers behind both heatmaps, in the order shown: `variable_*` and `lesional_*` sheets with z-scores, NPX values and the sample order (`*_samples`). |

### `isf_serum/` – do skin fluid and blood go together? (step 06)
`isf_serum_correlation.csv`: for each protein and skin site, the correlation of dISF with serum taken
at the same visit (MicroAD only). Column `site`: L = lesion site, NL = non-lesional skin (the healthy
volunteers' skin counts as NL). `r_within` = within the same person over the visits (in practice
the AD patients, the only ones with several visits); `r_between` = between people, which can partly
reflect AD vs healthy differences. `significant_proteins.csv` lists the proteins with FDR < 0.05
(used by step 08), `matched_pairs.csv` the number of dISF–serum pairs, and `top_correlations.png`
shows the strongest correlations.

### `matrix_comparison/` – dISF vs serum proteome (step 10, aim 2)
MicroAD samples only.

| File | Content |
|---|---|
| `matrix_comparison.xlsx` | All tables of this step (sheets `detection`, `relative_enrichment`, `disease_concordance`, `lesion_signal_vs_serum`) |
| `detection_by_matrix.csv` | Measurable (≥ 50 % of samples above LOD) in dISF only, serum only, both, or neither |
| `relative_enrichment.csv` | Proteins relatively **enriched in skin fluid** compared with blood (candidates for local production in the skin), separately for AD lesional, ex-lesional and non-lesional skin and for healthy skin. Only proteins measurable in both matrices are tested; column `direction` says "enriched in dISF", "enriched in serum" or "not different". Proteins measurable in dISF only are in `detection_by_matrix.csv`. Full lists per skin type are also in `Executive_summary_tables.xlsx` |
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

**Under review:** the pooled (all visits) relapse comparison also includes the lesional V1 samples.
Treat pooled relapse results with extra caution.

### `trajectories/` – disease course and relapse (step 11, aim 3)
Does the remaining lesional signal fade after clearing? Does it rise before a relapse? Results
are in `trajectory_results.csv` (number of significant proteins per model: `trajectory_summary.csv`).
`plots/` holds one figure for each top protein (the strongest results of steps 04, 05 and 11),
showing each patient over time (red = lesion site, blue = non-lesional, yellow = serum, dashed line
= relapse). If a severity file (SCORAD/EASI; see `data/README.md`) is present, the proteins that
follow the severity are listed too.

### `leip_reference/` – the population reference (step 08)
For the proteins where skin fluid and blood are correlated (step 06): the normal range in the LEIP
population, where the in-study serum samples (AD and healthy; MicroAD, RELAD and RELAD2) fall in
it, and whether the protein depends on age, sex, BMI, CRP, lipids … in healthy people.
`leip_reference.xlsx` holds the summary, the clinical associations and the check of LEIP against
the in-study healthy controls; the position of each sample is in `samples_vs_leip.csv`.

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

### `export/` – the data as CSV (step 19)
| File | Content |
|---|---|
| `ISF_NPX_wide.csv`, `Serum_NPX_wide.csv` | The Olink result: one row per sample, one column per protein (NPX as delivered) |
| `*_PCNormalizedNPX_wide.csv` | The same with the values used in the analysis |
| `samples.csv`, `proteins.csv` | Sample information and protein information (UniProt, LOD, detection) |
| `NPX_long.csv.gz` | Everything in one long table, including LOD and QC flags |
| `RELAD2/` | The RELAD2 samples only |

Tip: if Excel shows everything in one column, set `sep: ";"` under `export:` in `config.yml`
and rerun step 19 (`start_at <- 19; source("run_all.R")`).

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
| `AD_vs_HC` (serum) | AD patients vs healthy controls of the same studies (which studies: see the model name below) |
| `AD_vs_Biobank`, `HC_vs_Biobank` | vs LEIP biobank serum (the latter checks for biobank handling effects) |
| `active_vs_cleared` (serum) | Serum at visits with an active lesion vs after clearing |

**Model names (column `model`).** dISF: `states_all_visits` = all visits together; `baseline_V1` =
V1 only; `V1`…`V6` = one visit; `time_ex_lesional`, `time_non_lesional` = change per week; models
with `relapse` in the name = relapse analyses. Serum: the contrast `AD_vs_HC` comes from three
models: `AD_vs_HC_in_study` (all in-study controls), `MicroAD_AD_vs_HC` (MicroAD only, AD at V1;
used when dISF is compared with serum) and `RELAD_AD_vs_HC` (RELAD/RELAD2). Relapse: `MicroAD_relapse`,
`RELAD_relapse` (RELAD and RELAD2 together), `RELAD_only_relapse`, `RELAD2_only_relapse` and
`RELAD_relapse_unflagged` (without conflicting labels). All serum models are explained in
[`docs/STEPS.md`](docs/STEPS.md) (step 05).

---

## 6. Please keep in mind

- **Relapse results are exploratory.** There are only 4 relapsing patients in MicroAD.
- **Few samples per visit** (6–11 patients) make the FDR strict. "Not significant" does **not** mean "no difference".
- **NPX is relative:** compare a protein between groups, never protein A with protein B.
- **LEIP biobank serum** was collected and stored differently. Trust AD-vs-LEIP differences only if they also appear against the in-study healthy controls.
- **Metadata issues** (see `metadata/data_flags.csv`) are flagged, not fixed.
- **Under review:** the pooled relapse comparison in step 14 (`serum_vs_disf/`) also includes lesional V1 samples.
- **Under review:** the wording of key questions Q2–Q5 and the choice of one primary test per question (step 15).

---

## 7. What is in the repository

| Folder / file | Content |
|---|---|
| `oli_Jule.Rproj` | Open this in RStudio |
| `run_all.R` | Runs everything in order |
| `install_packages.R` | Installs the R packages (once) |
| `config.yml` | All settings: file paths, thresholds, focus proteins, mast cell markers, TNFRSF9 targets, CSV separator |
| `scripts/01_…` to `scripts/19_…` | One script per analysis step. The number is also the step in `start_at`. |
| `R/` | Shared functions used by the scripts |
| `tools/download_hpa.R` | Optional: downloads the public Human Protein Atlas table for step 17 |
| `data/` | Your input files (not uploaded to GitHub) – see `data/README.md` |
| `output/` | All results (not uploaded to GitHub) |
| `tests/` | Test with invented data |
| `docs/STEPS.md` | What each step does, in plain words |
| `docs/TECHNICAL.md` | Statistical methods and design decisions |

| Step | Script | What it does | Output (in `output/`) |
|---|---|---|---|
| 01 | `01_metadata.R` | Read and check the manifest | `metadata/` |
| 02 | `02_import_qc.R` | Read the Olink files, LOD, quality control | `qc/` |
| 03 | `03_explore.R` | PCA, sources of variation | `explore/` |
| 04 | `04_isf_models.R` | dISF comparisons | `models/ISF_*` |
| 05 | `05_serum_models.R` | Serum comparisons | `models/Serum_*`, `relad/` |
| 06 | `06_isf_vs_serum.R` | dISF–serum correlation | `isf_serum/` |
| 07 | `07_enrichment.R` | Pathways | `enrichment/` |
| 08 | `08_leip_reference.R` | LEIP population reference | `leip_reference/` |
| 09 | `09_isf_profile.R` | dISF proteome profile (aim 1) | `isf_profile/` |
| 10 | `10_matrix_comparison.R` | dISF vs serum proteome (aim 2) | `matrix_comparison/` |
| 11 | `11_trajectories.R` | Disease course (aim 3) | `trajectories/` |
| 12 | `12_focus_proteins.R` | Focus proteins (CD137 …) | `focus/` |
| 13 | `13_visit_course.R` | Per visit and time course | `visit_course/`, `models/ISF_by_visit_*` |
| 14 | `14_serum_vs_disf.R` | Serum vs dISF overlap per visit | `serum_vs_disf/` |
| 15 | `15_key_questions.R` | Answers to the key questions (mast cells, CD137, relapse, dISF vs serum) | `key_questions/` |
| 16 | `16_tnfrsf9_correlation.R` | Which proteins correlate with TNFRSF9 (CD137) in dISF (IL33, IL4 …) | `tnfrsf9_correlation/` |
| 17 | `17_disf_serum_signatures.R` | Do serum and dISF carry the same or different signatures? (MicroAD only) | `signatures/` |
| 18 | `18_summary_report.R` | Executive summary | `Executive_summary.pdf`, `Executive_summary_tables.xlsx` |
| 19 | `19_export_data.R` | CSV export of the data | `export/` |

---

## 8. Something went wrong?

| Message | Fix |
|---|---|
| `cannot open file 'R/utils.R'` / "Working directory must be the project folder" | Open `oli_Jule.Rproj`, or use `setwd()` to go to the repository folder |
| `[MISSING] manifest` or `[MISSING] NPX parquet files` | Check the file names and places in `data/README.md` (the `.parquet` files go into `data/npx/`) |
| `[--] Olink fixed LOD file` at the start, or `WARNING: no Olink fixed LOD file` in step 02 | The Olink fixed LOD file was not found, so the LOD comes from the negative controls (less precise). Put it in `data/` as `Explore_HT_Fixed_LOD.csv`, or keep Olink's own name (e.g. `Explore HT_Fixed LOD.csv`): it must contain "Fixed LOD" or "Fixed_LOD", end in `.csv`, and be the only such file in the folder. |
| A step fails | Read the last lines of the message, fix the problem, and restart from that step with `start_at <- N; source("run_all.R")` |
| The run skips steps you wanted, or writes to `output_sim/` | `start_at` or the test settings are still set from earlier: restart R (RStudio: Session → Restart R) and run again |
| Package error | Run `source("install_packages.R")` again. It only installs missing packages; to update an old one, use `install.packages("name")` or `BiocManager::install("name")`. |

When asking for help, send the error message and the names of files, columns or samples,
**not the data itself**.
