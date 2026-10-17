# O-MicroAD – Olink proteomics of skin fluid and blood

This repository contains the R scripts that analyse the **Olink Explore HT** data of the O-MicroAD
study: about 5,400 proteins measured in **dermal interstitial fluid (dISF)** and in **serum** of
patients with atopic dermatitis (AD) and control persons.

This page tells you how to run the analysis and **where to find which result**. Two more pages:
- [`docs/STEPS.md`](docs/STEPS.md) – what each step does (01–19, plus 02b, and the dedicated analyses), in plain words;
- [`docs/TECHNICAL.md`](docs/TECHNICAL.md) – technical details (models, design decisions);
- [`docs/CHANGES.md`](docs/CHANGES.md) – what changed in each version and which steps to rerun.

> **Which script do I run? Only `source("run_all.R")`.** It runs everything: all steps, the QC overview, visit 3,
> the V1 relapse screen, the relapse follow-up, the per-visit analysis V1–V4 and the threshold check (these right after the
> import, so their folders appear in the first minutes), and at the end the executive summary.
>
> **Sample manifest:** `data/2026-10-15_manifest.xlsx` (set in `config.yml`). Every analysis takes the samples from it
> (step 01); if the file is missing, the run stops – an older manifest is never used instead.
>
> **Lost in the output folder?** Open **`output/00_FOLDER_GUIDE.pdf`**: one page, one line per folder.
>
> **Start here:** after a run, open **`output/Executive_summary.pdf`**. It summarises all
> results (the answers to the key questions come right after the overview page) and tells you
> which folder holds the details.

**Dedicated analyses** – all included in `run_all.R`. Each also has its own script, **only** as a shortcut to rerun that one
analysis quickly (import / QC plus that analysis, a few minutes) – you never need these for a full run:

| Question | Shortcut (optional) | Open first |
|---|---|---|
| How many proteins are measurable in serum and dISF, which only in dISF, and in which pathways? | `source("run_qc_overview.R")` | `output/qc_overview/qc_overview.pdf`, `dISF_only_proteins.xlsx`, `dISF_only_pathways.pdf`, `dISF_only_KEGG.pdf` |
| Visit 3: lesional vs healthy / non-lesional skin, AD vs healthy serum – and what does dISF show that serum does not? | `source("run_visit3.R")` | `output/visit3/V3_dISF_vs_serum.pdf`, `V3_DEGs_dISF.xlsx`, `V3_DEGs_serum.xlsx` |
| Which markers at V1 predict a later relapse – in dISF, in serum, and in both? | `source("run_relapse_v1.R")` | **`output/relapse_V1/V1_relapse_executive_summary.pdf`** |
| Relapse vs non-relapse, follow-up: is it calculated right, CD45 / T cells, IL4 · CCL27 · CCL5 · CCL8 over the visits, volcano plots with every p < 0.05 labelled, druggable targets? | `source("run_relapse_followup.R")` | **`output/relapse_followup/Relapse_followup_executive_summary.pdf`** |
| Relapse vs non-relapse and lesion site vs non-lesional at each visit V1–V4 (DEG + GSEA, volcano plots for PowerPoint), AKR1B10 and the mTORC1 gene set over the visits? | `source("run_visits_v1_v4.R")` | **`output/visits_V1_V4/Visits_V1-V4_summary.pdf`**, `1_relapse_volcano_V1.png` |
| Are the thresholds too strict? What do the LOD filter, FDR, minimum group sizes and sample QC remove – measured on the data at other values | `source("run_thresholds.R")` | **`output/thresholds/Thresholds_summary.pdf`**, `thresholds.xlsx` |

Details: [section 4](#qc_overview--how-many-proteins-are-measurable-step-02b) (QC overview), [section 8](#8-separate-analysis-visit-3-run_visit3r) (visit 3), [section 9](#9-separate-analysis-v1-relapse-screen-run_relapse_v1r) (relapse), [section 10](#10-separate-analysis-relapse-follow-up-run_relapse_followupr) (relapse follow-up), [section 11](#11-separate-analysis-relapse-and-lesional-signatures-per-visit-v1v4-run_visits_v1_v4r) (per visit V1–V4), [section 12](#12-separate-analysis-thresholds-run_thresholdsr) (thresholds).

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
| **Relapser** | An AD patient whose lesion site, after clearing, becomes lesional again (MicroAD, manifest of 2026-10-15: 5 patients, relapse at V3, V4, V5, V6 and V6; they leave the study at relapse; 5 non-relapsers). |
| **Leave one out** | Each patient is predicted by a rule (marker, cut-off) set on all other patients only – an honest test of prediction in small groups. |
| **Exact relabelling / chance check** | The analysis is repeated for every possible way to choose the relapse group (252 for 5 of 10 patients). It shows how many "hits" random labels give, and an exact p-value (smallest possible 1/252). |
| **GSEA / over-representation** | Gene set analyses: GSEA asks whether a gene set (pathway) sits at the top or bottom of all proteins ranked by their change; over-representation asks whether a list (e.g. the dISF-only proteins) contains more members of a set than expected. Gene sets: MSigDB Hallmark, Reactome, GO biological processes, and KEGG. |

---

## 3. How to run it

1. Open **`oli_Jule.Rproj`** in RStudio. This sets the working folder correctly.
2. Put the data files into the `data/` folder, as described in [`data/README.md`](data/README.md).
   **Required:** the manifest, and the two Olink NPX files (`.parquet`) in the subfolder `data/npx/`.
   **Recommended:** the Olink fixed LOD file. **Optional:** the LEIP clinical data and the severity scores.
   **Data files are never uploaded to GitHub.**
3. The first time only, run `source("install_packages.R")` (needs R 4.3 or newer).
4. Optional, once (needs internet): `source("tools/download_hpa.R")` downloads the public Human Protein
   Atlas table. Step 02b uses it for the protein classes, step 17 for the tissue origin of proteins, the relapse
   screen for drug targets. Without it, these steps still run but leave these parts out.
   `source("tools/download_kegg.R")` downloads the current KEGG human pathways (free for academic use) for the
   KEGG analyses; without it the older MSigDB KEGG collection is used. Both files stay in `data/reference/`.
5. Run `source("run_all.R")`. This is the only script you need: it runs every step and every dedicated analysis
   (QC overview, visit 3, V1 relapse screen, relapse follow-up, per visit V1–V4, thresholds). The full run can take an hour or more,
   and the results appear in `output/`. Optional shortcut to rerun just one dedicated analysis: its own script (`run_qc_overview.R`,
   `run_visit3.R`, `run_relapse_v1.R`, `run_relapse_followup.R`, `run_visits_v1_v4.R`, `run_thresholds.R`; see the table at the top) –
   a few minutes each.
6. **To resume** from a later step, e.g. after changing a setting for step 12: `start_at <- 12; source("run_all.R")`.
   Step numbers are the script numbers; for step 02b use `start_at <- "02b"`.
   Earlier results are reused, and all later steps run again (including the summary, step 18). Start
   from the **earliest** step you changed (step numbers are in [section 7](#7-what-is-in-the-repository)).
   After a successful run `start_at` is removed automatically, so the next `source("run_all.R")` is a full run.
   (If a run stops with an error, `start_at` stays set: fix the problem and run again.)
   Resuming only works with results made from the **same manifest**: every script run is logged with the manifest it
   used (`output/metadata/steps_run.csv`), and a resume stops if earlier results come from another manifest. After a
   new manifest, run everything once with `source("run_all.R")`.

The first lines of the run list the input files: `[ok]` = found, `[MISSING]` = required file
missing (the run stops), `[--]` = optional or recommended file not found (the run continues).

To check that R and the packages work, run `source("tests/test_pipeline.R")`. It uses invented data
only, takes about 15 minutes, and must end with "All pipeline checks passed". The dedicated analyses have their own
quick tests: `tests/test_visit3.R`, `tests/test_relapse_v1.R`, `tests/test_relapse_followup.R`, `tests/test_visits_v1_v4.R` and `tests/test_thresholds.R` (a few minutes each). Its results go to
`output_sim/` and are **invented**: never use them as study results. The test switches its settings
off again when it finishes; the first lines of every run say which settings file and output folder are used.

### Updating to a new version
Pull the new version (RStudio: **Git** tab → **Pull**), then check [`docs/CHANGES.md`](docs/CHANGES.md):
it says whether new packages are needed (`source("install_packages.R")`) and from which step to rerun.

---

## 4. What is where – the `output/` folder

Each script writes into its own folder. Most folders contain an **`.xlsx` file** that collects
all their tables. Open that one first.

### `Executive_summary.pdf` – the overview (step 18)
First an overview of the samples and the quality control, the proteins above / below LOD by protein class
(step 02b), then the answers to the key questions
with the data behind them, the TNFRSF9 correlations (step 16) and the serum vs dISF signatures
(step 17), key numbers and findings for each aim, all comparisons in one table, the most important
figures (including the relapse predictors over the visits: `run_relapse_v1.R` adds them and rebuilds this summary), and the methods and caveats. **Read this first.**

### `Executive_summary_tables.xlsx` – the full lists behind the summary (step 18)
Wherever the PDF shows a shortened list ("up: TNC, LAIR2 …"), the complete list is in this
workbook. The first sheet, `index`, says what each sheet contains:
- proteins above / below LOD per protein class, serum and dISF (step 02b, sheet `detection_by_class`);
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
(public data, saved in `data/reference/`), then rerun from step 02b (`start_at <- "02b"; source("run_all.R")`),
which also adds the protein classes of step 02b and updates the summary. Without the file, the tissue origin is left out.

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
| `manifest_used.txt` | Which manifest the results come from (file name, checksum, date of step 01). Later steps stop if the manifest has changed since. |
| `steps_run.csv` | Every script run, with the manifest it used and when. A resume (`start_at`) stops if earlier results come from another manifest. |

### `qc_overview/` – how many proteins are measurable? (step 02b)
Only this figure: `source("run_qc_overview.R")` (runs steps 01, 02 and 02b, nothing else).
One page (`qc_overview.pdf` / `.png`): ring charts with the share of Olink proteins above and below LOD
(total in the middle; serum left, dISF right), and below them the distribution across 14 protein classes
(enzymes, transcription factors … extracellular / secreted). A protein counts as above LOD if it is above
LOD in at least 50 % of the samples of that matrix; serum = MicroAD samples (settings: `config.yml` →
`qc_overview`). The classes come from the Human Protein Atlas; a protein can belong to several classes.
`qc_overview.xlsx` has the numbers, one row per protein with its classes, and the counts under other
above-LOD definitions (e.g. the analysis filter of step 02).
**Which proteins are measurable in dISF only?** Based on matched samples: every dISF sample with a serum sample of
the same person and visit, and those serum samples (dISF against serum only, skin types pooled).
`dISF_only_proteins.xlsx` lists every protein above LOD in ≥ 50 % of the matched dISF but < 50 % of the
matched serum samples (sheet `dISF_only`; `by_class`: all of them grouped by protein class; one sheet per class; `serum_only`);
`dISF_only_proteins.pdf` shows them per class (% of serum vs % of dISF samples above LOD, the clearest labelled).
**Which pathways?** Sheet `pathways` and `dISF_only_pathways.pdf`: gene sets (MSigDB Hallmark, Reactome, GO:BP)
over-represented among the dISF-only proteins, tested against all proteins measured in both matrices (not the genome).
The same within transcription factors, enzymes and intracellular proteins (background = the measured proteins of that
class): sheets `pathways_transcription_factors`, `pathways_enzymes`, `pathways_intracellular` and
`dISF_only_pathways_<class>.pdf` (classes: `config.yml` → `qc_overview$pathway_classes`).
**KEGG:** the same analyses with the current KEGG pathways (372 human pathways; download once with
`source("tools/download_kegg.R")`, which stores the KEGG release date; without it the older MSigDB KEGG collection
with 186 pathways is used) in their own files:
`dISF_only_KEGG.xlsx` (sheets `all_classes`, `transcription_factors`, `enzymes`, `intracellular`) and
`dISF_only_KEGG.pdf`, `dISF_only_KEGG_<class>.pdf`.

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
| `ISF_results.xlsx`, `Serum_results.xlsx`, `ISF_by_visit_results.xlsx` | **The same results in Excel: one sheet per volcano panel** (e.g. `V1 Lsite_vs_HC`), all proteins sorted by p-value; the `index` sheet lists all sheets with the number of significant proteins |
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
| `top_variable_heatmap.png` | The 50 dISF proteins that vary most between samples, chosen **without** using any group information. Rows and samples are clustered (similar ones side by side); colour bars show skin state, patient and visit. Use it to see which proteins move together and whether samples group by skin state or by patient. |
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
| `fonts/` | Nimbus Sans font files (free, URW) used for the QC overview figure |
| `tools/download_hpa.R` | Optional: downloads the public Human Protein Atlas table for steps 02b and 17 |
| `run_visit3.R`, `visit3/` | Separate visit 3 analysis (see section 8) |
| `run_relapse_v1.R`, `relapse_v1/` | Separate V1 relapse screen (see section 9) |
| `run_relapse_followup.R`, `relapse_followup/` | Separate relapse follow-up: eight questions on relapse vs non-relapse (see section 10) |
| `run_visits_v1_v4.R`, `visits_v1_v4/` | Separate per-visit analysis V1–V4: relapse vs non-relapse, lesion site vs non-lesional, AKR1B10, mTORC1 (see section 11) |
| `run_thresholds.R`, `thresholds/` | Separate threshold check: how strict each threshold is and what it removes (see section 12) |
| `tools/download_kegg.R` | Optional: downloads the current KEGG human pathways (free for academic use) for the KEGG analysis of step 02b |
| `data/` | Your input files (not uploaded to GitHub) – see `data/README.md` |
| `output/` | All results (not uploaded to GitHub) |
| `tests/` | Tests with invented data: `test_pipeline.R` (all steps), `test_visit3.R`, `test_relapse_v1.R`, `test_relapse_followup.R`, `test_visits_v1_v4.R`, `test_thresholds.R` |
| `docs/STEPS.md` | What each step does, in plain words |
| `docs/CHANGES.md` | What changed, and what to rerun after an update |
| `docs/TECHNICAL.md` | Statistical methods and design decisions |

| Step | Script | What it does | Output (in `output/`) |
|---|---|---|---|
| 01 | `01_metadata.R` | Read and check the manifest | `metadata/` |
| 02 | `02_import_qc.R` | Read the Olink files, LOD, quality control | `qc/` |
| 02b | `02b_qc_overview.R` | Proteins above / below LOD in serum and dISF, by protein class | `qc_overview/` |
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

## 8. Separate analysis: visit 3 (`run_visit3.R`)
Part of `run_all.R`; the shortcut `source("run_visit3.R")` runs only the import / QC and this analysis. MicroAD, visit 3
(`config.yml` → `visit3`). Healthy volunteers are sampled once, so their samples are the healthy reference.

| Matrix | Comparisons (log2FC) |
|---|---|
| dISF | lesion site vs healthy skin (lesional − healthy); lesion site vs non-lesional skin of the same patients (lesional − non-lesional, paired) |
| serum | AD vs healthy (serum has no skin site, so no lesional vs non-lesional comparison) |

The lesion site is taken as sampled at visit 3, whatever its clinical state then (the Excel README gives how many
are still lesional and how many cleared). Results in `output/visit3/`:
- `V3_DEGs_dISF.xlsx`, `V3_DEGs_serum.xlsx` – differential proteins (limma), one sheet per comparison, all tested
  proteins sorted by p with log2FC, FDR, means and % above LOD per group.
- `volcano/` – one volcano plot per comparison, the 25 most regulated proteins labelled.
- `V3_GSEA.xlsx`, `gsea/` – GSEA per comparison (MSigDB Hallmark, Reactome, GO:BP, and KEGG).
- `V3_dISF_vs_serum.xlsx` / `.pdf` – AD vs healthy in both matrices: which proteins change only in dISF (and whether
  they are measurable in serum at all), the fold changes side by side, and which gene sets are significant only in dISF.

## 9. Separate analysis: V1 relapse screen (`run_relapse_v1.R`)
Which markers at V1 (baseline, before the topical steroids at V2) predict a later relapse – in dISF, in serum, and in
both? The same screen runs on the dISF lesion site and on the serum of the same patients; markers found in both
matrices make the strongest case. `source("run_relapse_v1.R")` runs only the import / QC and this screen → `output/relapse_V1/`, and then rebuilds the main
`Executive_summary.pdf` with the predictor plots (if `run_all.R` has been run once before).
MicroAD has 5 relapsers (relapse at V3, V4, V5, V6, V6) and 5 non-relapsers (manifest of 2026-10-15), so this is a screen to pick candidates:
- `V1_relapse_targets.xlsx`: README – the answer in one sentence; `can_V1_predict` – can V1 dISF predict relapse
  at all: patients predicted correctly when the marker is picked and its cut-off set without the predicted patient
  (leave one out), against all 252 random labellings, and how many proteins separate perfectly by chance;
  `targets_up` – the ranked candidate list (per marker: patients predicted correctly leave-one-out, sensitivity,
  specificity, the V1 cut-off in NPX, correlation with how soon the relapse comes) (V1 lesion site p < 0.05 and exact
  relabelling p < 0.05 over all 252 ways to choose 5 of the 10 patients; ranked by AUC, still higher at V2 after
  treatment, stable when any one patient is left out), with flags: also higher in V1 serum, approved / potential drug
  target, secreted / membrane (Human Protein Atlas); `targets_down`; `chance_check` – how many such proteins random
  labels give; full results for V1 lesion site, V1 non-lesional skin, V1 lesion minus non-lesional, V2 cleared
  lesion site and V1 serum. Serum: `serum_targets_up` / `serum_targets_down`, `serum_can_V1_predict`;
  **`dISF_and_serum`** – markers in both matrices (and markers in one supported in the other), with `overlap_chance`
  (how many proteins random labels make "higher in both"). `V1_relapse_GSEA.xlsx` also lists the gene sets
  changed in the same direction in dISF and serum (`gene_sets_dISF_and_serum`).
- **`V1_relapse_executive_summary.pdf`** – start here: the answer, the prediction test, the ranked predictors
  (cut-off, sensitivity / specificity, V2, serum, drug target), their course over the visits, the gene sets higher
  in relapsers at V1 and V2, and the methods.
- `V1_relapse_dISF_vs_serum.pdf` (fold changes in both matrices, the shared markers labelled),
  `V1_relapse_markers_dISF_and_serum_course.pdf`, `V1_relapse_volcano_serum.pdf`,
  `V1_relapse_volcano.pdf`, `V1_relapse_targets_course.pdf` (each target per patient over V1–V6, relapse marked),
  `V1_relapse_predictors_higher_over_visits.pdf` / `V1_relapse_predictors_lower_over_visits.pdf` (the predictors higher / lower in
  relapsers over V1–V6: mean ± SE at the lesion site and in non-lesional skin, relapsers and non-relapsers, grey band = healthy skin;
  data in `V1_relapse_predictors_over_visits.xlsx`),
  `V1_relapse_chance_check.pdf`, `V1_relapse_prediction_test.pdf`, `V1_relapse_GSEA.xlsx` / `.pdf`.

## 10. Separate analysis: relapse follow-up (`run_relapse_followup.R`)
Eight follow-up questions on relapse vs non-relapse (the same patients as section 9: MicroAD, 5 relapsers, 5 non-relapsers).
`source("run_relapse_followup.R")` runs only the import / QC and this analysis → `output/relapse_followup/`. Start with
**`Relapse_followup_executive_summary.pdf`**: page 1 answers all eight questions in short, then one section per question.
Settings: `config.yml` → `relapse_followup` (the selected markers, the gene sets checked, the p-value for the volcano labels).

| # | Question | Files |
|---|---|---|
| 1 | Was relapse vs non-relapse calculated correctly – and is the lower immune signalling in relapsers (TCR, death receptor, NF-κB / kinase gene sets) real? Labels vs skin course, direction of log2FC / NES, the proteins behind these gene sets (each tested on its own), overall sample level, AD lesional vs healthy skin as positive control | `1_check_relapse_vs_nonrelapse.xlsx`, `1_immune_leading_edge_proteins.pdf`, `1_sample_level_check.pdf` |
| 2 | How the dISF lesion site comparison at V1 (5 later relapse vs 5 non-relapse) was calculated, in plain words, with a worked example | executive summary, section 2 |
| 3 | Leukocyte content: CD45 (PTPRC) and other leukocyte proteins, each on its own | `3_4_leukocytes_T_cells.xlsx`, `3_leukocyte_markers.pdf`, `3_CD45_all_samples.pdf` |
| 4 | T cells: LCK, CD3 and the other selected T-cell and activation markers, each on its own; each T-cell marker relative to CD45 (T-cell share of leukocytes) | `3_4_leukocytes_T_cells.xlsx`, `4_T_cell_markers.pdf`, `4_T_cell_markers_per_CD45.pdf`, `4_T_cell_activation_markers.pdf` |
| 5 | IL4, CCL27, CCL5, CCL8 at each visit: heatmap of the group means (y = protein, x = visit, relapse and non-relapse side by side) and per patient (clustered) | `5_cytokines_over_visits.xlsx`, `5_cytokines_heatmap_group_means.pdf`, `5_cytokines_heatmap_per_patient.pdf` |
| 6 | Volcano dISF lesion site V1, every protein with p < 0.05 labelled, with the plotted data – and whether the GSEA uses it (yes: all tested proteins ranked by the `t` of the same fit) | `6_volcano_dISF_lesion_V1.pdf`, `6_volcano_dISF_lesion_V1_data.xlsx` |
| 7 | Volcano serum V1, likewise | `7_volcano_serum_V1.pdf`, `7_volcano_serum_V1_data.xlsx` |
| 8 | Druggable targets that could lengthen the relapse-free time: higher in relapsers at V1 and higher in lesional **and** non-lesional AD skin than in healthy skin; drug-target flags (Human Protein Atlas) and known drugs (`relapse_followup/ad_drugs.csv`) | `8_druggable_targets.xlsx`, `8_druggable_targets.pdf` |

Olink measures proteins in the fluid, not cells: CD45, LCK and CD3 are proxies of cell content, not cell counts.

## 11. Separate analysis: relapse and lesional signatures per visit, V1–V4 (`run_visits_v1_v4.R`)
dISF of the MicroAD AD patients, at each visit V1–V4 (`config.yml` → `visits_v1_v4`). Part of `run_all.R`; the shortcut
`source("run_visits_v1_v4.R")` runs only the import / QC and this analysis → `output/visits_V1_V4/`. Start with
**`Visits_V1-V4_summary.pdf`**. The figures are **PNG files for PowerPoint** (16:9, 300 dpi, large fonts – insert them as
they are) and the same as PDF.

| # | Question | Files |
|---|---|---|
| 1 | Relapse vs non-relapse at the lesion site, at each visit (DEG, limma; exact relabelling check per visit): the relapse signature at V1, and which proteins differ at V2, V3, V4 – and at ≥ 2 visits | `1_relapse_vs_nonrelapse_DEG_V1-V4.xlsx` (one sheet per visit, `signature_over_visits`), `1_relapse_volcano_V1.png` … `V4`, `1b_relapse_signature_over_visits.png` |
| 2 | GSEA of the same comparisons, and the gene sets over the visits | `2_relapse_vs_nonrelapse_GSEA_V1-V4.xlsx`, `2_relapse_GSEA_V1.png` …, `2b_relapse_GSEA_over_visits.png` |
| 3 | Lesion site vs non-lesional skin of the same patients, at each visit (paired DEG) | `3_lesional_vs_nonlesional_DEG_V1-V4.xlsx`, `3_lesional_volcano_V1.png` … `V4` |
| 4 | AKR1B10 at the lesion site and in non-lesional skin over the visits (mean ± SE, paired test per visit, grey band = healthy skin) | `4_AKR1B10_lesional_vs_nonlesional_over_visits.png`; sheets `AKR1B10`, `AKR1B10_per_patient` of file 3 |
| 5 | The mTORC1 gene set (Hallmark mTORC1 signalling), lesion site vs non-lesional at each visit: GSEA (NES, p, FDR) and each of its proteins | `5_lesional_vs_nonlesional_GSEA_V1-V4.xlsx` (`signature`, `signature_proteins`), `5_mTORC1_lesional_vs_nonlesional_by_visit.png`, `5b_mTORC1_proteins_by_visit.png` |

Relapsers' samples at their relapse visit (lesional again) are left out of the relapse comparison (`include_relapse_visit: false`),
so from V2 on all compared lesion sites are in remission: with the manifest of 2026-10-15 the groups are 5 vs 5 at V1 and V2,
4 vs 5 at V3 and 3 vs 5 at V4 (a visit is tested only with ≥ 3 patients per group). The volcano plots label the 20 proteins
with the lowest p (`volcano_labels`); every protein is in the Excel files. The proteins and gene sets tracked in 4 and 5
are set in `track_proteins` and `signature_sets`.

## 12. Separate analysis: thresholds (`run_thresholds.R`)
Are the thresholds too strict – could one of them remove the signal? Each threshold that removes samples or proteins or decides
a result is measured on the data, at the current value and at other values (`config.yml` → `thresholds`, which only sets the values
compared; the analyses use `qc`, `stats` and `enrichment`). Part of `run_all.R`; the shortcut `source("run_thresholds.R")` runs only
the import / QC and this check → `output/thresholds/`. Start with **`Thresholds_summary.pdf`**: page 1 has one line per threshold –
its setting, its effect on these data and a verdict.

| What | How it is checked | Files |
|---|---|---|
| Sample QC | samples removed (Olink SampleQC FAIL – a FAIL in one block removes the whole sample: how much of them really failed), kept but flagged (WARN, outliers at 2–4 SD, low volume); samples in the NPX files but not in the manifest (left out of every analysis); proteins dropped from the main models for > 20 % missing values | sheets `samples`, `outliers_by_sd`, `removed_samples`, `sample_id_mismatches`, `missing_values` |
| Detection (LOD) filter | assays kept per matrix at 0–90 % above LOD (in at least one detection group, as step 02, or of all samples); the proteins of interest of `config.yml` near or below the threshold | `1_assays_kept_by_detection_threshold.png`, `4_proteins_of_interest_above_LOD.png`; sheets `detection`, `detection_groups`, `proteins_of_interest` |
| Effect on results | five comparisons (relapse V1, lesion site vs non-lesional V1 and V3, AD vs healthy in dISF and serum) re-run at each detection threshold: proteins tested, p < 0.05, FDR < 0.05 / 0.1 / 0.2, chance level for relapse. Key test: do the proteins the filter removes carry more p < 0.05 than the 5 % chance gives? | `2_comparisons_by_detection_threshold.png`, `3_removed_proteins_vs_chance.png`; sheets `comparisons`, `removed_by_filter`, `removed_p_below_0.05` |
| Gene sets | GSEA at three detection thresholds (gene sets with FDR < 0.05 / 0.1 / 0.25, the mTORC1 set); gene sets eligible at other minimum sizes | sheets `gene_sets`, `signature_gene_set`, `gene_set_sizes` |
| Group sizes | patients per visit and group against `stats$min_group_n` and `visit_min_subjects`; the smallest exact p the group sizes allow | `5_group_sizes_per_visit.png`; sheet `group_sizes` |
| Relapse criteria (V1) | proteins left after each criterion of the V1 relapse screen (p, exact p, LOD in the relapsers), and with looser criteria | `6_relapse_criteria_V1.png`; sheet `relapse_funnel` |
| Inventory | every threshold in the code: value, `config.yml` key or hard-coded, where, what it does | sheet `inventory` |
| Inconsistencies | the same kind of threshold defined differently in different steps (e.g. what "detected" means, assays without LOD, minimum group sizes, before-relapse samples) – listed, not changed | sheet `inconsistencies`, last page of the PDF |

## 13. Something went wrong?

| Message | Fix |
|---|---|
| `cannot open file 'R/utils.R'` / "Working directory must be the project folder" | Open `oli_Jule.Rproj`, or use `setwd()` to go to the repository folder |
| `[MISSING] manifest` or `[MISSING] NPX parquet files` | Check the file names and places in `data/README.md` (the `.parquet` files go into `data/npx/`) |
| `[--] Olink fixed LOD file` at the start, or `WARNING: no Olink fixed LOD file` in step 02 | The Olink fixed LOD file was not found, so the LOD comes from the negative controls (less precise). Put it in `data/` as `Explore_HT_Fixed_LOD.csv`, or keep Olink's own name (e.g. `Explore HT_Fixed LOD.csv`): it must contain "Fixed LOD" or "Fixed_LOD", end in `.csv`, and be the only such file in the folder. |
| A step fails | Read the last lines of the message, fix the problem, and restart from that step with `start_at <- N; source("run_all.R")` |
| `[MISSING] manifest` although a manifest is in `data/` | The manifest must have exactly the name in `config.yml` (`data/2026-10-15_manifest.xlsx`); another file (e.g. an older `manifest.xlsx`) is never used instead. Rename the file or set `paths$manifest`. |
| `The results in output/ were made from the manifest …` or `Cannot resume: these results were made from another sample manifest` | The manifest changed since the earlier results were made. Run everything again: `rm(start_at); source("run_all.R")`. |
| The run skips steps you wanted, or writes to `output_sim/` | `start_at` or the test settings are still set from an interrupted run: restart R (RStudio: Session → Restart R) and run again. The first lines of the run show the settings file and output folder. |
| `there is no package called '…'` | A new version needs a new package: run `source("install_packages.R")`, restart R, run again |
| Package error | Run `source("install_packages.R")` again. It only installs missing packages; to update an old one, use `install.packages("name")` or `BiocManager::install("name")`. |

When asking for help, send the error message and the names of files, columns or samples,
**not the data itself**.
