# What each step does – in plain words

The pipeline has 17 steps (`scripts/01_…` to `scripts/17_…`). `run_all.R` runs them in this
order, and each step uses what the earlier ones saved. For every step this page explains the
**question** it answers, **what it does**, **what comes out**, and **what to keep in mind**.

Folder names refer to the `output/` folder. Words such as NPX, LOD, FDR and logFC are explained in
the [README](../README.md#2-words-you-will-meet). The statistical details are in
[`TECHNICAL.md`](TECHNICAL.md).

| Part | Steps |
|---|---|
| A. Preparing the data | 01 metadata · 02 import and QC · 03 first look |
| B. Group comparisons, protein by protein | 04 dISF models · 05 serum models · 07 pathways |
| C. Skin fluid and blood together | 06 correlation · 08 LEIP reference · 10 dISF vs serum proteome · 14 overlap per visit |
| D. The three study aims | 09 aim 1 (dISF profile) · 10 aim 2 · 11 aim 3 (disease course) · 13 visit by visit |
| E. Specific proteins and questions | 12 focus proteins · 15 key questions |
| F. Reporting | 16 executive summary · 17 data export |

---

## A. Preparing the data

### Step 01 – Sample information (`01_metadata.R` → `metadata/`)
**Question:** Which sample is what?

**What it does**
- Reads the sample manifest.
- Gives each sample its patient, visit, matrix (dISF or serum), cohort and group: AD, healthy, CPUO or LEIP biobank.
- For dISF, records the skin site (tracked lesion site **L** or non-lesional site **NL**) and the skin state (lesional, ex-lesional, non-lesional).
- Unifies spelling variants, e.g. "remisison" → remission, "GS"/"H" → healthy.
- Works out the course of each AD patient:
  - days since V1;
  - at which visit the lesion relapsed (the first time it is lesional again after it had cleared);
  - which visits were before the relapse.
- Adds sex and the sample volume, plus the LEIP clinical data (age, BMI, CRP …).

**Checks**
- It stops if a sample ID appears twice, or if a sample sits in a different well than in the plate layout.
- It **flags but does not change** contradictions, such as a relapse label that doesn't fit, a missing date or a low volume.

**Output:** `sample_metadata.csv` (one row per sample), `data_flags.csv` (what looks inconsistent).

### Step 02 – Import and quality control (`02_import_qc.R` → `qc/`)
**Question:** Are the Olink measurements usable, and which proteins can we measure at all?

**What it does**
1. Reads the two Olink files (ISF and serum). It checks that every sample is in the manifest and in the right file.
2. Works out the **limit of detection (LOD)** for every value, using Olink's fixed LOD file (the same method as Olink's own software).
3. **Sample QC:**
   - Samples that Olink marks as *FAIL* are removed.
   - *WARN* samples and statistical outliers (median or spread of the sample far from the others) are **flagged, not removed**.
4. Removes proteins that Olink excluded (no values at all).
5. **Detection filter:** a protein is analysed in a matrix only if it is above LOD in at least 50 % of the samples of **at least one group**. In dISF a group is, e.g., "AD lesional"; in serum, e.g., "AD". So a protein that is only measurable in lesional skin is kept.
6. Checks the reproducibility of the control samples (CV) and draws a PCA by plate.

**Output:**
- `sample_qc.csv` and `assay_detection.csv` (which proteins are kept);
- `lod.csv` and `pca_by_plate.png`;
- the cleaned data for all later steps (`data/*.rds`).

**Keep in mind:** the analysis uses **PCNormalizedNPX**. Plates were filled by matrix (plate 1 dISF, plates 3–4 serum), so the plates *will* differ in the PCA; that is expected.

### Step 03 – First look at the data (`03_explore.R` → `explore/`)
**Question:** What drives the differences between samples: the person, the skin state, the plate or the cohort?

**What it does**
- **PCA** for each matrix, coloured by skin state, group, visit, plate, cohort or relapse. Samples that look alike lie close together.
- **Variance partitioning:** for every protein, the share of its variation explained by:
  - in dISF: the person, the skin state and the plate;
  - in serum: the group, the cohort and the plate.

**How to read it:** if "person" explains most of the variation, every comparison must account for repeated samples of the same patient, and later steps do.

---

## B. Group comparisons, protein by protein

All comparisons follow the same scheme. For **each protein separately**, a statistical model
compares two groups:
- *logFC* is the difference;
- *P.Value* is the p-value for that one protein;
- *adj.P.Val* is the FDR, corrected for testing thousands of proteins.

Repeated samples of the same person are taken into account (patient as a "random effect"), and the plate is always included in the model.

### Step 04 – dISF comparisons (`04_isf_models.R` → `models/ISF_*`)
**Question:** Which proteins differ between skin states in the skin fluid?

| Model | Compares |
|---|---|
| `states_all_visits` | All visits together: lesional vs non-lesional, ex-lesional vs non-lesional, lesional vs ex-lesional, AD skin vs healthy skin, AD vs CPUO, CPUO lesional vs non-lesional |
| `baseline_V1` | The same at V1 only |
| `time_ex_lesional`, `time_non_lesional` | Change per week in cleared and in non-lesional AD skin |
| `relapse_ex_lesional` | Cleared skin of patients who relapse later vs those who don't (exploratory) |
| `relapse_delta_xL_minus_NL` | The same, on the difference cleared minus non-lesional skin at the same visit (removes plate and day-to-day effects) |

**Output:**
- `ISF_results.csv`: every protein × every comparison;
- `ISF_summary.csv`: the number of significant proteins;
- volcano plots in `models/volcano/`.

### Step 05 – Serum comparisons (`05_serum_models.R` → `models/Serum_*`)
**Question:** Which proteins differ in blood?

| Model | Compares |
|---|---|
| `AD_vs_HC_in_study` | AD vs the healthy controls of the studies (one sample per person) |
| `AD_vs_Biobank`, `HC_vs_Biobank` | vs LEIP biobank serum. The second shows what differs just because biobank serum was handled differently. |
| `MicroAD_active_vs_cleared` | The same patients' serum when their lesion is active vs cleared |
| `MicroAD_relapse`, `RELAD_relapse` | Relapsers vs non-relapsers (MicroAD; RELAD + RELAD2), plus a version without samples with contradictory labels |

**Output:**
- `Serum_results.csv` and `Serum_summary.csv`;
- `Serum_AD_vs_controls_agreement.csv`: proteins that differ in AD against **both** the in-study and the biobank controls. These are the most trustworthy.

### Step 07 – Pathways (`07_enrichment.R` → `enrichment/`)
**Question:** Which biological processes stand behind the protein changes?

**What it does:** for every comparison of steps 04 and 05, it sorts all proteins from "most up" to "most down". It then checks whether the proteins of a known pathway sit at the top or the bottom of that list (GSEA). The pathways come from:
- Hallmark, Reactome and GO;
- a custom set of AD/Th2 proteins (CCL17, CCL22, IL13, POSTN …).

**Output:** `gsea_results.csv`. `NES > 0` means the pathway is higher in the first group; it counts if `padj < 0.05`.

---

## C. Skin fluid and blood together

### Step 06 – Do skin fluid and blood go together? (`06_isf_vs_serum.R` → `isf_serum/`)
**Question:** If a protein is high in a patient's skin fluid, is it also high in their blood?

**What it does:** pairs every dISF sample with the serum of the same patient and visit, then computes two correlations per protein:
- **within patients** (`r_within`): do dISF and serum rise and fall together over the visits?
- **between patients** (`r_between`): do patients with high dISF levels also have high serum levels?

Both are done separately for the lesion site and the non-lesional site.

**Output:** `isf_serum_correlation.csv`; `significant_proteins.csv` (used by step 08).

### Step 08 – The population reference LEIP (`08_leip_reference.R` → `leip_reference/`)
**Question:** For proteins linked between skin and blood, what is "normal", and what else influences them?

**What it does:** for each protein from step 06, in the 35 LEIP biobank samples:
- the normal range (5th–95th percentile) and where AD patients fall in it;
- whether the protein depends on age, sex, BMI, CRP, lipids … in healthy people. Such factors could confound the AD results.
- whether LEIP differs from the in-study healthy controls (biobank handling).

**Output:** `leip_reference.xlsx`.

### Step 10 – The dISF proteome vs the serum proteome (`10_matrix_comparison.R` → `matrix_comparison/`)
**Question (aim 2):** What does skin fluid contain that blood doesn't, and do skin and blood show the same disease signals?

**What it does**
1. **Detection:** which proteins are measurable in dISF only, serum only, both, or neither.
2. **Relative enrichment:** for each matched dISF–serum pair it computes dISF minus serum and centres that on the typical protein, because dISF is more dilute overall. Proteins far above the typical protein are relatively enriched in skin fluid, i.e. candidates for **local production in the skin**. This is done separately for AD lesional, ex-lesional and non-lesional skin and for healthy skin. "Enriched" requires FDR < 0.05 and at least a 2-fold difference.
3. **Disease signals:** it compares the effects of step 04 (skin) with step 05 (blood): same direction, skin only, or blood only?
4. Does the serum level follow the lesional-minus-non-lesional skin difference from visit to visit?

**Output:** `matrix_comparison.xlsx`, `relative_enrichment.csv`, `disease_signal_concordance.*`.

### Step 14 – Serum vs dISF per visit (`14_serum_vs_disf.R` → `serum_vs_disf/`)
**Question:** If we ask the same question in skin fluid and in blood, how much do the answers overlap, and what does only dISF show?

**What it does:** asks two questions in both matrices, at each visit and pooled over all visits:
- **AD vs healthy:** dISF lesion site or non-lesional skin vs healthy skin; serum AD vs healthy.
- **relapse vs non-relapse:** using only visits before the relapse.

Each protein then falls into one category:
- significant in both (same or opposite direction);
- dISF only;
- dISF only because the protein isn't measurable in serum;
- serum only.

This is done at FDR < 0.05 and at p < 0.05 (exploratory).

**Output:**
- Venn diagrams and stacked bars ("red = information only dISF provides");
- coloured volcano plots;
- `overlap_summary.csv`.

**Keep in mind:** at p < 0.05 about 5 % of all proteins are significant by chance alone. The pooled relapse comparison is currently under review, because it also includes lesional V1 samples.

---

## D. The three study aims

### Step 09 – What is in skin fluid? (`09_isf_profile.R` → `isf_profile/`)
**Question (aim 1):** Which proteins are measurable in dermal ISF, and how does that depend on the skin state?

**What it does**
- For each protein: how often it is above LOD in dISF, overall and per skin state. It is then classed as *robust* (≥ 90 %), *detected*, *sporadic* or *not detected*.
- **Lesion-restricted proteins:** detectable in lesional AD skin but not in non-lesional or healthy skin.
- Which pathways the detectable dISF proteome covers, compared with the whole panel.
- A heatmap of the 50 proteins that vary most between dISF samples.

**Keep in mind:** NPX values of different proteins are not comparable, so this step does not rank proteins by "amount".

### Step 11 – Disease course and relapse (`11_trajectories.R` → `trajectories/`)
**Question (aim 3):** How does the skin fluid change after the lesion clears, and before a relapse?

**What it does**
- **A.** Cleared skin still differs from non-lesional skin (a "molecular scar"). Does that difference fade with the weeks since clearing?
- **B.** In relapsers, does it rise in the weeks before the relapse? (exploratory)
- **C.** The same question for serum.
- **D.** If a severity file (SCORAD, EASI …) is present: which proteins follow the severity?
- A figure per top protein showing every patient over time (red = lesion site, blue = non-lesional, yellow = serum, dashed line = relapse).

**Output:** `trajectory_results.csv`, `plots/`.

### Step 13 – Visit by visit (`13_visit_course.R` → `visit_course/`)
**Question:** What does the lesion site look like at each visit, and which proteins stay changed throughout?

**What it does:** at each visit V1–V6 (with at least 5 patients) it compares:
- lesion site vs healthy skin (`Lsite_vs_HC`);
- lesion site vs non-lesional skin of the same patients (`Lsite_vs_NL`);
- non-lesional vs healthy skin (`NL_vs_HC`).

Healthy skin (one visit) is the reference at every visit. **"Regulated at all visits"** means significant at every analysed visit, in the same direction.

**Output:**
- a volcano plot per visit (`models/volcano/ISF_by_visit_V*.png`);
- time-course plots and heatmaps;
- `visit_course.xlsx`.

**Keep in mind:**
- At V1 the lesion site is lesional; after clearing it is ex-lesional (see `lesion_site_state_per_visit.csv`).
- With 6–11 patients per visit, "not significant" doesn't mean "no difference".
- The fallback list "p < 0.05 at every visit" is exploratory.

---

## E. Specific proteins and questions

### Step 12 – Focus proteins (`12_focus_proteins.R` → `focus/`)
**Question:** What can we say about the proteins we chose in advance? These are CD137 (TNFRSF9), CD137L (TNFSF9), KITLG, CPA4, FCER1A, TPSAB1, MS4A2, TPSD1, PNOC and POSTN; the list is in `config.yml`.

**What it does:** for each protein, even if it failed the detection filter:
- how often it is measurable, per matrix and group;
- the same comparisons as steps 04 and 05, for this protein alone;
- correlation with serum, severity (if available) and LEIP clinical data;
- every result the other steps produced for it;
- figures: skin states, V1 lesional vs non-lesional per patient, serum groups, the course per patient, and dISF vs serum.

**How to read it:** because these proteins were chosen beforehand, their own **p-value** is the main result. The proteome-wide FDR is shown alongside.

### Step 15 – Key questions (`15_key_questions.R` → `key_questions/`)
**Question:**
- Q1 Mast cell markers in AD and/or relapse?
- Q2 CD137 / CD137L as a mast cell marker?
- Q3 CD137 / CD137L and relapse?
- Q4 Marker or predictor of relapse?
- Q5 Is dISF superior to serum?

**What it does:** builds a **mast cell score**: each marker is standardised, and the score is a sample's average over all measured markers. It then answers each question with pre-specified tests and writes a one-line verdict plus the numbers behind it. An evidence figure per protein shows the data behind each statement.

**Output:** `answers.csv`, `key_questions.xlsx`, `evidence/*.png`.

**Keep in mind:** the verdicts are generated automatically from p < 0.05 and must be read with the evidence. Relapse results rest on 4 vs 6 patients in MicroAD. The wording of Q2–Q5 and the choice of one primary test per question are currently under review.

---

## F. Reporting

### Step 16 – Executive summary (`16_summary_report.R` → `Executive_summary.pdf`, `Executive_summary_tables.xlsx`)
Collects everything in one PDF:
1. answers to the key questions with evidence;
2. data and QC;
3. key findings per aim;
4. all comparisons in one table;
5. visit course, volcano plots, pathways, dISF vs serum, serum and focus proteins;
6. methods and caveats.

Where the PDF shows a shortened list, the complete list is in `Executive_summary_tables.xlsx`. The first sheet, `index`, explains each sheet. If an earlier step did not run, its page says so instead of stopping the report.

### Step 17 – Data export (`17_export_data.R` → `export/`)
Writes the Olink data as CSV files for use in Excel, Prism, SPSS …:
- all samples × all proteins, as delivered (NPX) and as analysed (PCNormalizedNPX);
- sample and protein information;
- a long table with LOD and QC flags;
- the RELAD2 samples separately.

Values below LOD are exported as measured, which is Olink's recommendation; the long table marks them.
