# What each step does – in plain words

The pipeline has 20 steps (`scripts/01_…` to `scripts/19_…`, plus `02b`). `run_all.R` runs them in this
order, and each step uses what the earlier ones saved. For every step this page explains the
**question** it answers, **what it does**, **what comes out**, and **what to keep in mind**.

Folder names refer to the `output/` folder. Words such as NPX, LOD, FDR and logFC are explained in
the [README](../README.md#2-words-you-will-meet). The statistical details are in
[`TECHNICAL.md`](TECHNICAL.md).

| Part | Steps |
|---|---|
| A. Preparing the data | 01 metadata · 02 import and QC · 02b QC overview · 03 first look |
| B. Group comparisons, protein by protein | 04 dISF models · 05 serum models · 07 pathways |
| C. Skin fluid and blood together | 06 correlation · 08 LEIP reference · 10 aim 2 (dISF vs serum proteome) · 14 overlap per visit |
| D. The three study aims | 09 aim 1 (dISF profile) · 11 aim 3 (disease course) · 13 visit by visit · *aim 2 = step 10, see part C* |
| E. Specific proteins and questions | 12 focus proteins · 15 key questions · 16 TNFRSF9 correlations · 17 serum vs dISF signatures |
| F. Reporting | 18 executive summary · 19 data export |

The steps are grouped by topic, so below they are not in number order; `run_all.R` always runs
them 01 → 19.

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
- It **stops** if a manifest column is missing, if a sample ID or a plate position appears twice, if a sample sits in a different well than in the plate layout, if SampleType is not dISF or Serum, or if a patient has two dISF samples from the same site at one visit.
- It **flags but does not change** contradictions, such as a relapse label that doesn't fit TimeToRelapse, dates that differ within one visit or are out of order, or a low volume (dISF < 20 µL, serum < 40 µL).

**Output:** `sample_metadata.csv` (one row per sample), `data_flags.csv` (what looks inconsistent).

### Step 02 – Import and quality control (`02_import_qc.R` → `qc/`)
**Question:** Are the Olink measurements usable, and which proteins can we measure at all?

**What it does**
1. Reads the two Olink files, one for dISF (its file name says "ISF") and one for serum. It checks that every sample is in the manifest and in the right file.
2. Works out the **limit of detection (LOD)** for every value, using Olink's fixed LOD file (the same method as Olink's own software). Without that file it falls back to the negative controls, which is less precise.
3. **Sample QC:**
   - Samples that Olink marks as *FAIL* are removed.
   - *WARN* samples and statistical outliers (median or spread of the sample far from the others) are **flagged, not removed**.
4. Removes proteins that Olink excluded (no values at all).
5. **Detection filter:** a protein is analysed in a matrix only if it is above LOD in at least 50 % of the samples of **at least one group**. In dISF a group is, e.g., "AD lesional"; in serum, e.g., "AD". So a protein that is only measurable in lesional skin is kept.
6. Checks the reproducibility of the control samples (CV) and draws a PCA for each matrix, coloured by plate.

**Output:**
- `sample_qc.csv` and `assay_detection.csv` (which proteins are kept);
- `lod.csv` and `pca_by_plate.png`;
- the cleaned data for all later steps (`data/*.rds`).

**Keep in mind:** the analysis uses **PCNormalizedNPX**, because the plates were filled by matrix: plate 1 holds dISF, plate 2 both matrices and plates 3–4 serum. `pca_by_plate.png` shows one PCA per matrix, so within a panel the points should not separate by plate. If they do, check which groups sit on which plate.

### Step 02b – QC overview: how many proteins are measurable? (`02b_qc_overview.R` → `qc_overview/`)
**Question:** How many of the Olink proteins are above LOD in serum and in dISF, and which kinds of proteins are they?

**What it does**
- Calls each protein **above LOD** in a matrix if at least 50 % of that matrix's samples are above its LOD
  (`qc$min_detect_frac`). Serum = MicroAD samples, the same people as dISF (settings in `config.yml` → `qc_overview`).
- Draws **ring charts** (serum left, dISF right): share of proteins above / below LOD, total in the middle.
- Below them, the **distribution across 14 protein classes** from the Human Protein Atlas: enzymes, transcription
  factors, nuclear receptors, GPCRs, voltage-gated ion channels, transporters, drug related, disease related,
  cancer related, immune related, essential (DepMap), intracellular, membrane, extracellular / secreted.
  A protein can belong to several classes, so the bars do not add up to the total.

**Output:** `qc_overview.pdf` and `.png` (font Nimbus Sans, Helmholtz Munich violet / pink), `qc_overview.xlsx`
(numbers, one row per protein with its classes, counts under other above-LOD definitions, proteins not found in the atlas).

**Keep in mind:** the protein classes need the Human Protein Atlas table (`source("tools/download_hpa.R")` once);
without it only the ring charts are drawn. The counts can differ slightly from the "proteins kept" of step 02,
which uses a more generous rule (≥ 50 % in at least one group); the Excel file shows both. If `qc_overview: rule` is set to
`analysis_filter`, the step 02 rule is used, and then serum covers all cohorts, because that filter is computed over all of them.

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

Repeated samples of the same person are taken into account (patient as a "random effect"). The
plate is included in every model except `relapse_delta_xL_minus_NL` (step 04): there the two
samples compared come from the same visit and normally the same plate, so their difference
already removes the plate effect.

### Step 04 – dISF comparisons (`04_isf_models.R` → `models/ISF_*`)
**Question:** Which proteins differ between skin states in the skin fluid?

| Model | Compares |
|---|---|
| `states_all_visits` | All visits together: lesional vs non-lesional, ex-lesional vs non-lesional, lesional vs ex-lesional, lesional and non-lesional AD skin vs healthy skin, AD lesional vs CPUO lesional, CPUO lesional vs non-lesional |
| `baseline_V1` | V1 only: lesional vs non-lesional, lesional and non-lesional AD skin vs healthy skin, CPUO lesional vs non-lesional |
| `time_ex_lesional`, `time_non_lesional` | Change per week in cleared and in non-lesional AD skin |
| `relapse_ex_lesional` | Cleared skin of patients who relapse later vs those who don't (exploratory) |
| `relapse_delta_xL_minus_NL` | The same, on the difference cleared minus non-lesional skin at the same visit (removes plate and day-to-day effects) |

**Output:**
- `ISF_results.csv`: every protein × every comparison of the first five models;
- `ISF_relapse_delta_results.csv`: the `relapse_delta_xL_minus_NL` model;
- `ISF_summary.csv`: the number of significant proteins;
- volcano plots in `models/volcano/`.

### Step 05 – Serum comparisons (`05_serum_models.R` → `models/Serum_*`)
**Question:** Which proteins differ in blood?

| Model | Compares |
|---|---|
| `AD_vs_HC_in_study` | AD vs the healthy controls of the studies (one sample per person) |
| `AD_vs_Biobank`, `HC_vs_Biobank` | vs LEIP biobank serum. The second shows what differs just because biobank serum was handled differently. |
| `MicroAD_active_vs_cleared` | The same patients' serum when their lesion is active vs cleared |
| `MicroAD_AD_vs_HC` | MicroAD only: AD (V1) vs healthy volunteers. **This is the serum reference whenever dISF is compared with serum** (step 10, the dISF–serum correlation in step 12, Q5 of step 15). The other dISF–serum analyses (steps 06, 14, 16, 17) also use MicroAD serum only. Steps 12 and 15 also show serum results that include RELAD/RELAD2 next to it (in step 15: Q1 and Q4). |
| `MicroAD_relapse` | MicroAD relapsers vs non-relapsers, at visits when the lesion is cleared (exploratory) |
| `RELAD_relapse`, `RELAD_relapse_unflagged` | RELAD + RELAD2 relapsers vs non-relapsers; the second leaves out samples with conflicting relapse labels |
| `RELAD_only_relapse`, `RELAD2_only_relapse`, `RELAD_AD_vs_HC` | Each RELAD cohort on its own; AD vs healthy within RELAD/RELAD2 |

**Output:**
- `Serum_results.csv` and `Serum_summary.csv`;
- `Serum_AD_vs_controls_agreement.csv`: proteins that differ in AD against **both** the in-study and the biobank controls. These are the most trustworthy.
- `relad/RELAD_RELAD2_serum_results.xlsx`: **all RELAD and RELAD2 results in one workbook**: every model, group means, % above LOD, sample labels and NPX values.

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

**What it does:** pairs every MicroAD dISF sample with the serum of the same patient and visit, then computes two correlations per protein:
- **within patients** (`r_within`): do dISF and serum rise and fall together over the visits?
- **between patients** (`r_between`): do patients with high dISF levels also have high serum levels?

Both are done separately for the lesion site (AD lesion site plus CPUO lesional skin) and the
non-lesional site (AD and CPUO non-lesional skin plus the healthy volunteers' skin). Only AD
patients have several visits, so only they contribute to `r_within`. `r_between` also includes
CPUO patients and healthy volunteers, so it can partly reflect group differences.

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
1. **Detection** (MicroAD samples only): a protein counts as measurable in a matrix if it is above LOD in ≥ 50 % of all MicroAD samples of that matrix (not per group as in step 02). Each protein is then classed as measurable in dISF only, serum only, both, or neither.
2. **Relative enrichment:** only proteins measurable in both dISF and serum are tested; proteins measurable in dISF only appear as "dISF only" in `detection_by_matrix.csv`. For each matched dISF–serum pair it computes dISF minus serum and centres that on the typical protein, because dISF is more dilute overall. Proteins far above the typical protein are relatively enriched in skin fluid, i.e. candidates for **local production in the skin**. This is done separately for AD lesional, ex-lesional and non-lesional skin and for healthy skin. "Enriched" requires FDR < 0.05 and at least a 2-fold difference.
3. **Disease signals:** it compares the effects of step 04 (skin) with step 05 (blood, MicroAD serum only): same direction, skin only, or blood only?
4. Does the serum level follow the lesional-minus-non-lesional skin difference from visit to visit?

**Output:** `matrix_comparison.xlsx`, `detection_by_matrix.csv`, `relative_enrichment.csv`, `disease_signal_concordance.*`.

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
**Question (aim 1):** Which proteins are measurable in dISF, and how does that depend on the skin state?

**What it does**
- For each protein: how often it is above LOD in dISF, overall and per skin state. It is then classed as *robust* (≥ 90 %), *detected*, *sporadic* or *not detected*.
- **Lesion-restricted proteins:** detectable in lesional AD skin but not in non-lesional or healthy skin.
- Which pathways the detectable dISF proteome covers, compared with the whole panel.
- Two heatmaps (z-score per protein; numbers in `isf_profile.xlsx`):
  - the 50 proteins that vary most between dISF samples, selected without group information. Proteins and samples are clustered, with colour bars for skin state, patient and visit. It shows which proteins move together and whether samples group by skin state or by patient.
  - the 25 strongest up and 25 strongest down proteins of lesional vs non-lesional skin (step 04), in lesional, ex-lesional, non-lesional and healthy skin. It shows the disease signal and whether ex-lesional skin still carries it.

**Keep in mind:** NPX values of different proteins are not comparable, so this step does not rank proteins by "amount".

### Step 11 – Disease course and relapse (`11_trajectories.R` → `trajectories/`)
**Question (aim 3):** How does the skin fluid change after the lesion clears, and before a relapse?

**What it does**
- **A.** Cleared skin still differs from non-lesional skin (a "molecular scar"). Does that difference fade with the weeks since clearing?
- **B.** In relapsers, does it rise in the weeks before the relapse? (exploratory)
- **C.** The same question for serum.
- **D.** If a severity file (SCORAD, EASI …) is present: which proteins follow the severity?
- A figure per top protein showing every patient over time (red = lesion site, blue = non-lesional, gold = serum, dashed line = relapse).

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

### Step 16 – Which proteins correlate with TNFRSF9 in dISF? (`16_tnfrsf9_correlation.R` → `tnfrsf9_correlation/`)
**Question:** Which measured proteins correlate with TNFRSF9 (4-1BB / CD137) in dISF? In particular, do IL33 (IL-33), IL4 (IL-4), CSF2, IL6, IL18, CXCL8, IL1RL1, KIT, KITLG, TPSAB1, TPSB2 and FCER1A? The list is in `config.yml`.

**What it does**
- **Detectability (reported first):** % of dISF samples above LOD per site and group. Proteins below LOD in more than half of the samples are flagged before any correlation. Values below LOD are used as measured, not replaced, and every targeted correlation is repeated on the samples where both proteins are above LOD.
- **Targeted correlations:** Spearman ρ, p and n for TNFRSF9 vs each target, separately for the lesion site, for non-lesional / healthy skin, and for all dISF together (AD and healthy; CPUO samples are not used in this step). Each is reported:
  - (a) pooled;
  - (b) within each group: AD, healthy, relapse, non-relapse, lesional, ex-lesional;
  - (c) per visit.

  A correlation is flagged "pooled only" when the pooled ρ has p < 0.05 but no within-group ρ (per-visit results not counted) reaches p < 0.05 in the same direction. It then probably comes from group differences rather than from the two proteins being linked; small groups can also miss a real link.
- **Proteome-wide:** TNFRSF9 against every dISF protein, adjusted for skin state / group, visit and plate.
  - The partial Spearman correlation gives the ranking: the correlation that is left after removing differences due to skin state, visit and plate.
  - A mixed model with the patient as random effect is a second method.
  - The output is a ranked list with FDR, showing where IL33, IL4 and the other targets fall.
- **Within patients** (AD patients): do changes of TNFRSF9 from visit to visit go with changes of the targets in the same patient? Two methods: repeated-measures correlation (each patient's values around their own average) and Δ–Δ correlation (the change from one visit to the next).
- **dISF vs serum:** dISF TNFRSF9 vs serum TNFRSF9 at the same visit, and whether the TNFRSF9–target relationship exists in serum. MicroAD only.

**Output:**
- `TNFRSF9_correlation.xlsx`;
- scatter plots per target (coloured by group, per site, open symbols = below LOD);
- a heatmap and a proteome-wide plot.

**Keep in mind:**
- Pooled correlations use several samples of the same patient, so trust the within-group and within-patient results more.
- The p-values of the targeted correlations are not corrected for testing several targets.
- TPSB2 is not on the Explore HT panel.

### Step 17 – Serum vs dISF signatures (`17_disf_serum_signatures.R` → `signatures/`)
**Question:** Do serum and dISF reflect the same biology at different sensitivity, or different processes? Does that change between groups (AD vs healthy; relapse vs non-relapse) and across V1–V6?

**Only MicroAD samples are used.** RELAD, RELAD2 and LEIP serum are excluded. "Measured in both" = passes the detection filter of step 02 in dISF and in serum, **and** is above LOD in at least 50 % of the samples of at least one MicroAD serum group (AD or healthy). dISF proteins not detected in serum are reported separately, so that "dISF only" can be told apart: is there no effect in serum, or is the protein simply not measurable there?

1. **Effect-size concordance:** for every protein measured in both, the serum effect is plotted against the dISF effect, per comparison, per visit and pooled. Reported:
   - Spearman ρ;
   - the slope: below 1 = the same signal, weaker in serum;
   - % same direction.

   This is done for all proteins and again for those significant in either compartment.
2. **Time-resolved models:** group × visit in each compartment.
   - Healthy controls have one visit, so the test asks whether the AD–healthy difference changes over V1–V6. Relapse × visit is tested within AD.
   - dISF-significant proteins are grouped by their time profile: high at V1 and resolving, persistent, late-rising, reversing or fluctuating.
   - For each profile, the step checks whether serum shows a weaker version of it.
3. **Compartment × group:** which proteins change with disease (AD vs healthy) differently in dISF than in serum. Tested on the dISF − serum difference of each paired sample (same person, same visit), for all visits together and for V1 alone.
4. **Paired correlation within patients** (from step 06):
   - proteins that follow serum are candidates for systemic spill-over;
   - proteins that don't follow serum and are relatively enriched in dISF are candidates for local production.
5. **Signature sets** (from the all-visits results): dISF-only, serum-only, shared-concordant and shared-discordant, at FDR < 0.05 (p < 0.05 as sensitivity analysis). Each set gets:
   - Reactome / GO enrichment, with the measured panel as background;
   - tissue origin from the Human Protein Atlas (skin / keratinocyte, immune, liver);
   - the full protein lists.
6. **Relapse, predictive:** dISF lesion site, dISF non-lesional skin, their difference (ex-lesional minus non-lesional) and serum, at the **visit before the relapse** vs the cleared visits of non-relapsers, one value per patient.

**Output:**
- `signatures.xlsx` (sheet `README` explains each sheet);
- concordance and time-profile figures;
- a relapse volcano plot.

**Keep in mind:**
- Relapse analyses rest on 4 vs 6 patients.
- The tissue origin needs the Human Protein Atlas table: download it once with `source("tools/download_hpa.R")`, then rerun from step 02b with `start_at <- "02b"; source("run_all.R")` (this also adds the protein classes of step 02b). Without the table the origin columns stay empty.

---

## F. Reporting

### Step 18 – Executive summary (`18_summary_report.R` → `Executive_summary.pdf`, `Executive_summary_tables.xlsx`)
Collects everything in one PDF:
1. overview: samples, data and QC, then the proteins above / below LOD by protein class (step 02b);
2. answers to the key questions with evidence, then the TNFRSF9 correlations (step 16) and the serum vs dISF signatures (step 17);
3. key findings per aim;
4. all comparisons in one table;
5. visit course, volcano plots, pathways, dISF vs serum, serum, serum vs dISF per visit (step 14) and focus proteins;
6. methods and caveats.

Where the PDF shows a shortened list, the complete list is in `Executive_summary_tables.xlsx`. The first sheet, `index`, explains each sheet. If an earlier step did not run, its page says so instead of stopping the report.

### Step 19 – Data export (`19_export_data.R` → `export/`)
Writes the Olink data as CSV files for use in Excel, Prism, SPSS …:
- all samples × all proteins, as delivered (NPX) and as analysed (PCNormalizedNPX);
- sample and protein information;
- a long table with LOD and QC flags;
- the RELAD2 samples separately.

Values below LOD are exported as measured, which is Olink's recommendation; the long table marks them.
