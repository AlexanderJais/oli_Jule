# LEIP cohort – galanin

An analysis of the **LEIP biobank sera** (population controls) that were measured on the Olink
Explore HT panel together with the O-MicroAD samples. It asks three questions about **galanin**, and which proteins go with each clinical parameter.

> **Start here:** after a run, open **`output/LEIP_galanin_report.pdf`**. Its first pages give the
> answers; the other pages show the evidence.

The dISF/serum analysis of the O-MicroAD study (skin fluid and serum of AD patients) is not part of
this branch; it is kept on the branch `claude/affectionate-davinci-8bnkoy`. Both use the same Olink
data files.

---

## 1. The questions

| | Question | How it is answered | Folder |
|---|---|---|---|
| **1** | **Does Olink confirm our galanin ELISA?** | Galanin was measured twice in the same sera: ELISA (pg/mL) and Olink (assay GAL). Do both rank the persons the same way (also within each sex and each plate)? Also: is GAL above its LOD; does the ELISA follow Olink GAL more than any other protein; do both correlate with the same proteins; do the values differ between plates; and, as a **benchmark**, how well do the other lab assays (ApoA-I, ApoB, Lp(a), A-FABP …) agree with Olink in the same sera – which also shows whether the clinical file and the Olink samples are correctly matched, overall and **per person** (a single swapped sample barely changes the overall correlations); and what the ELISA follows instead (proteins, gene sets, clinical parameters) | `output/1_elisa_validation/` |
| **2** | **Olink galanin on its own (the ELISA is ignored)** | **2a** Olink GAL vs every clinical parameter – unadjusted, adjusted for age, sex and plate, within each sex, leaving out one person at a time; which associations are robust; technical factors (plate, flagged samples). **2b** Olink GAL vs every measurable Olink protein; pre-specified: proteins stored and released together with galanin (chromogranins, secretogranins, NPY …); pathways (GSEA); the main axes of the serum proteome (principal components); do GAL's strongest partners form groups (heatmap)? | `output/2_olink_galanin/` |
| **3** | **Galanin and HDL – could galanin bind to HDL?** | Galanin vs HDL-C, ApoA-I and the other lipids – for all persons, **adjusted for sex** (women have higher HDL) and within each sex; is it specific to HDL (vs LDL, triglycerides, ApoB)?; galanin vs the HDL proteins measured by Olink (APOA1, APOA2, APOM, LCAT, PON1/3 …); does GAL behave like an HDL-associated protein?; **do the ELISA and Olink disagree more when HDL is high** (as expected if one assay does not see HDL-bound galanin, or HDL interferes with one assay)? Added after the first results (exploratory): does the **agreement weaken as HDL rises** (ELISA × HDL interaction), and is the disagreement larger in either direction? | `output/3_galanin_hdl/` |
| **4** | **Which Olink proteins go with each clinical parameter?** | Every measurable protein vs every analysed clinical parameter, unadjusted and adjusted for age, sex and plate; for each parameter the list of proteins at FDR < 0.05; how much of it is sex and age (many parameters differ between women and men or change with age and then share those proteins); which associations hold after adjusting. Also a data check: known associations (leptin–BMI, cystatin C–eGFR …) must show up | `output/4_clinical_screen/` |

Correlations can show whether the data **fit** the idea that galanin binds HDL; they cannot prove
it. The report lists experiments that could (galanin in lipoprotein fractions, ApoA-I pull-down,
ELISA vs Olink in HDL-depleted serum).

---

## 2. Data

Put the files into the `data/` folder as described in [`data/README.md`](data/README.md):
the Olink serum NPX file (it contains the LEIP samples), the Olink fixed LOD file, the LEIP clinical
file and, optionally, the manifest. **The column `Olink_SampleID` of the clinical file defines which
samples are LEIP samples.** Data files are never uploaded to GitHub.

Column names of the clinical file are used in plain ASCII: `Ins0_µU_ml` becomes `Ins0_uU_ml`,
`lipämisch` becomes `lipamisch`. `sex_MF` becomes `sex_male` (1 = man, 0 = woman) and the ELISA
column becomes `galanin_elisa`.

---

## 3. How to run

1. Open **`oli_Jule.Rproj`** in RStudio (sets the working folder).
2. The first time only: `source("install_packages.R")`.
3. Run **`source("run_all.R")`**. It takes a few minutes. The first lines show whether all input files were found (`[ok]` / `[MISSING]`).
4. To restart from a later step, e.g. after changing a setting for question 3: `start_at <- 4; source("run_all.R")`.

To check that everything works, with invented data: `source("tests/test_leip.R")` (2–4 minutes;
must end with "All LEIP galanin checks passed"). The invented data contain known effects (e.g.
galanin follows HDL within each sex) that the analysis must find.

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_data.R` | Olink values of the LEIP samples (LOD, QC), clinical data, sample checks |
| 2 | `scripts/02_elisa_validation.R` | Question 1 |
| 3 | `scripts/03_olink_galanin.R` | Question 2 (Olink galanin only) |
| 4 | `scripts/04_galanin_hdl.R` | Question 3 |
| 5 | `scripts/05_clinical_screen.R` | Question 4 (proteins vs all clinical parameters) and sanity checks |
| 6 | `scripts/06_report.R` | PDF report and `answers.csv` |

---

## 4. What is where – the `output/` folder

| File / folder | Content |
|---|---|
| `LEIP_galanin_report.pdf` | Answers, figures, key tables, methods and caveats |
| `answers.csv` | One line per answer: the verdict and the numbers behind it |
| `0_data/` | `samples.csv` (Olink QC + clinical data per sample), `sample_checks.csv` (failed QC, no clinical data, plate/ID mismatches …), `detection.csv` (per protein: share above LOD in LEIP), `parameters.csv` (which clinical variables are analysed and why the others are not) |
| `1_elisa_validation/` | `GAL_Olink_vs_ELISA.png`, `agreement_z_scores.png` (where do the methods disagree?), `lab_vs_Olink_benchmark.png`, `ELISA_vs_all_proteins.png`, `protein_profiles_GAL_vs_ELISA.png`; `agreement.csv`, `plate_effects.csv`, `sample_identity.csv` (per person: does the Olink sample fit the person's lab values?), `ELISA_vs_clinical.csv` and `ELISA_gene_sets.csv` (what the ELISA follows instead); all tables in `elisa_validation.xlsx` |
| `2_olink_galanin/` | 2a: `olink_galanin_vs_clinical.png`, `olink_galanin_top_clinical.png`, `olink_galanin_vs_clinical.csv`, `technical_factors.csv`; 2b: `olink_galanin_vs_all_proteins.png`, `olink_galanin_top_proteins.png`, `top_proteins_heatmap.png`, `neuroendocrine_proteins.png`, `gene_sets.png`, `proteome_axes.png`, `olink_galanin_vs_proteins.csv`; all tables in `olink_galanin.xlsx` |
| `3_galanin_hdl/` | `galanin_vs_HDL_by_sex.png`, `galanin_vs_lipids.png`, `galanin_vs_HDL_proteins.png`, `HDL_vs_all_proteins.png`, `assay_discordance_vs_HDL.png`, `agreement_by_HDL.png` (Olink vs ELISA at lower and higher HDL); `discordance_vs_HDL.csv`, `agreement_by_HDL.csv`; tables in `galanin_hdl.xlsx` |
| `4_clinical_screen/` | `per_parameter.csv` (for each clinical parameter: every protein at FDR < 0.05, unadjusted and adjusted), `n_significant_per_parameter.png`, `heatmap_key_parameters.png`, `volcano_key_parameters.png`; `associations_all.csv.gz` (every protein × every parameter), `clinical_screen.xlsx` (sheet `significant`: every significant pair with all statistics) |

---

## 5. How to read the results

| Term | Meaning |
|---|---|
| **rho** | Spearman correlation, from −1 to +1. +0.5: the protein tends to be high where galanin is high. Ranks are used, so outliers and skewed values do no harm. |
| **p** | How surprising the correlation would be by chance, for **one** test. The galanin questions were chosen in advance, so their p-values are the main result. |
| **FDR** | p-value corrected for testing many proteins at once (Benjamini–Hochberg). For the proteome-wide lists: significant = FDR < 0.05. |
| **adjusted** | Partial Spearman correlation: the same correlation after removing the influence of age, sex and Olink plate (for the ELISA: age and sex). |
| **robust** (2a) | p < 0.05 unadjusted and adjusted in the same direction, the same direction in women and in men, and still p < 0.05 whichever person is left out. An exploratory filter, not a proof. |
| **principal component** | One of the main axes along which the serum proteome varies between the persons (e.g. plate, sex, age, inflammation). If GAL follows one, it shares that source of variation. |
| **women only / men only** | The correlation within one sex: a relation that holds in both sexes is not a sex effect. |
| **NPX** | Olink's relative value (log2). It cannot be compared with pg/mL; ELISA and Olink can only agree in **ranking** the persons. |
| **discordance** | z(ELISA) − z(Olink GAL): positive when the ELISA is relatively higher than Olink for that person. Its size (absolute value) is the disagreement in either direction. |
| **interaction** | Does the agreement of ELISA and Olink change with HDL? Negative: the two agree less where HDL is high. The slopes show the agreement at low (−1 SD) and high (+1 SD) HDL, on ranks (close to rho). |
| **exploratory** | Added after seeing the first results: a hint to be confirmed, not a test of the original question. |

With 34 persons only fairly strong correlations are detectable: about |rho| ≥ 0.34 for p < 0.05,
and about |rho| ≥ 0.6 to pass the FDR over ~3,000 proteins. With thousands of proteins, single
correlations of |rho| ≈ 0.5 also occur by chance. "Not significant" does not mean "no correlation".

---

## 6. Please keep in mind

- **Sex:** women have higher HDL and may have higher galanin. Always look at the sex-adjusted and within-sex results before concluding that galanin follows HDL.
- **Different forms of galanin:** the Olink assay targets the galanin precursor (UniProt P22466, which also contains GMAP); a galanin ELISA may detect the mature peptide or fragments. Poor agreement does not automatically mean that one method is wrong.
- **Biobank serum** was not collected for peptide measurements (no protease inhibitors); galanin is a small, easily degraded peptide.
- **Benchmark:** if other lab assays agree well with Olink but galanin does not, the problem is specific to galanin; if all agree poorly, the samples themselves are the issue.
- The sample checks (`output/0_data/sample_checks.csv`) are flagged, not corrected. Samples failing Olink QC are left out (`drop_failed_samples` in the config).

---

## 7. What is in the repository

| Folder / file | Content |
|---|---|
| `run_all.R` | Runs all steps in order |
| `config.yml` | All settings: file paths, galanin assay and ELISA column, covariates, HDL proteins and lipids, pre-specified neuroendocrine proteins, lab-vs-Olink pairs, sanity checks, which clinical variables are tested |
| `scripts/` | One script per step (the number is also the step in `start_at`) |
| `R/utils.R` | Settings, output files, logging |
| `R/olink.R` | Reading the Olink files and the LOD (the LOD routine of the O-MicroAD analysis, unchanged) |
| `R/leip.R` | Clinical file, choice of clinical variables, (partial) Spearman correlations, bootstrap, answers and figures |
| `R/report.R` | PDF pages |
| `data/` | Your input files (not uploaded to GitHub) – see `data/README.md` |
| `output/` | All results (not uploaded to GitHub) |
| `tests/` | `make_test_data.R` (invented data with known effects) and `test_leip.R` (runs everything and checks the effects are found) |

---

## 8. Something went wrong?

| Message | Fix |
|---|---|
| "Working directory must be the project folder" | Open `oli_Jule.Rproj`, or use `setwd()` to go to the repository folder |
| `[MISSING] Olink NPX parquet file(s)` / `LEIP clinical file` | Check the file names and places in `data/README.md` |
| "None of the LEIP samples of the clinical file is in the NPX files" | The `Olink_SampleID` values of the clinical file must match the SampleIDs of the Olink file |
| A step fails | Read the last lines of the message, fix the problem, and restart from that step with `start_at <- N; source("run_all.R")` |
| Package error | Run `source("install_packages.R")` again |

When asking for help, send the error message and the names of files, columns or samples,
**not the data itself**.
