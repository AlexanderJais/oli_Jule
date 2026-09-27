# LEIP cohort – galanin

A dedicated analysis of the **LEIP biobank sera** (population controls) that were measured on the
Olink Explore HT panel together with the O-MicroAD samples. It is independent of the dISF/serum
analysis in the rest of this repository: it only uses the same Olink data files.

> **Start here:** after a run, open **`leip_galanin/output/LEIP_galanin_report.pdf`**. Its first
> pages give the answers; the other pages show the evidence.

---

## 1. The questions

| | Question | How it is answered | Folder |
|---|---|---|---|
| **1** | **Does Olink confirm our galanin ELISA?** | Galanin was measured twice in the same sera: ELISA (pg/mL) and Olink (assay GAL). Do both rank the persons the same way? Also: is GAL above its LOD; does the ELISA follow Olink GAL more than any other protein; do both correlate with the same proteins; and, as a **benchmark**, how well do the other lab assays (ApoA-I, ApoB, Lp(a), A-FABP …) agree with Olink in the same sera | `output/1_elisa_validation/` |
| **2** | **What does galanin correlate with?** | Olink GAL (and the ELISA) vs every measurable Olink protein; pre-specified: proteins stored and released together with galanin (chromogranins, secretogranins, NPY …); pathways (GSEA); clinical parameters | `output/2_galanin_correlates/` |
| **3** | **Galanin and HDL – could galanin bind to HDL?** | Galanin vs HDL-C, ApoA-I and the other lipids – for all persons, **adjusted for sex** (women have higher HDL) and within each sex; is it specific to HDL (vs LDL, triglycerides, ApoB)?; galanin vs the HDL proteins measured by Olink (APOA1, APOA2, APOM, LCAT, PON1/3 …); does GAL behave like an HDL-associated protein?; **do the ELISA and Olink disagree more when HDL is high** (as expected if one assay does not see HDL-bound galanin)? | `output/3_galanin_hdl/` |
| S | Supplementary: every protein vs every clinical parameter | Context (what else relates to HDL, sex, BMI …) and a data check: known associations (leptin–BMI, cystatin C–eGFR …) must show up | `output/4_clinical_screen/` |

Correlations can show whether the data **fit** the idea that galanin binds HDL; they cannot prove
it. The report lists experiments that could (galanin in lipoprotein fractions, ApoA-I pull-down,
ELISA vs Olink in HDL-depleted serum).

---

## 2. Data

The same files as for the main analysis, in the `data/` folder of the project
(paths in `leip_galanin/config.yml`):

| File | Content |
|---|---|
| `data/npx/*.parquet` | Olink NPX files. The LEIP samples are in the serum file; files without LEIP samples are skipped. |
| `data/Explore_HT_Fixed_LOD.csv` | Olink fixed LOD file (recommended). |
| `data/LEIP_clinical_parameters_n35.xlsx` | Sheet `Key_parameters` (incl. `Galanin [pg/mL]`) and sheet `All_SORB_parameters` (all other SORB variables). **Its column `Olink_SampleID` defines which samples are LEIP samples.** |
| `data/manifest.xlsx` | Optional: cross-check of the SubjectIDs. |

Column names of the clinical file are used in plain ASCII: `Ins0_µU_ml` becomes `Ins0_uU_ml`,
`lipämisch` becomes `lipamisch`. `sex_MF` becomes `sex_male` (1 = man, 0 = woman) and the ELISA
column becomes `galanin_elisa`. Data files are never uploaded to GitHub.

---

## 3. How to run

1. Open **`oli_Jule.Rproj`** in RStudio (sets the working folder).
2. The first time only: `source("install_packages.R")`.
3. Run **`source("leip_galanin/run.R")`**. It takes a few minutes. The first lines show whether all input files were found (`[ok]` / `[MISSING]`).
4. To restart from a later step, e.g. after changing a setting for question 3: `leip_start_at <- 4; source("leip_galanin/run.R")`.

To check that everything works, with invented data: `source("leip_galanin/tests/test_leip.R")`
(2–4 minutes; must end with "All LEIP galanin checks passed"). The invented data contain known
effects (e.g. galanin follows HDL within each sex) that the analysis must find.

| Step | Script | What it does |
|---|---|---|
| 1 | `01_data.R` | Olink values of the LEIP samples (LOD, QC), clinical data, sample checks |
| 2 | `02_elisa_validation.R` | Question 1 |
| 3 | `03_galanin_correlates.R` | Question 2 |
| 4 | `04_galanin_hdl.R` | Question 3 |
| 5 | `05_clinical_screen.R` | Supplementary screen and sanity checks |
| 6 | `06_report.R` | PDF report and `answers.csv` |

---

## 4. What is where – `leip_galanin/output/`

| File / folder | Content |
|---|---|
| `LEIP_galanin_report.pdf` | Answers, figures, key tables, methods and caveats |
| `answers.csv` | One line per answer: the verdict and the numbers behind it |
| `0_data/` | `samples.csv` (Olink QC + clinical data per sample), `sample_checks.csv` (failed QC, no clinical data, plate/ID mismatches …), `detection.csv` (per protein: share above LOD in LEIP), `parameters.csv` (which clinical variables are analysed and why the others are not) |
| `1_elisa_validation/` | `GAL_Olink_vs_ELISA.png`, `agreement_z_scores.png` (where do the methods disagree?), `lab_vs_Olink_benchmark.png`, `ELISA_vs_all_proteins.png`, `protein_profiles_GAL_vs_ELISA.png`; all tables in `elisa_validation.xlsx` |
| `2_galanin_correlates/` | `GAL_vs_all_proteins.png`, `GAL_top_proteins.png`, `neuroendocrine_proteins.png`, `gene_sets.png`, `galanin_vs_clinical.png` (+ the same for the ELISA); tables in `galanin_correlates.xlsx` |
| `3_galanin_hdl/` | `galanin_vs_HDL_by_sex.png`, `galanin_vs_lipids.png`, `galanin_vs_HDL_proteins.png`, `HDL_vs_all_proteins.png`, `assay_discordance_vs_HDL.png`; tables in `galanin_hdl.xlsx` |
| `4_clinical_screen/` | `n_significant_per_parameter.png`, `heatmap_key_parameters.png`, `volcano_key_parameters.png`; `associations_all.csv.gz` (every protein × every parameter), `clinical_screen.xlsx` |

---

## 5. How to read the results

| Term | Meaning |
|---|---|
| **rho** | Spearman correlation, from −1 to +1. +0.5: the protein tends to be high where galanin is high. Ranks are used, so outliers and skewed values do no harm. |
| **p** | How surprising the correlation would be by chance, for **one** test. The galanin questions were chosen in advance, so their p-values are the main result. |
| **FDR** | p-value corrected for testing many proteins at once (Benjamini–Hochberg). For the proteome-wide lists: significant = FDR < 0.05. |
| **adjusted** | Partial Spearman correlation: the same correlation after removing the influence of age, sex and Olink plate (for the ELISA: age and sex). |
| **women only / men only** | The correlation within one sex: a relation that holds in both sexes is not a sex effect. |
| **NPX** | Olink's relative value (log2). It cannot be compared with pg/mL; ELISA and Olink can only agree in **ranking** the persons. |
| **discordance** | z(ELISA) − z(Olink GAL): positive when the ELISA is relatively higher than Olink for that person. |

With 34 persons only fairly strong correlations are detectable: about |rho| ≥ 0.34 for p < 0.05,
and about |rho| ≥ 0.6 to pass the FDR over ~3,000 proteins. With thousands of proteins, single
correlations of |rho| ≈ 0.5 also occur by chance. "Not significant" does not mean "no correlation".

---

## 6. Please keep in mind

- **Sex:** women have higher HDL and may have higher galanin. Always look at the sex-adjusted and within-sex results before concluding that galanin follows HDL.
- **Different forms of galanin:** the Olink assay targets the galanin precursor (UniProt P22466, which also contains GMAP); a galanin ELISA may detect the mature peptide or fragments. Poor agreement does not automatically mean that one method is wrong.
- **Biobank serum** was not collected for peptide measurements (no protease inhibitors); galanin is a small, easily degraded peptide.
- **Benchmark:** if other lab assays agree well with Olink but galanin does not, the problem is specific to galanin; if all agree poorly, the samples themselves are the issue.
- The sample checks (`0_data/sample_checks.csv`) are flagged, not corrected. Samples failing Olink QC are left out (`drop_failed_samples` in the config).

---

## 7. Settings (`leip_galanin/config.yml`)

`galanin` (Olink assay, ELISA column), `covariates` (adjusted analyses), `hdl` (lipids, HDL
proteins, comparisons), `neuroendocrine_proteins`, `lab_vs_olink` (benchmark pairs: clinical
column → Olink assay), `expected_associations` (sanity checks), `clinical` (which variables are
tested), `gsea`, `bootstrap`, `permutations`, `seed`.

The folder uses three generic helpers of the repository: `R/utils.R` (settings, output files),
`R/qc.R` (Olink LOD, as in the main analysis) and `R/report.R` (PDF pages). Nothing else of the
dISF/serum analysis is used.
