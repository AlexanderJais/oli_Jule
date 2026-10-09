# Reference models and marker lists (Leipzig analysis, `run_leipzig.R`)

Public, published information only – no study data. Used by `leipzig/02_case.R` to estimate the
clinical values of one Leipzig sample from its proteins. Every model is fixed here in advance; the pipeline
only re-scales it (offset and slope, plus sex where listed) on the Leipzig samples with clinical data.

| File | What | Source | Licence |
|---|---|---|---|
| `age_clock_goeminne2025.csv` | Protein age clock: weight per protein (mean of the 5 published cross-validation fits) and its SD in UK Biobank. Score = Σ weight × UKB SD × z. Three names changed to the Explore HT names: NTproBNP → NPPB, CERT → CERT1, WARS → WARS1. | Goeminne LJE et al. (2025) *Plasma protein-based organ-specific aging and mortality models unveil diseases as accelerated aging of organismal systems.* Cell Metabolism 37:205–222. doi:10.1016/j.cmet.2024.10.005; github.com/ludgergoeminne/organAging, commit 5147b03 (`Conventional_coefs_GTEx_4x_FC.csv`, `standard_deviations.rds`) | Academic / non-commercial use only, with citation – see `LICENSE_age_clock_goeminne2025.txt` |
| `bmi_score_watanabe2023.csv` | Protein BMI score (67 proteins; LASSO weights, mean of 10 models) | Watanabe K et al. (2023) *Multiomic signatures of body mass index identify heterogeneous health phenotypes and responses to a lifestyle intervention.* Nature Medicine 29:996–1008 (Supplementary Data 3); 67-protein Olink Explore version as in Wang et al. (2024) Diabetes | cite the papers |
| `leip_case_markers.csv` | Marker sets with weights: sex, lipids, CRP, kidney, insulin resistance, body fat, apolipoproteins, age (fallback), and the sample-handling checks. `source` says where each set comes from; UK Biobank weights are partial correlations from the UK Biobank Olink proteome–phenome atlas (Deng et al. 2025, Cell). | see column `source` | cite the papers |

Not on the Explore HT panel, so not used: PSA (KLK3), CGA, CRP, SAA1/2, insulin (only C-peptide), adiponectin.
Olink's APOB assay is not used for LDL (Spearman 0.08 with clinical ApoB; Sun et al. 2023, Nature).
