# What changed – and what to rerun

Newest first. After pulling a new version, read the entries since your last run. **"Rerun from"**
is the earliest step you need: `start_at <- N; source("run_all.R")`. If in doubt, rerun everything
with `source("run_all.R")`.

| Date | Change | Rerun from | New packages? |
|---|---|---|---|
| 2026-10-13 | **Code audit (qc-overview) and README update.** Relapse screen: the serum volcano title now gives the serum patient numbers; the same-plate check also covers serum markers and lowered markers; GSEA also for V1 serum, with the gene sets changed in the same direction in dISF and serum (`V1_relapse_GSEA.xlsx`, sheet `gene_sets_dISF_and_serum`, and in the executive summary). README: table of the dedicated analyses at the top, KEGG download, new terms. | – (`run_relapse_v1.R`) | – |
| 2026-10-13 | **V1 relapse screen: serum and dISF + serum.** The same screen on the V1 serum of the same patients; markers in both matrices (sheet `dISF_and_serum`, with a chance check using the same random labels for both), `V1_relapse_dISF_vs_serum.pdf`, the shared markers over the visits, serum volcano, and serum / shared sections in the executive summary. | – (`run_relapse_v1.R`) | – |
| 2026-10-13 | **V1 relapse screen: executive summary PDF** `relapse_V1/V1_relapse_executive_summary.pdf` – the answer, the prediction test, the ranked predictors, the gene sets higher in relapsers (V1, V2) and the methods. | – (`run_relapse_v1.R`) | – |
| 2026-10-13 | **Code audit** of the new analyses: the relapse prediction test counts a left-out patient without a value as not predicted (it could stop the count); the dISF-only class table no longer uses `dplyr::pick()`; KEGG download keeps working if KEGG changes its release line; model workbook sheet names always fit Excel's 31 characters; the folder guide lists `dISF_only_proteins.xlsx`, `visit3/` and `relapse_V1/`. | – | – |
| 2026-10-13 | **New separate V1 relapse screen**, `source("run_relapse_v1.R")` → `output/relapse_V1/`: dISF proteins already higher at V1 in patients who later relapse, ranked as candidate targets with an exact relabelling check, leave-one-out, AUC, persistence at V2, serum and drug-target flags; volcano, per-patient course plots, GSEA. | – (`run_relapse_v1.R`) | – |
| 2026-10-11 | **New separate visit 3 analysis**, `source("run_visit3.R")` → `output/visit3/`: differential proteins in dISF (lesion site vs healthy, lesion site vs non-lesional) and serum (AD vs healthy) as Excel files, volcano plots with the 25 most regulated proteins labelled, GSEA (incl. KEGG), and dISF vs serum (proteins and gene sets seen only in dISF). | – (`run_visit3.R`) | `ggrepel` |
| 2026-10-11 | **Current KEGG** for the KEGG analysis of the dISF-only proteins: run `source("tools/download_kegg.R")` once (372 human pathways; the release date is in the files). Without it the older MSigDB KEGG collection (186 pathways) is used. | 02b (`source("run_qc_overview.R")`) | – |
| 2026-10-10 | **Which proteins are measurable in dISF only?** On matched dISF / serum samples (same person and visit, skin types pooled): `qc_overview/dISF_only_proteins.xlsx` (overall, per protein class, serum only) and `dISF_only_proteins.pdf`; their pathways (over-representation against all measured proteins) in sheet `pathways` and `dISF_only_pathways.pdf`, also within transcription factors, enzymes and intracellular proteins (`dISF_only_pathways_<class>.pdf`); the same with KEGG pathways in `dISF_only_KEGG.xlsx` / `dISF_only_KEGG*.pdf`. `run_qc_overview.R` runs only steps 01, 02 and 02b. | 02b (`source("run_qc_overview.R")`) | – |
| 2026-10-08 | **QC overview figure fixed:** Nimbus Sans now comes from the project's `fonts/` folder and is drawn as outlines (before, the PDF only referred to the font and other computers showed a substitute; now it looks the same everywhere, but the figure text can't be selected). Ring-chart labels keep an even distance from the rings; rings always form full circles. | 02b | `showtext` |
| 2026-10-08 | **New QC step 02b: proteins above / below LOD** → `qc_overview/`: ring charts (serum, dISF) and the distribution across 14 protein classes (Human Protein Atlas), Nimbus Sans, Helmholtz violet / pink. Also a page in the executive summary. `start_at` now uses script numbers (`start_at <- "02b"` works). For the classes run `source("tools/download_hpa.R")` once. | 02b | – (font Nimbus Sans optional) |
| 2026-10-04 | **Easier navigation.** `output/00_FOLDER_GUIDE.pdf` (one page: what is in each folder). Model results also as Excel with one sheet per volcano panel: `models/ISF_results.xlsx`, `Serum_results.xlsx`, `ISF_by_visit_results.xlsx`. Heatmaps without plate bar. | 04 | – |
| 2026-10-04 | **Third code audit.** Step 12 no longer shows CPUO serum as AD; step 16 serum correlations exclude CPUO; step 17 uses only proteins detected in MicroAD serum for serum effects and compares relapsers only at cleared visits; Q3 counts each relapse data set once; RELAD workbook robust to duplicate assay names; summary tables add `serum_MicroAD_AD_vs_HC`. `start_at` is now removed after a successful run, and the test no longer leaves its settings active. | 05 | – |
| 2026-10-04 | **Step 09 heatmaps.** The most-variable-proteins heatmap is now clustered, with colour bars for skin state, patient and visit. There is a new heatmap of the strongest lesional vs non-lesional proteins. The numbers behind both are in `isf_profile.xlsx` (sheets `variable_*`, `lesional_*`). | 09 | `pheatmap` |
| 2026-10-04 | **New step 16: proteins correlating with TNFRSF9 (CD137)** in dISF, including IL-33, IL-4, CSF2, IL6, IL18, CXCL8, IL1RL1, KIT, KITLG, TPSAB1 and FCER1A → `tnfrsf9_correlation/`. | 01 (full run) | – |
| 2026-10-04 | **New step 17: serum vs dISF signatures**, MicroAD only → `signatures/`. For the tissue origin, run `source("tools/download_hpa.R")` once. | 01 (full run) | – |
| 2026-10-04 | **RELAD/RELAD2 workbook** with all their serum results → `relad/RELAD_RELAD2_serum_results.xlsx` (step 05). dISF vs serum comparisons now use **MicroAD serum only** (new model `MicroAD_AD_vs_HC`). | 01 (full run) | – |
| 2026-10-04 | **Steps renumbered:** the summary is now step **18**, the export step **19**. | – | – |
| 2026-10-03 | Second code audit: date handling, pathway analysis on re-runs, severity file visit format, "regulated at all visits" rule, Q4 pre-relapse predictor. | 01 | `lmerTest` (listed explicitly) |
| 2026-10-03 | Executive summary review: evidence plots for the key questions, explanation of the mast cell score, full lists in `Executive_summary_tables.xlsx`, dISF vs serum enrichment per skin type. | 10 | – |
| 2026-10-02 | Manifest v4 (Sex, NoRELAD2, clinical states). CPA3, CMA1, KIT and HDC added to the mast cell markers. Key questions step added. | 01 | – |

**Still under review** (results available, interpret with care):
- step 14: the pooled relapse comparison includes lesional V1 samples;
- step 15: the wording and the choice of one primary test for key questions Q2–Q5.

## How to update

1. Get the new version. In RStudio: **Git** tab → **Pull**. In a terminal: `git pull`.
2. If the table lists new packages, run `source("install_packages.R")` again. It installs only missing
   packages.
3. Restart R (RStudio: Session → Restart R), so no old settings such as `start_at` remain.
4. Rerun from the step given in the table.
