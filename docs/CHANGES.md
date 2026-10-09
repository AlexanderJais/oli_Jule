# What changed – and what to rerun

Newest first. After pulling a new version, read the entries since your last run. **"Rerun from"**
is the earliest step you need: `start_at <- N; source("run_all.R")`. If in doubt, rerun everything
with `source("run_all.R")`.

| Date | Change | Rerun from | New packages? |
|---|---|---|---|
| 2026-10-09 | **New: separate Leipzig biobank analysis**, `source("run_leipzig.R")` (not part of `run_all.R`) → `output/leipzig/`. Reads only the Leipzig serum samples. For LEIP_35 (no clinical data): which proteins differ from the other 34, checks (duplicate, sample handling), and what its proteins say about sex, age, BMI, body fat, lipids, CRP, kidney function … (estimate with 80 % range, low / middle / high, or "not predictable"). Published models in `reference/`. Also: QC overview PNG works on a Mac without XQuartz. | – (`run_leipzig.R`) | `patchwork` |
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
