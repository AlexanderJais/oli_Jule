# O-MicroAD – Olink Explore HT analysis

R pipeline for the Olink Explore HT run of O-MicroAD: dermal interstitial fluid (dISF) and serum
from the MicroAD study (atopic dermatitis patients followed over up to 6 visits, healthy volunteers,
CPUO), plus serum from RELAD / RELAD2 and LEIP biobank controls.

ISF and serum are analysed separately. They are only combined for the ISF–serum correlation
(step 06) and the LEIP population check (step 08).

## Quick start

```bash
Rscript install_packages.R        # once
# put the data files in data/ (see below), then
Rscript run_all.R                 # runs steps 01-08; results in output/
```

To run a single step, use `Rscript scripts/0X_....R`. Each step reads what the previous one saved.

## Data (never committed – `data/` and `output/` are git-ignored)

| file | what | where it is set |
|---|---|---|
| `data/manifest.xlsx` | Olink sample submission sheet; the `manifest` sheet is the master | `paths$manifest` |
| `data/npx/*.parquet` | Olink NPX file(s); several files (e.g. one per matrix) are combined | `paths$npx_dir` |
| `data/LEIP_clinical_parameters_n35.xlsx` | LEIP clinical data, sheet `Key_parameters` (optional) | `paths$leip_clinical` |
| `data/Explore_HT_Fixed_LOD.csv` | Olink fixed LOD file for Explore HT, version ≥ 6.0.0, from olink.com (recommended) | `paths$fixed_lod` |

All settings (thresholds, FDR, number of cores) are in `config.yml`.

## Steps

| script | does | main output |
|---|---|---|
| `01_metadata.R` | Cleans the manifest: harmonised groups (AD / HC / CPUO / Biobank), skin site and state, visit and days since V1, relapse visit and visits before relapse, sample volume, and LEIP clinical data. Checks IDs against the plate layout, and flags inconsistencies without dropping anything. | `metadata/sample_metadata.csv`, `metadata/data_flags.csv` |
| `02_import_qc.R` | Reads the parquet, checks sample IDs and normalisation, and works out LOD. Sample QC uses Olink flags plus a median/IQR outlier check. Keeps an assay in a matrix if it is detected in ≥ 50 % of at least one group, and reports control CVs and low-volume samples. | `qc/*`, `data/npx_clean.rds`, `data/npx_wide.rds` |
| `03_explore.R` | Per matrix: PCA coloured by state, group, visit, plate and cohort; variance partitioning. | `explore/*` |
| `04_isf_models.R` | ISF models per protein (see below). | `models/ISF_*` |
| `05_serum_models.R` | Serum models per protein (see below). | `models/Serum_*` |
| `06_isf_vs_serum.R` | On matched MicroAD visits: within-subject (repeated-measures) and between-subject correlation of ISF and serum, separately for the lesional and non-lesional site. | `isf_serum/*` |
| `07_enrichment.R` | GSEA (fgsea) for every contrast: MSigDB Hallmark, Reactome, GO:BP, plus a custom AD/Th2 set. | `enrichment/*` |
| `08_leip_reference.R` | For the proteins significant in step 06: LEIP normal range, where AD patients fall in it, clinical associations in LEIP, detectability, and LEIP vs in-study controls. | `leip_reference/*` (incl. `.xlsx`) |

### Models

The models are limma/dream (variancePartition) with empirical Bayes moderation. Subject is a random
effect wherever a person contributes several samples. NPX differences are on the log2 scale, and
BH FDR is applied within each contrast.

**ISF (04)**, all models adjusted for plate:
- `states_all_visits`: lesional, ex-lesional and non-lesional AD skin vs each other and vs healthy skin; CPUO lesional vs non-lesional.
- `baseline_V1`: the same comparisons at V1 only.
- `time_ex_lesional`, `time_non_lesional`: change per week in AD skin.
- `relapse_ex_lesional`, exploratory: cleared skin of relapsers vs non-relapsers.
- `relapse_delta_xL_minus_NL`, exploratory: the same question on the within-visit ex-lesional minus non-lesional difference. Plate and systemic day-to-day variation cancel out.

**Serum (05):**
- AD vs controls, fitted twice: against the in-study healthy controls (adjusted for cohort and plate) and against LEIP biobank controls (adjusted for plate).
  - `Serum_AD_vs_controls_agreement.csv` marks proteins that agree in both comparisons.
  - `HC_vs_Biobank` shows proteins affected by the biobank source.
- `MicroAD_active_vs_cleared`: serum when the tracked lesion is active vs cleared.
- `MicroAD_relapse`: relapsers vs non-relapsers at cleared visits.
- `RELAD_relapse`: RELAD and RELAD2, adjusted for cohort and plate, plus a sensitivity analysis without the samples with conflicting relapse labels.

## Design decisions

- **`PCNormalizedNPX` is analysed, not intensity-normalised NPX.** Plate 1 is ISF only, plate 2 is mixed, and plates 3–4 are serum. Intensity normalisation assumes randomised samples of one matrix and would distort plate 2. Step 02 reports the `Normalization` column of the delivered file.
- **LOD:** OlinkAnalyze's negative-control LOD requires ≥ 10 negative controls, and this run has 8 (2 per plate).
  - Preferred: Olink's fixed LOD file (`paths$fixed_lod`).
  - Fallback: the same formula (median + max(0.2, 3 SD)) on the 8 negative controls pooled.
  - `qc/lod.csv` records which source was used for each assay.
- **Relapse in ISF is exploratory.** All four ISF relapsers are on plate 1, and plate 2 holds only non-relapsers.
- **Metadata issues are flagged, not fixed**, in `metadata/data_flags.csv`. As of manifest v3 this covers the RELAD / RELAD2 label conflicts, 5 low-volume ISF samples, and LEIP_35 without clinical data.
- **Age and sex** are currently only available for LEIP, so the serum models are not adjusted for them.

## Testing without real data

```bash
Rscript tests/test_pipeline.R
```

This simulates an Explore HT parquet with the real sample layout and known effects, runs all
steps with `data_sim/config_sim.yml` (output in `output_sim/`), and checks that the effects are
recovered and false positives stay rare.
