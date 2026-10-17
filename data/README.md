# Data folder

Put the input files here, under the names set in `config.yml`. For the fixed LOD file, the LEIP data and the severity
scores the original names also work (e.g. `Explore HT_Fixed LOD.csv`), as long as only one matching file is in the
folder. The manifest must have exactly the name in `config.yml`.

```
data/
├── 2026-10-15_manifest.xlsx               required: sample manifest of 2026-10-15 (exactly this name, see below)
├── npx/                                   required: the Olink NPX files
│   ├── O-MicroAD_ISF_NPX_2026-09-24.parquet
│   └── O-MicroAD_Serum_NPX_2026-09-24.parquet
├── Explore_HT_Fixed_LOD.csv               recommended: Olink fixed LOD file (Explore HT, version >= 6.0.0)
├── LEIP_clinical_parameters_n35.xlsx      optional: LEIP clinical data (sheet "Key_parameters")
├── severity.xlsx                          optional: SubjectID, Visit (V1–V6 or 1–6), SCORAD / EASI / NRS ...
└── reference/
    ├── proteinatlas.tsv.zip               optional: Human Protein Atlas table (public; tools/download_hpa.R) - protein classes (02b), tissue origin (17)
    └── kegg_hsa_pathways.tsv              optional: current KEGG human pathways (tools/download_kegg.R) - KEGG analysis (02b)
```

What each file needs:

- **`2026-10-15_manifest.xlsx`** (the name set in `config.yml` → `paths$manifest`; it is never replaced by another
  file such as an older `manifest.xlsx` – if it is missing, the run stops):
  - the master table, one row per sample: sheet `manifest`, or else the one sheet with these columns (e.g. `Sheet1`):
    `SampleID`, `SubjectID`, `Visit`, `Skin`, `SampleType` (`dISF` or `Serum`), `Group`, `Study`, `Date`, `SampleName`,
    `Relapse`, `TimeToRelapse`, `ClinicalStateSkin`, `plate`, `well` and `Note`. `Sex`, `NoRELAD2` and `ClinicalPhase`
    are used when present; columns without a header are ignored.
  - optional `Sample Info`: sample volumes (a header row starting with `SampleID`, and a `SampleVolume` column).
  - optional `Plate Layout`: the plate maps; step 01 then checks that each sample sits in the same well as in the table.
- **`npx/`**: all `.parquet` files in this folder are read. With one file per matrix, the file name
  should contain "ISF" or "Serum", so that step 02 can check that each file's samples have that
  matrix in the manifest. Plate 2 holds both matrices, so its control wells appear in both files.
  This is expected: the LOD is worked out per file, and the controls are counted once in the control QC.
- **Fixed LOD file**: from olink.com. Without it, the LOD comes from the negative controls, which is
  less precise.
- **LEIP clinical data**: sheet `Key_parameters`, one row per LEIP sample, with the columns
  `Olink_SampleID` (matching `SampleID` in the manifest) and `sex_MF`. Further columns such as `age` and
  `BMI` are used in steps 08 and 12.
- **Severity scores** (`.xlsx` or `.csv`): one row per patient and visit, with `SubjectID`, `Visit`
  and one numeric column per score. Steps 11 and 12 then model the proteins against each score.
- **Human Protein Atlas table**: download it once with `source("tools/download_hpa.R")`, which
  creates `data/reference/`. Step 02b uses it for the protein classes, step 17 for the tissue origin of proteins.

The data files are git-ignored, so they are not uploaded to GitHub. Never commit them: the repository is public.
When asking for help, share error messages and file, column or sample names – not the data.
