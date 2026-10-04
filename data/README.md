# Data folder

Put the input files here, under the names set in `config.yml`. The original names from Olink or the
lab also work (e.g. `Explore HT_Fixed LOD.csv`, `... Sample Submission Sheet.xlsx`), as long as only
one matching file is in the folder.

```
data/
├── manifest.xlsx                          required: Olink sample submission sheet (3 sheets, see below)
├── npx/                                   required: the Olink NPX files
│   ├── O-MicroAD_ISF_NPX_2026-09-24.parquet
│   └── O-MicroAD_Serum_NPX_2026-09-24.parquet
├── Explore_HT_Fixed_LOD.csv               recommended: Olink fixed LOD file (Explore HT, version >= 6.0.0)
├── LEIP_clinical_parameters_n35.xlsx      optional: LEIP clinical data (sheet "Key_parameters")
├── severity.xlsx                          optional: SubjectID, Visit (V1–V6 or 1–6), SCORAD / EASI / NRS ...
└── reference/
    └── proteinatlas.tsv.zip               optional: Human Protein Atlas table (public; tools/download_hpa.R)
```

What each file needs:

- **`manifest.xlsx`** needs three sheets; step 01 stops with an error if one of them is missing.
  - `manifest`: the master table, one row per sample, with the columns `SampleID`, `SubjectID`,
    `Visit`, `Skin`, `SampleType` (`dISF` or `Serum`), `Group`, `Study`, `Date`, `SampleName`,
    `Relapse`, `TimeToRelapse`, `ClinicalStateSkin`, `plate`, `well` and `Note`. `Sex` and
    `NoRELAD2` are used when present.
  - `Sample Info`: sample volumes (a header row starting with `SampleID`, and a `SampleVolume` column).
  - `Plate Layout`: the plate maps; step 01 checks that each sample sits in the same well as in `manifest`.
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
  creates `data/reference/`. Step 17 uses it for the tissue origin of proteins.

The data files are git-ignored, so they are not uploaded to GitHub. Never commit them: the repository is public.
When asking for help, share error messages and file, column or sample names – not the data.
