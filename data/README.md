# Data folder

Put the input files here, with exactly these names (paths are set in `config.yml`):

```
data/
├── manifest.xlsx                          Olink sample submission sheet (sheet "manifest" is used)
├── LEIP_clinical_parameters_n35.xlsx      LEIP clinical data (sheets "Key_parameters" and "All_SORB_parameters")
├── Explore_HT_Fixed_LOD.csv               Olink fixed LOD file (Explore HT, version >= 6.0.0)
├── severity.xlsx                          optional: SubjectID, Visit, SCORAD / EASI / NRS ...
└── npx/
    ├── O-MicroAD_ISF_NPX_2026-09-24.parquet
    └── O-MicroAD_Serum_NPX_2026-09-24.parquet
```

All `.parquet` files in `data/npx/` are read. With one file per matrix, the file name must contain
"ISF" or "Serum"; step 02 checks that each file's samples have that matrix in the manifest. Controls
of plate 2 (in both files) are handled per file for the LOD and counted once for control QC.

The data files are git-ignored, so they are not uploaded to GitHub.
