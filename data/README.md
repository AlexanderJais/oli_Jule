# Data folder

Put the input files here. These are the names set in `config.yml`; if a file has a different name,
the usual original name is also found (e.g. `Explore HT_Fixed LOD.csv`) as long as only one candidate is in the folder.

```
data/
├── manifest.xlsx                          Olink sample submission sheet (sheet "manifest" is used)
├── LEIP_clinical_parameters_n35.xlsx      LEIP clinical data (sheet "Key_parameters")
├── Explore_HT_Fixed_LOD.csv               Olink fixed LOD file (Explore HT, version >= 6.0.0)
├── reference/proteinatlas.tsv.zip          optional: Human Protein Atlas table (public; tools/download_hpa.R)
├── severity.xlsx                          optional: SubjectID, Visit (V1 or 1), SCORAD / EASI / NRS ...
└── npx/
    ├── O-MicroAD_ISF_NPX_2026-09-24.parquet
    └── O-MicroAD_Serum_NPX_2026-09-24.parquet
```

All `.parquet` files in `data/npx/` are read. With one file per matrix, the file name must contain
"ISF" or "Serum"; step 02 checks that each file's samples have that matrix in the manifest. Controls
of plate 2 (in both files) are handled per file for the LOD and counted once for control QC.

The data files are git-ignored, so they are not uploaded to GitHub. Never commit them: the repository is public.
When asking for help, share error messages and file, column or sample names - not the data.
