# Data folder

Put the input files here, with exactly these names (paths are set in `config.yml`):

```
data/
├── manifest.xlsx                          Olink sample submission sheet (sheet "manifest" is used)
├── LEIP_clinical_parameters_n35.xlsx      LEIP clinical data (sheet "Key_parameters")
├── Explore_HT_Fixed_LOD.csv               Olink fixed LOD file (Explore HT, version >= 6.0.0)
├── severity.xlsx                          optional: SubjectID, Visit, SCORAD / EASI / NRS ...
└── npx/
    └── <your Olink NPX file>.parquet      one or more .parquet files from the Olink delivery
```

The data files are git-ignored, so they are not uploaded to GitHub.
