# Data folder

Put the input files here (paths are set in `config.yml`; the usual Olink and lab file names are
also found when they differ slightly, e.g. `Explore HT_Fixed LOD.csv`):

```
data/
├── LEIP_clinical_parameters_n35.xlsx      LEIP clinical data: sheet "Key_parameters" (incl. "Galanin [pg/mL]")
│                                          and sheet "All_SORB_parameters". Its column Olink_SampleID defines
│                                          which Olink samples are LEIP samples.
├── Explore_HT_Fixed_LOD.csv               Olink fixed LOD file (Explore HT, version >= 6.0.0) - recommended
├── manifest.xlsx                          optional: Olink sample submission sheet, to cross-check the SubjectIDs
└── npx/
    └── O-MicroAD_Serum_NPX_2026-09-24.parquet   the Olink serum file that contains the LEIP samples
```

All `.parquet` files in `data/npx/` are looked at; files without LEIP samples (e.g. the dISF file of
the same Olink run) are skipped.

The data files are git-ignored, so they are not uploaded to GitHub.
