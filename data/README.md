# Data folder

Put the input files here, with exactly these names (paths are set in `config.yml`):

```
data/
├── manifest.xlsx                          Olink sample submission sheet (sheet "manifest" is used; manifest Ver2)
├── hpa_annotation.tsv                     optional: Human Protein Atlas annotation (step 20 downloads it if missing)
├── LEIP_clinical_parameters_n35.xlsx      LEIP clinical data (sheet "Key_parameters")
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

The manifest is read with the codes of manifest **Ver2**: column `Sex` (MicroAD), `NoRELAD2`
(ignored), `Relapse = activeAD`, `TimeToRelapse = relapse_<1w / relapse_>1w / non-relapse` (RELAD2) or a
number (RELAD), and the spelling variants `remisison` / `helthy` in `ClinicalStateSkin`, which are
harmonised to `remission` / `healthy`.
