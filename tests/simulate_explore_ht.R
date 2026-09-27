# Simulate an Olink Explore HT parquet file for the O-MicroAD sample layout.
#
# Uses the synthetic manifest from tests/make_synthetic_manifest.R (same layout and design as
# the real one, no study data) so the whole pipeline can be run end to end without real data. Known effects are built in and
# written to data_sim/truth.csv, so tests/test_pipeline.R can check they are recovered.
#
# Run from the repository root:  Rscript tests/simulate_explore_ht.R

source("R/utils.R")
source("R/metadata.R")
set.seed(20260926)

# synthetic manifest + LEIP file (tests/make_synthetic_manifest.R) - no study data needed
cfg  <- load_config("config.yml")
cfg$paths$manifest <- "data_sim/manifest.xlsx"
cfg$paths$leip_clinical <- "data_sim/LEIP_clinical.xlsx"
if (!file.exists(cfg$paths$manifest)) stop("Run tests/make_synthetic_manifest.R first.")
sim_dir <- "data_sim"
dir.create(file.path(sim_dir, "npx"), recursive = TRUE, showWarnings = FALSE)

meta   <- build_metadata(cfg$paths$manifest, cfg$paths$leip_clinical)
layout <- read_plate_layout(cfg$paths$manifest)

# ---- assays -----------------------------------------------------------------------
n_assays <- as.integer(Sys.getenv("SIM_N_ASSAYS", "400"))
ht <- readRDS(system.file("extdata", "OlinkID_HT_mapping.rds", package = "OlinkAnalyze")) |>
  distinct(OlinkID, .keep_all = TRUE)
# put well-known AD / Th2 proteins first so enrichment of the lesional signal is meaningful
th2 <- c("CCL17", "CCL22", "CCL18", "IL13", "IL4", "POSTN", "CCL26", "IL5", "TSLP", "IL31",
         "CCL11", "CCL13", "CCL24", "IL4R", "IL13RA2", "MMP12", "PI3", "SERPINB4", "S100A7", "S100A8",
         "S100A9", "IL19", "IL22", "IL36G", "TNFRSF9")   # TNFRSF9 = CD137 (focus protein)
focus_sim <- c("TNFSF9", "KITLG", "CPA4", "FCER1A", "TPSAB1", "TPSD1", "KIT")   # other focus proteins (step 12)
ht <- bind_rows(ht |> filter(Gene %in% th2), ht |> filter(Gene %in% focus_sim),
                ht |> filter(!Gene %in% c(th2, focus_sim)) |> slice_sample(prop = 1))
assays <- ht |> slice_head(n = n_assays) |>
  transmute(OlinkID, UniProt, Assay = Gene, Block = as.character(Block), idx = row_number())

eff <- function(from, to) assays$idx >= from & assays$idx <= to
role <- case_when(
  eff(1, 25)    ~ "ISF_lesional",        # lesional > ex-lesional > non-lesional
  eff(26, 35)   ~ "ISF_AD_vs_HC",        # all AD skin (also non-lesional) > healthy
  eff(36, 45)   ~ "ISF_relapse",         # ex-lesional skin of relapsers > non-relapsers
  eff(46, 60)   ~ "ISF_serum_coupled",   # shared subject-visit signal in ISF and serum
  eff(61, 75)   ~ "Serum_AD_vs_HC",
  eff(76, 85)   ~ "Serum_relapse",       # RELAD/RELAD2 relapse
  eff(86, 95)   ~ "Biobank_shift",       # LEIP pre-analytical shift
  eff(96, 100)  ~ "ISF_lesion_restricted", # below LOD in ISF except in lesional skin
  eff(n_assays - 19, n_assays) ~ "ISF_undetected",
  TRUE ~ "null"
)
assays$role <- role
# coupled proteins 46-50 also depend on BMI in the population
assays$bmi_slope <- if_else(eff(46, 50), 0.12, 0)
# LEIP (step 18): Olink GAL follows the galanin ELISA and shares a factor with 3 partner proteins
# (positions 114-116); FABP4 follows the lab A-FABP; IL10 is not linked to the lab IL-10
leip_sim <- c("GAL", "FABP4", "IL10")
if (n_assays >= 130) {
  for (k in seq_along(leip_sim)) if (!leip_sim[k] %in% assays$Assay) {
    m <- ht |> filter(Gene == leip_sim[k]) |> slice(1)
    if (nrow(m)) assays[110 + k, c("OlinkID", "UniProt", "Assay")] <- list(m$OlinkID, m$UniProt, m$Gene)
  }
  assays$role[assays$Assay %in% c("GAL", "FABP4") & assays$role == "null"] <- "LEIP_lab_linked"
  assays$role[114:116][assays$role[114:116] == "null"] <- "LEIP_GAL_partner"
}

base   <- rnorm(n_assays, 6, 1.5)          # serum level
isf_off <- rnorm(n_assays, -1.2, 0.8)      # ISF relative to serum
isf_off[assays$role == "ISF_undetected"] <- -6
isf_off[assays$role == "ISF_lesion_restricted"] <- -6.5
isf_off[eff(101, 110)] <- 2        # clearly enriched in ISF relative to serum
nc_lvl <- base - 5.5                       # negative control background

# ---- samples and controls ------------------------------------------------------------
ctrl <- layout |>
  filter(content %in% c("PC", "SC", "Neg Ctrl")) |>
  mutate(SampleType = recode(content, PC = "PLATE_CONTROL", SC = "SAMPLE_CONTROL",
                             `Neg Ctrl` = "NEGATIVE_CONTROL"),
         SampleID = recode(content, `Neg Ctrl` = "NC")) |>
  group_by(plate, SampleID) |> mutate(SampleID = paste0(SampleID, row_number())) |> ungroup()   # PC1..PC5 on every plate, like Olink
samples <- bind_rows(
  meta |> transmute(SampleID, SampleType = "SAMPLE", plate, well),
  ctrl |> select(SampleID, SampleType, plate, well)
) |> left_join(meta |> select(-plate, -well), by = "SampleID")

# latent terms
subj <- unique(na.omit(samples$SubjectID))
# person-level effects, independent between ISF and serum (only the coupled proteins share signal)
u_subj <- list(ISF = matrix(rnorm(length(subj) * n_assays, 0, 0.3), length(subj), dimnames = list(subj, NULL)),
               Serum = matrix(rnorm(length(subj) * n_assays, 0, 0.3), length(subj), dimnames = list(subj, NULL)))
sv_key <- with(samples, ifelse(is.na(visit), SubjectID, paste(SubjectID, visit)))
sv <- unique(na.omit(sv_key))
z_sv <- matrix(rnorm(length(sv) * n_assays, 0, 0.9), length(sv), dimnames = list(sv, NULL))
plate_off <- matrix(rnorm(4 * n_assays, 0, 0.08), 4, dimnames = list(paste("Plate", 1:4), NULL))

npx <- matrix(NA_real_, nrow(samples), n_assays)
for (i in seq_len(nrow(samples))) {
  s <- samples[i, ]
  if (s$SampleType == "NEGATIVE_CONTROL") { npx[i, ] <- nc_lvl + rnorm(n_assays, 0, 0.3); next }
  if (s$SampleType == "PLATE_CONTROL")    { npx[i, ] <- base + rnorm(n_assays, 0, 0.08); next }
  if (s$SampleType == "SAMPLE_CONTROL")   { npx[i, ] <- base + 0.2 + rnorm(n_assays, 0, 0.12); next }
  isf <- s$matrix == "ISF"
  x <- base + if (isf) isf_off else 0
  x <- x + u_subj[[s$matrix]][s$SubjectID, ] + plate_off[s$plate, ]
  if (isf && s$group == "AD") {
    x <- x + 0.8 * (assays$role == "ISF_AD_vs_HC")
    if (identical(s$state, "lesional"))    x <- x + 1.5 * (assays$role == "ISF_lesional") +
                                              5 * (assays$role == "ISF_lesion_restricted")
    if (identical(s$state, "ex-lesional")) {
      x <- x + 0.6 * (assays$role == "ISF_lesional")
      if (identical(s$relapse, "relapse")) x <- x + 1.0 * (assays$role == "ISF_relapse")
    }
  }
  if (!is.na(s$visit) && s$cohort == "MicroAD") x <- x + z_sv[sv_key[i], ] * (assays$role == "ISF_serum_coupled")
  if (!isf) {
    if (s$group == "AD") x <- x + 0.7 * (assays$role == "Serum_AD_vs_HC")
    if (s$cohort %in% c("RELAD", "RELAD2") && identical(s$relapse, "relapse"))
      x <- x + 0.8 * (assays$role == "Serum_relapse")
    if (s$cohort == "LEIP") {
      x <- x + 0.6 * (assays$role == "Biobank_shift")
      if (!is.na(s$BMI)) x <- x + assays$bmi_slope * (s$BMI - 24)
    }
  }
  npx[i, ] <- x + rnorm(n_assays, 0, 0.35)
}
# LEIP lab-linked proteins (step 18), drawn from a separate random stream so that the rest of the
# simulation stays the same
rng <- .Random.seed; set.seed(18)
leip_i <- which(samples$cohort %in% "LEIP")
z_gal <- rnorm(length(leip_i))
sorb <- readxl::read_excel(cfg$paths$leip_clinical, sheet = "All_SORB_parameters")
lab_effect <- \(v, centre, slope) slope * coalesce(log2(v / centre), 0)
for (a in which(assays$role == "LEIP_lab_linked")) {
  npx[leip_i, a] <- npx[leip_i, a] + switch(assays$Assay[a],
    GAL   = lab_effect(samples[["Galanin [pg/mL]"]][leip_i], 140, 1.8) + 0.8 * z_gal,
    FABP4 = lab_effect(sorb$c_AFABP4[match(samples$SampleID[leip_i], sorb$Olink_SampleID)], 12, 1.0))
}
for (a in which(assays$role == "LEIP_GAL_partner")) npx[leip_i, a] <- npx[leip_i, a] + 1.0 * z_gal
assign(".Random.seed", rng, envir = globalenv())
npx <- pmax(npx, nc_lvl[col(npx)] - 0.5)     # floor near background, like real data

# ---- long parquet in Explore HT format ------------------------------------------------
dar_id <- "SIM_DAR_1"
# Olink data model: ExtNPX = log2(Count / extension-control count); PC-normalised NPX = ExtNPX minus
# the median ExtNPX of the plate controls (per plate and assay). Count-based LODs depend on this.
blocks <- sort(unique(assays$Block))
ext_count <- expand_grid(i = seq_len(nrow(samples)), Block = blocks) |>
  mutate(ExtCount = as.integer(round(5000 * 2^rnorm(n(), 0, 0.1))))
long <- expand_grid(i = seq_len(nrow(samples)), a = seq_len(n_assays)) |>
  mutate(
    SampleID = samples$SampleID[i], SampleType = samples$SampleType[i],
    WellID = samples$well[i], PlateID = str_replace(samples$plate[i], "Plate ", "Plate"),
    DataAnalysisRefID = dar_id,
    OlinkID = assays$OlinkID[a], UniProt = assays$UniProt[a], Assay = assays$Assay[a],
    AssayType = "assay", Panel = "Explore_HT", Block = assays$Block[a],
    ExtNPX = npx[cbind(i, a)] - 4
  ) |>
  left_join(ext_count, by = c("i", "Block")) |>
  group_by(PlateID, OlinkID) |>
  mutate(PCNormalizedNPX = ExtNPX - median(ExtNPX[SampleType == "PLATE_CONTROL"])) |>
  ungroup() |>
  mutate(NPX = PCNormalizedNPX, Count = as.integer(round(ExtCount * 2^ExtNPX)),
         Normalization = "Plate control", AssayQC = "PASS", SampleQC = "PASS", ExploreVersion = "SIM") |>
  select(-a)

# internal control assays, one set per block (removed by the pipeline); ext_ctrl carries ExtCount
ctrl_assays <- expand_grid(i = seq_len(nrow(samples)), Block = blocks, ct = c("ext_ctrl", "inc_ctrl", "amp_ctrl")) |>
  left_join(ext_count, by = c("i", "Block")) |>
  mutate(SampleID = samples$SampleID[i], SampleType = samples$SampleType[i], WellID = samples$well[i],
         PlateID = str_replace(samples$plate[i], "Plate ", "Plate"), DataAnalysisRefID = dar_id,
         OlinkID = sprintf("OID9%d%03d", match(ct, c("ext_ctrl", "inc_ctrl", "amp_ctrl")), as.integer(Block)),
         UniProt = NA_character_, Assay = paste0(ct, "_", Block), AssayType = ct, Panel = "Explore_HT",
         ExtNPX = if_else(ct == "ext_ctrl", 0, rnorm(n(), 0, 0.2)), NPX = ExtNPX, PCNormalizedNPX = ExtNPX,
         Count = if_else(ct == "ext_ctrl", ExtCount, 5000L), Normalization = "Plate control",
         AssayQC = "PASS", SampleQC = "PASS", ExploreVersion = "SIM") |>
  select(-ct)
long <- bind_rows(long, ctrl_assays) |> select(-i, -ExtCount)

# a few QC flags to exercise the QC code
bad <- meta$SampleID[meta$matrix == "ISF"][c(3, 40)]
long$SampleQC[long$SampleID == bad[1]] <- "FAIL"
long$NPX[long$SampleID == bad[1]] <- long$NPX[long$SampleID == bad[1]] - 3
long$PCNormalizedNPX[long$SampleID == bad[1]] <- long$NPX[long$SampleID == bad[1]]
long$SampleQC[long$SampleID == bad[2]] <- "WARN"
long$AssayQC[long$OlinkID == assays$OlinkID[200] & long$PlateID == "Plate2"] <- "WARN"

# like the real delivery: intensity-normalised NPX (PCNormalizedNPX kept), and two assays that
# Olink excluded (Normalization = EXCLUDED, NPX empty for all samples)
long <- long |>
  group_by(PlateID, OlinkID) |>
  mutate(NPX = if_else(AssayType == "assay", PCNormalizedNPX - median(PCNormalizedNPX[SampleType == "SAMPLE"]), NPX)) |>
  ungroup() |>
  mutate(Normalization = "Intensity")
excl <- assays$OlinkID[c(n_assays - 1, n_assays)]
long <- long |> mutate(excluded = OlinkID %in% excl,
                       NPX = if_else(excluded, NA_real_, NPX), PCNormalizedNPX = if_else(excluded, NA_real_, PCNormalizedNPX),
                       Normalization = if_else(excluded, "EXCLUDED", Normalization)) |> select(-excluded)

# one file per matrix, like the Olink delivery: each file has its samples plus all control wells
# of the plates they sit on (so the mixed plate's controls appear in both files)
unlink(list.files(file.path(sim_dir, "npx"), full.names = TRUE))
smx <- meta |> select(SampleID, matrix)
for (mx in c("ISF", "Serum")) {
  plates_mx <- str_replace(unique(meta$plate[meta$matrix == mx]), "Plate ", "Plate")
  part <- long |> left_join(smx, by = "SampleID") |>
    filter(matrix %in% mx | (SampleType != "SAMPLE" & PlateID %in% plates_mx)) |> select(-matrix)
  tbl <- arrow::arrow_table(part)
  tbl$metadata <- list(FileVersion = "NA", ProjectName = "O-MicroAD", SampleMatrix = mx,
                       DataFileType = "NPX File", Product = "ExploreHT")
  arrow::write_parquet(tbl, file.path(sim_dir, "npx", sprintf("O-MicroAD_%s_NPX_SIM.parquet", mx)))
}

# Olink-style fixed LOD file
# same layout as Olink's file; every 5th assay uses a count-based LOD (LODCount on the
# extension-count scale, converted per sample by OlinkAnalyze)
lod <- assays |>
  transmute(OlinkID, AssayType = "assay", UniProt, Assay, Panel = "Explore_HT", Block,
            DataAnalysisRefID = dar_id, BimodalDistribution = FALSE,
            LODNPX = nc_lvl + 0.9 - base,
            LODCount = as.integer(pmax(150, round(5000 * 2^(nc_lvl + 0.9 - 4)))),
            LODMethod = if_else(idx %% 5 == 0, "lod_count", "lod_npx"), Version = "10.2.0")
write.table(lod, file.path(sim_dir, "fixed_lod.csv"), sep = ";", row.names = FALSE, quote = FALSE)

write_csv(assays |> mutate(isf_offset = isf_off) |> select(OlinkID, Assay, role, bmi_slope, isf_offset),
          file.path(sim_dir, "truth.csv"))

# clinical severity per AD visit: high when the tracked lesion is active
sev <- meta |> filter(cohort == "MicroAD", group == "AD", matrix == "ISF", site == "L") |>
  transmute(SubjectID, Visit = visit, EASI = round(pmax(0, if_else(state == "lesional", 12, 3) + rnorm(n(), 0, 2)), 1))
write_csv(sev, file.path(sim_dir, "severity.csv"))

# step 18 settings: lab/Olink pairs and sanity checks that exist in the simulation (+ one that does not)
leip_cfg <- cfg$leip_biobank
leip_cfg$lab_vs_olink <- list(galanin_elisa = "GAL", c_AFABP4 = "FABP4", IL10 = "IL10", c_CRP = "CRP")
leip_cfg$expected_associations <- c(map(assays$Assay[assays$bmi_slope > 0], \(a) list(protein = a, parameter = "BMI", direction = "positive")),
                                    list(list(protein = "NOT_ON_PANEL", parameter = "BMI", direction = "positive")))
yaml::write_yaml(list(
  paths = list(manifest = cfg$paths$manifest, npx_dir = file.path(sim_dir, "npx"),
               severity = file.path(sim_dir, "severity.csv"),
               leip_clinical = cfg$paths$leip_clinical, fixed_lod = file.path(sim_dir, "fixed_lod.csv"),
               output = "output_sim"),
  npx_column = cfg$npx_column, qc = cfg$qc, stats = cfg$stats, enrichment = cfg$enrichment,
  focus_proteins = cfg$focus_proteins, leip_biobank = leip_cfg
), file.path(sim_dir, "config_sim.yml"))

msg("Simulated %d samples/controls x %d assays -> %s", nrow(samples), n_assays, sim_dir)
