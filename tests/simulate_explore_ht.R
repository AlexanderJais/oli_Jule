# Simulate an Olink Explore HT parquet file for the O-MicroAD sample layout.
#
# Uses the real manifest (sample IDs, plates, wells, control wells) so the whole pipeline
# can be run end to end before the real data arrive. Known effects are built in and
# written to data_sim/truth.csv, so tests/test_pipeline.R can check they are recovered.
#
# Run from the repository root:  Rscript tests/simulate_explore_ht.R

source("R/utils.R")
source("R/metadata.R")
set.seed(20260926)

cfg  <- load_config("config.yml")       # real manifest + LEIP file
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
         "S100A9", "IL19", "IL22", "IL36G", "CXCL10")
ht <- bind_rows(ht |> filter(Gene %in% th2), ht |> filter(!Gene %in% th2) |> slice_sample(prop = 1))
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
         SampleID = paste(recode(content, `Neg Ctrl` = "NC"), str_remove(plate, "Plate "), well, sep = "_"))
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
npx <- pmax(npx, nc_lvl[col(npx)] - 0.5)     # floor near background, like real data

# ---- long parquet in Explore HT format ------------------------------------------------
dar_id <- "SIM_DAR_1"
long <- expand_grid(i = seq_len(nrow(samples)), a = seq_len(n_assays)) |>
  mutate(
    SampleID = samples$SampleID[i], SampleType = samples$SampleType[i],
    WellID = samples$well[i], PlateID = str_replace(samples$plate[i], "Plate ", "Plate"),
    DataAnalysisRefID = dar_id,
    OlinkID = assays$OlinkID[a], UniProt = assays$UniProt[a], Assay = assays$Assay[a],
    AssayType = "assay", Panel = "Explore_HT", Block = assays$Block[a],
    NPX = npx[cbind(i, a)],
    Count = as.integer(round(pmax(0, 2^(NPX - nc_lvl[a]) * 60 + rnorm(n(), 0, 10)))),
    ExtNPX = NPX - 1, Normalization = "Plate control", PCNormalizedNPX = NPX,
    AssayQC = "PASS", SampleQC = "PASS", ExploreVersion = "SIM"
  ) |> select(-i, -a)

# internal control assays (removed by the pipeline)
ctrl_assays <- expand_grid(SampleID = samples$SampleID, ct = c("ext_ctrl", "inc_ctrl", "amp_ctrl")) |>
  left_join(long |> distinct(SampleID, SampleType, WellID, PlateID), by = "SampleID") |>
  mutate(DataAnalysisRefID = dar_id, OlinkID = paste0("OID9", match(ct, c("ext_ctrl", "inc_ctrl", "amp_ctrl")), "000"),
         UniProt = NA_character_, Assay = paste0(ct, "_1"), AssayType = ct, Panel = "Explore_HT", Block = "1",
         NPX = rnorm(n(), 0, 0.2), Count = 5000L, ExtNPX = NPX, Normalization = "Plate control",
         PCNormalizedNPX = NPX, AssayQC = "PASS", SampleQC = "PASS", ExploreVersion = "SIM") |>
  select(-ct)
long <- bind_rows(long, ctrl_assays)

# a few QC flags to exercise the QC code
bad <- meta$SampleID[meta$matrix == "ISF"][c(3, 40)]
long$SampleQC[long$SampleID == bad[1]] <- "FAIL"
long$NPX[long$SampleID == bad[1]] <- long$NPX[long$SampleID == bad[1]] - 3
long$PCNormalizedNPX[long$SampleID == bad[1]] <- long$NPX[long$SampleID == bad[1]]
long$SampleQC[long$SampleID == bad[2]] <- "WARN"
long$AssayQC[long$OlinkID == assays$OlinkID[200] & long$PlateID == "Plate2"] <- "WARN"

tbl <- arrow::arrow_table(long)
tbl$metadata <- list(FileVersion = "NA", ProjectName = "SIM", SampleMatrix = "NA",
                     DataFileType = "NPX File", Product = "ExploreHT")
arrow::write_parquet(tbl, file.path(sim_dir, "npx", "sim_explore_ht.parquet"))

# Olink-style fixed LOD file
lod <- assays |>
  transmute(OlinkID, DataAnalysisRefID = dar_id, LODNPX = nc_lvl + 0.9, LODCount = 150L,
            LODMethod = "lod_npx", Panel = "Explore_HT", Version = "6.0.0")
write.table(lod, file.path(sim_dir, "fixed_lod.csv"), sep = ";", row.names = FALSE, quote = FALSE)

write_csv(assays |> mutate(isf_offset = isf_off) |> select(OlinkID, Assay, role, bmi_slope, isf_offset),
          file.path(sim_dir, "truth.csv"))

# clinical severity per AD visit: high when the tracked lesion is active
sev <- meta |> filter(cohort == "MicroAD", group == "AD", matrix == "ISF", site == "L") |>
  transmute(SubjectID, Visit = visit, EASI = round(pmax(0, if_else(state == "lesional", 12, 3) + rnorm(n(), 0, 2)), 1))
write_csv(sev, file.path(sim_dir, "severity.csv"))

yaml::write_yaml(list(
  paths = list(manifest = cfg$paths$manifest, npx_dir = file.path(sim_dir, "npx"),
               severity = file.path(sim_dir, "severity.csv"),
               leip_clinical = cfg$paths$leip_clinical, fixed_lod = file.path(sim_dir, "fixed_lod.csv"),
               output = "output_sim"),
  npx_column = cfg$npx_column, qc = cfg$qc, stats = cfg$stats, enrichment = cfg$enrichment
), file.path(sim_dir, "config_sim.yml"))

msg("Simulated %d samples/controls x %d assays -> %s", nrow(samples), n_assays, sim_dir)
