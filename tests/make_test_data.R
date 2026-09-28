# Invented test data for the LEIP galanin study (no study data needed). Writes tests/data/:
#   npx/LEIP_TEST_Serum_NPX.parquet  Olink Explore HT format: 35 LEIP + 20 other serum samples on 3 plates
#   npx/LEIP_TEST_ISF_NPX.parquet    a file without LEIP samples (must be skipped)
#   fixed_lod.csv                    Olink fixed LOD file format (every 5th assay count-based)
#   LEIP_clinical.xlsx               sheets Key_parameters and All_SORB_parameters, like the real file
#   manifest.xlsx, truth.csv, config_test.yml
# Built-in truth:
#   - the galanin ELISA rises with HDL (also within each sex) and is a little higher in women
#   - Olink GAL follows the ELISA, but less where HDL is high (as if HDL-bound galanin escaped the
#     Olink assay), shares a factor with CHGA, NPY and SCG2, and falls with age (the ELISA does not)
#   - the Olink HDL proteins (APOA1, APOA2, APOM, LCAT, PON1, PON3, CLU) follow HDL, APOB follows LDL
#   - LEP follows BMI and sex; FABP4, LPA, GRN and RARRES2 follow their lab values (A-FABP, Lp(a),
#     progranulin, chemerin); 15 proteins are below LOD; all else is noise
#   - the Olink samples of LEIP_05 and LEIP_06 are swapped (the per-person identity check must find it)
#   - platelet proteins (PF4, PPBP, CXCL5, ... and PEAR1) share a release factor, and Olink GAL follows it
#     in part (the paper figures of scripts/07 must find it); CD40LG is below LOD

source("R/utils.R")
set.seed(20260927)
td <- "tests/data"
unlink(td, recursive = TRUE); dir.create(file.path(td, "npx"), recursive = TRUE)

# ---- LEIP persons and clinical values ---------------------------------------------------------------------------
n <- 35
p <- tibble(SampleID = sprintf("S%03d", 181 + seq_len(n)), SubjectID = sprintf("LEIP_%02d", seq_len(n)),
            plate = rep(c("Plate 2", "Plate 3", "Plate 4"), c(7, 11, 17)), sex = rep(c("M", "F"), length.out = n)) |>
  mutate(female = sex == "F", age = round(runif(n, 18, 80), 1), BMI = round(rnorm(n, 25, 3.5), 1),
         C_HDL = round(1.45 + 0.35 * female - 0.03 * (BMI - 25) + rnorm(n, 0, 0.3), 2),
         c_apo = round(1.0 + 0.45 * C_HDL + rnorm(n, 0, 0.1), 2),
         C_LDL = round(rnorm(n, 3, 0.8), 2), C_TRIGLY = round(rlnorm(n, 0, 0.4), 2),
         C_CHOL = round(C_LDL + C_HDL + C_TRIGLY / 2.2 + rnorm(n, 0, 0.2), 2),
         C_APO_B = round(0.3 + 0.18 * C_LDL + rnorm(n, 0, 0.08), 2), C_LIPO = round(rlnorm(n, log(0.15), 1), 3),
         hdl_z = as.numeric(scale(C_HDL)), ldl_z = as.numeric(scale(C_LDL)), age_z = as.numeric(scale(age)),
         log2_elisa = 7.1 + 0.25 * hdl_z + 0.15 * female + rnorm(n, 0, 0.3),
         galanin = round(2^log2_elisa, 1),
         c_AFABP4 = round(rlnorm(n, log(12), 0.45), 2), IL10 = round(rlnorm(n, log(5), 0.8), 2),
         z_ne = rnorm(n))                                   # shared neuroendocrine factor (GAL, CHGA, NPY, SCG2)
# lab values of two more proteins that Olink also measures; own random stream, so the other test data stay the same
p <- bind_cols(p, local({
  seed <- .Random.seed; on.exit(assign(".Random.seed", seed, envir = globalenv()))
  set.seed(4242)
  tibble(Progranulin = round(rlnorm(n, log(100), 0.35), 1), c_chemerin = round(rlnorm(n, log(180), 0.35), 1),
         z_plt = rnorm(n))                                  # platelet release factor (platelet proteins, GAL)
}))
swap <- c(S186 = "S187", S187 = "S186")                 # Olink samples of LEIP_05 and LEIP_06 swapped

# ---- assays ------------------------------------------------------------------------------------------------------
ht <- readRDS(system.file("extdata", "OlinkID_HT_mapping.rds", package = "OlinkAnalyze")) |>
  distinct(OlinkID, .keep_all = TRUE) |> distinct(Gene, .keep_all = TRUE)
named <- c("GAL", "CHGA", "NPY", "SCG2", "VGF", "CHGB", "APOA1", "APOA2", "APOM", "LCAT", "PON1", "PON3", "CLU",
           "APOB", "LPA", "LEP", "FABP4", "IL10", "GRN", "RARRES2")
assays <- bind_rows(ht |> filter(Gene %in% named), ht |> filter(!Gene %in% named) |> slice_sample(n = 300 - length(named))) |>
  transmute(OlinkID, UniProt, Assay = Gene, Block = as.character(Block), idx = row_number())
na <- nrow(assays)
base <- rnorm(na, 6, 1.5); nc_lvl <- base - 5.5
undetected <- assays$idx %in% sample(which(!assays$Assay %in% named), 15)
role <- case_when(assays$Assay == "GAL" ~ "galanin", assays$Assay %in% c("CHGA", "NPY", "SCG2") ~ "GAL_partner",
                  assays$Assay %in% c("APOA1", "APOA2", "APOM", "LCAT", "PON1", "PON3", "CLU") ~ "HDL_protein",
                  assays$Assay == "APOB" ~ "LDL_protein", assays$Assay %in% c("LEP", "FABP4", "LPA", "GRN", "RARRES2") ~ "clinical_linked",
                  undetected ~ "below_LOD", assays$Assay == "PEAR1" ~ "platelet", TRUE ~ "null")

# ---- samples and controls on the plates ------------------------------------------------------------------------------
others <- tibble(SampleID = sprintf("S%03d", 1:20), plate = rep(c("Plate 3", "Plate 4"), 10))
isf <- tibble(SampleID = sprintf("I%03d", 1:5), plate = "Plate 1")
ctrl <- expand_grid(plate = paste("Plate", 1:4), type = c(rep("PLATE_CONTROL", 5), rep("SAMPLE_CONTROL", 3), rep("NEGATIVE_CONTROL", 2))) |>
  group_by(plate, type) |>
  mutate(SampleID = paste0(recode(type, PLATE_CONTROL = "PC", SAMPLE_CONTROL = "SC", NEGATIVE_CONTROL = "NC"), row_number())) |>
  ungroup()
wells <- bind_rows(p |> transmute(SampleID, plate, SampleType = "SAMPLE", file = "Serum"),
                   others |> mutate(SampleType = "SAMPLE", file = "Serum"),
                   isf |> mutate(SampleType = "SAMPLE", file = "ISF"),
                   ctrl |> transmute(SampleID, plate, SampleType = type, file = if_else(plate == "Plate 1", "ISF", "Serum"))) |>
  group_by(plate) |> mutate(WellID = paste0(LETTERS[(row_number() - 1) %/% 12 + 1], (row_number() - 1) %% 12 + 1)) |> ungroup()
plate_off <- matrix(rnorm(4 * na, 0, 0.08), 4, dimnames = list(paste("Plate", 1:4), NULL))

npx <- matrix(NA_real_, nrow(wells), na)
for (i in seq_len(nrow(wells))) {
  w <- wells[i, ]
  if (w$SampleType == "NEGATIVE_CONTROL") { npx[i, ] <- nc_lvl + rnorm(na, 0, 0.3); next }
  if (w$SampleType == "PLATE_CONTROL")    { npx[i, ] <- base + rnorm(na, 0, 0.08); next }
  if (w$SampleType == "SAMPLE_CONTROL")   { npx[i, ] <- base + 0.2 + rnorm(na, 0, 0.12); next }
  e <- rnorm(na, 0, 0.35); e[role == "galanin"] <- 0.5 * e[role == "galanin"]   # GAL: far above LOD, measured precisely
  x <- base + plate_off[w$plate, ] + e - if (w$file == "ISF") 1 else 0
  k <- match(coalesce(swap[w$SampleID], w$SampleID), p$SampleID)
  if (!is.na(k)) {
    q <- p[k, ]
    x <- x + (role == "galanin") * (1.5 * (q$log2_elisa - 7.1) - 0.2 * q$hdl_z + 0.5 * q$z_ne - 0.4 * q$age_z + 0.35 * q$z_plt) +
      (role == "platelet") * 0.9 * q$z_plt +
      (role == "GAL_partner") * 0.9 * q$z_ne + (role == "HDL_protein") * 0.5 * q$hdl_z + (role == "LDL_protein") * 0.5 * q$ldl_z +
      (assays$Assay == "LEP") * (0.15 * (q$BMI - 25) + 0.5 * q$female) + (assays$Assay == "FABP4") * log2(q$c_AFABP4 / 12) +
      (assays$Assay == "LPA") * 0.5 * log2(q$C_LIPO / 0.15) + (assays$Assay == "GRN") * log2(q$Progranulin / 100) +
      (assays$Assay == "RARRES2") * log2(q$c_chemerin / 180)
  }
  x[undetected] <- nc_lvl[undetected] + 0.3 + rnorm(sum(undetected), 0, 0.3)
  npx[i, ] <- x
}

# ---- more platelet proteins: own random stream, so all other test data stay the same -----------------------------------
extra <- ht |> filter(Gene %in% c("PF4", "PPBP", "CXCL5", "CCL5", "TGFB1", "BDNF", "EGF", "SELP", "GP5", "GP6", "CD40LG"),
                      !Gene %in% assays$Assay) |>
  transmute(OlinkID, UniProt, Assay = Gene, idx = na + row_number())
extra$Block <- sort(unique(assays$Block))[(seq_len(nrow(extra)) - 1) %% n_distinct(assays$Block) + 1]   # existing blocks only
ex <- local({
  seed <- .Random.seed; on.exit(assign(".Random.seed", seed, envir = globalenv()))
  set.seed(777)
  k <- nrow(extra); b <- rnorm(k, 6, 1.5); dead <- extra$Assay == "CD40LG"
  off <- matrix(rnorm(4 * k, 0, 0.08), 4, dimnames = list(paste("Plate", 1:4), NULL))
  m <- matrix(NA_real_, nrow(wells), k)
  for (i in seq_len(nrow(wells))) {
    w <- wells[i, ]
    if (w$SampleType == "NEGATIVE_CONTROL") { m[i, ] <- b - 5.5 + rnorm(k, 0, 0.3); next }
    if (w$SampleType == "PLATE_CONTROL")    { m[i, ] <- b + rnorm(k, 0, 0.08); next }
    if (w$SampleType == "SAMPLE_CONTROL")   { m[i, ] <- b + 0.2 + rnorm(k, 0, 0.12); next }
    x <- b + off[w$plate, ] + rnorm(k, 0, 0.35) - if (w$file == "ISF") 1 else 0
    j <- match(coalesce(swap[w$SampleID], w$SampleID), p$SampleID)
    if (!is.na(j)) x <- x + 0.9 * p$z_plt[j]
    x[dead] <- b[dead] - 5.5 + 0.3 + rnorm(sum(dead), 0, 0.3)
    m[i, ] <- x
  }
  list(npx = m, base = b, role = if_else(dead, "below_LOD", "platelet"))
})
assays <- bind_rows(assays, extra |> select(names(assays)))
npx <- cbind(npx, ex$npx); base <- c(base, ex$base); nc_lvl <- base - 5.5; role <- c(role, ex$role); na <- nrow(assays)

# ---- Olink parquet files (Explore HT layout) -------------------------------------------------------------------------
blocks <- sort(unique(assays$Block))
ext <- expand_grid(i = seq_len(nrow(wells)), Block = blocks) |> mutate(ExtCount = as.integer(round(5000 * 2^rnorm(n(), 0, 0.1))))
long <- expand_grid(i = seq_len(nrow(wells)), a = seq_len(na)) |>
  mutate(SampleID = wells$SampleID[i], SampleType = wells$SampleType[i], WellID = wells$WellID[i],
         PlateID = str_remove(wells$plate[i], " "), file = wells$file[i], DataAnalysisRefID = "TEST_DAR_1",
         OlinkID = assays$OlinkID[a], UniProt = assays$UniProt[a], Assay = assays$Assay[a], AssayType = "assay",
         Panel = "Explore_HT", Block = assays$Block[a], ExtNPX = npx[cbind(i, a)] - 4) |>
  left_join(ext, by = c("i", "Block")) |>
  group_by(PlateID, OlinkID) |> mutate(PCNormalizedNPX = ExtNPX - median(ExtNPX[SampleType == "PLATE_CONTROL"])) |> ungroup() |>
  mutate(NPX = PCNormalizedNPX, Count = as.integer(round(ExtCount * 2^ExtNPX)), Normalization = "Plate control",
         AssayQC = "PASS", SampleQC = "PASS", ExploreVersion = "TEST") |> select(-a)
ctrl_assays <- expand_grid(i = seq_len(nrow(wells)), Block = blocks, ct = c("ext_ctrl", "inc_ctrl", "amp_ctrl")) |>
  left_join(ext, by = c("i", "Block")) |>
  mutate(SampleID = wells$SampleID[i], SampleType = wells$SampleType[i], WellID = wells$WellID[i],
         PlateID = str_remove(wells$plate[i], " "), file = wells$file[i], DataAnalysisRefID = "TEST_DAR_1",
         OlinkID = sprintf("OID9%d%03d", match(ct, c("ext_ctrl", "inc_ctrl", "amp_ctrl")), as.integer(Block)),
         UniProt = NA_character_, Assay = paste0(ct, "_", Block), AssayType = ct, Panel = "Explore_HT",
         ExtNPX = if_else(ct == "ext_ctrl", 0, rnorm(n(), 0, 0.2)), NPX = ExtNPX, PCNormalizedNPX = ExtNPX,
         Count = if_else(ct == "ext_ctrl", ExtCount, 5000L), Normalization = "Plate control",
         AssayQC = "PASS", SampleQC = "PASS", ExploreVersion = "TEST") |> select(-ct)
long <- bind_rows(long, ctrl_assays) |> select(-i, -ExtCount)
long$SampleQC[long$SampleID == "S216"] <- "FAIL"   # LEIP_35 fails Olink QC (excluded), LEIP_09 gets a warning (kept)
long$SampleQC[long$SampleID == "S190"] <- "WARN"
for (mx in c("Serum", "ISF")) {
  tbl <- arrow::arrow_table(long |> filter(file == mx) |> select(-file))
  tbl$metadata <- list(FileVersion = "NA", ProjectName = "LEIP-TEST", SampleMatrix = mx, DataFileType = "NPX File", Product = "ExploreHT")
  arrow::write_parquet(tbl, file.path(td, "npx", sprintf("LEIP_TEST_%s_NPX.parquet", mx)))
}
lod <- assays |>
  transmute(OlinkID, AssayType = "assay", UniProt, Assay, Panel = "Explore_HT", Block, DataAnalysisRefID = "TEST_DAR_1",
            BimodalDistribution = FALSE, LODNPX = nc_lvl + 0.9 - base, LODCount = as.integer(pmax(150, round(5000 * 2^(nc_lvl + 0.9 - 4)))),
            LODMethod = if_else(idx %% 5 == 0, "lod_count", "lod_npx"), Version = "10.2.0")
write.table(lod, file.path(td, "fixed_lod.csv"), sep = ";", row.names = FALSE, quote = FALSE)

# ---- clinical file, like the real one (non-ASCII column names included) --------------------------------------------------
mu <- intToUtf8(0xb5); ae <- intToUtf8(0xe4)
key <- p |>
  transmute(Olink_SampleID = SampleID, SubjectID, SORB_barcode = 4556000 + row_number(), age, sex_MF = sex, BMI,
            WHR = round(rnorm(n, 0.85, 0.07), 2), c_fett = round(rnorm(n, 20, 6), 1), NGT_IGT_IFG_DIAB = 0, t2d = 0,
            Gluc0_mg_dl = round(rnorm(n, 92, 7), 1), ins0 = round(rlnorm(n, log(4), 0.5), 2), HOMA_IR = round(ins0 * Gluc0_mg_dl / 405, 2),
            c_CRP = round(pmax(0.5, rlnorm(n, 0, 0.8)), 2), C_CHOL, C_HDL, C_LDL, C_TRIGLY, c_apo,
            MDRD_kurz = round(rnorm(n, 98, 12), 1), SMOKING_current = 0, low_serum = c(1, rep(0, n - 1)), lip = 0,
            `Galanin [pg/mL]` = galanin)
names(key)[names(key) == "ins0"] <- paste0("Ins0_", mu, "U_ml")
names(key)[names(key) == "lip"] <- paste0("lip", ae, "misch")
sorb <- key |>
  transmute(Olink_SampleID, SubjectID, Olink_plate = p$plate, Olink_well = wells$WellID[match(p$SampleID, wells$SampleID)],
            SORB_barcode, `Galanin [pg/mL]`, Galanin_ELISA_plate = sample(1:12, n, TRUE), sex = if_else(sex_MF == "M", 1, 0), sex_MF,
            age, BMI, ln_BMI = round(log(BMI), 3), Gluc0_mg_dl, gluk_0 = round(Gluc0_mg_dl / 18.016, 3),
            C_APO_B = p$C_APO_B, C_LIPO = p$C_LIPO, c_AFABP4 = p$c_AFABP4, IL10 = p$IL10, RESTRAINT = sample(0:15, n, TRUE),
            RE_BIN = as.integer(RESTRAINT > 7), c_tsh = round(rlnorm(n, log(1.8), 0.4), 2), Progranulin = p$Progranulin,
            c_chemerin = p$c_chemerin)
key[n, -(1:2)] <- NA; sorb[n, -c(1:4)] <- NA          # like LEIP_35: no clinical data
writexl::write_xlsx(list(README = tibble(Item = "Purpose", Note = "Invented test data"), Key_parameters = key,
                         All_SORB_parameters = sorb), file.path(td, "LEIP_clinical.xlsx"))
man <- wells |> filter(SampleType == "SAMPLE") |>
  transmute(SampleID, SubjectID = coalesce(p$SubjectID[match(SampleID, p$SampleID)], paste0("X", SampleID)),
            Study = if_else(SampleID %in% p$SampleID, "LEIP", "OTHER"), plate)
writexl::write_xlsx(list(manifest = man), file.path(td, "manifest.xlsx"))
write_csv(assays |> mutate(role) |> select(OlinkID, Assay, role), file.path(td, "truth.csv"))

cfg <- yaml::read_yaml("config.yml")
cfg$paths <- list(npx_dir = file.path(td, "npx"), fixed_lod = file.path(td, "fixed_lod.csv"),
                  leip_clinical = file.path(td, "LEIP_clinical.xlsx"), manifest = file.path(td, "manifest.xlsx"),
                  output = "tests/output")
cfg$lab_vs_olink <- list(galanin_elisa = "GAL", c_apo = "APOA1", c_AFABP4 = "FABP4", C_LIPO = "LPA", Progranulin = "GRN",
                         c_chemerin = "RARRES2", IL10 = "IL10", c_CRP = "CRP")
cfg$expected_associations <- list(list(protein = "LEP", parameter = "BMI", direction = "positive"),
                                  list(protein = "APOA1", parameter = "C_HDL", direction = "positive"),
                                  list(protein = "NOT_ON_PANEL", parameter = "BMI", direction = "positive"))
cfg$bootstrap <- 500; cfg$permutations <- 200; cfg$platelets$random_sets <- 2000
yaml::write_yaml(cfg, file.path(td, "config_test.yml"))
message("Test data: ", nrow(p), " LEIP samples, ", na, " assays -> ", td)
