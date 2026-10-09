# Synthetic stand-in for the O-MicroAD manifest and LEIP clinical file.
# Same sheet layout and study design as the real submission sheet (344 samples: 128 dISF on
# plates 1-2, 216 serum on plates 2-4; MicroAD patients with 2-6 visits, relapse / non-relapse,
# healthy volunteers, CPUO, RELAD, RELAD2, LEIP) but no real IDs, dates or clinical values.
# Lets tests/test_pipeline.R run without any study data.
#   writes data_sim/manifest.xlsx and data_sim/LEIP_clinical.xlsx

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(purrr); library(stringr) })
set.seed(42)
dir.create("data_sim", showWarnings = FALSE)
E <- "EMPTY"

# ---- MicroAD patients: visits and lesion course ---------------------------------------------
ad <- tibble(SubjectID = sprintf("AD%02d", 1:11),
             n_visits = c(5, 3, 6, 5, 6, 4, 6, 6, 6, 2, 6),
             relapse  = c("relapse", "relapse", "relapse", "non-relapse", "non-relapse", "relapse",
                          "non-relapse", "non-relapse", "non-relapse", "DROPOUT", "non-relapse"),
             start = as.Date("2026-01-12") + 14 * (0:10))
ad_visits <- ad |> rowwise() |>
  reframe(SubjectID, relapse, v = seq_len(n_visits), n_visits,
          Date = start + c(0, cumsum(sample(7:21, n_visits - 1, replace = TRUE)))) |>
  mutate(Visit = paste0("V", v),
         L_state = case_when(v == 1 ~ "lesional",
                             relapse == "relapse" & v == n_visits ~ "lesional",
                             relapse == "DROPOUT" ~ "lesional",
                             TRUE ~ "ex-lesional"))
isf_ad <- bind_rows(
  ad_visits |> transmute(SubjectID, Visit, Date, Relapse = relapse, Skin = "L", ClinicalStateSkin = L_state),
  ad_visits |> transmute(SubjectID, Visit, Date, Relapse = relapse, Skin = "NL", ClinicalStateSkin = "non-lesional")
) |> mutate(Group = "AD")
gs   <- tibble(SubjectID = sprintf("GS%02d", 1:12), Visit = "V1", Date = as.Date("2026-01-05") + 10 * (0:11))
cpuo <- tibble(SubjectID = sprintf("CPUO%02d", 1:3), Visit = "V1", Date = as.Date("2026-05-20") + 7 * (0:2))
isf <- bind_rows(
  isf_ad,
  gs |> mutate(Skin = "H", ClinicalStateSkin = "healthy", Group = "GS", Relapse = E),
  cpuo |> crossing(Skin = c("L", "NL")) |>
    mutate(ClinicalStateSkin = if_else(Skin == "L", "lesional", "non-lesional"), Group = "CPOU", Relapse = E)
) |> arrange(SubjectID, Visit, Skin) |>
  mutate(SampleID = sprintf("I%03d", row_number()), SampleType = "dISF", Study = "MicroAD")

serum_micro <- bind_rows(
  ad_visits |> transmute(SubjectID, Visit, Date, Relapse = relapse, Group = "AD"),
  gs |> mutate(Group = "GS", Relapse = E), cpuo |> mutate(Group = "CPOU", Relapse = E)
) |> mutate(Study = "MicroAD")
relad <- tibble(SubjectID = sprintf("RELAD_%02d", 1:35), Study = "RELAD",
                Group = c(rep("H", 7), rep("AD", 28)),
                Relapse = c(rep("healthy", 7), rep("lesional", 4),
                            sample(c("relapse", "non-relapse", "DROPOUT"), 24, TRUE, c(0.45, 0.45, 0.1))))
relad2 <- tibble(SubjectID = sprintf("RELAD2_%02d", 1:76), Study = "RELAD2",
                 Group = c(rep("AD", 68), rep("H", 8)),
                 Relapse = c(sample(c("relapse", "non-relapse", "DROPOUT", "active AD"), 68, TRUE, c(0.45, 0.4, 0.07, 0.08)),
                             rep("healthy", 8)))
leip <- tibble(SubjectID = sprintf("LEIP_%02d", 1:35), Study = "LEIP", Group = E, Relapse = E)
serum <- bind_rows(serum_micro, relad, relad2, leip) |>
  mutate(SampleID = sprintf("S%03d", row_number()), SampleType = "Serum", Skin = E, ClinicalStateSkin = E,
         Visit = coalesce(Visit, E),
         TimeToRelapse = case_when(Study %in% c("RELAD", "RELAD2") & Relapse == "relapse" ~ sample(c("relapse <1w", "relapse >1w"), n(), TRUE),
                                   Study %in% c("RELAD", "RELAD2") & Relapse == "non-relapse" ~ "no relapse",
                                   TRUE ~ E))

# ---- plates: 86 samples + 10 controls each; ISF on plates 1-2, serum on 2-4 --------------------
ctrl <- tibble(well = c(paste0(LETTERS[1:5], 12), paste0(LETTERS[6:8], 12), "G11", "H11"),
               content = c(rep("PC", 5), rep("SC", 3), rep("Neg Ctrl", 2)))
sample_wells <- setdiff(paste0(rep(LETTERS[1:8], each = 12), 1:12), ctrl$well)
# block subjects by plate (whole subjects on one plate), random order within plate
isf_order <- isf |> distinct(SubjectID) |> slice_sample(prop = 1) |> pull(SubjectID)
isf <- isf |> arrange(match(SubjectID, isf_order))
all <- bind_rows(isf |> mutate(TimeToRelapse = E), serum |> slice_sample(prop = 1))
all$plate <- paste("Plate", rep(1:4, each = 86))
all <- all |> group_by(plate) |> slice_sample(prop = 1) |> mutate(well = sample_wells) |> ungroup() |>
  mutate(row = str_sub(well, 1, 1), column = paste("Column", str_sub(well, 2)))

manifest <- all |>
  arrange(SampleType, SampleID) |>
  transmute(SampleID, SubjectID, Visit, Skin, SampleType, Group, SampleNumber = E, Study,
            Date = if_else(is.na(Date), E, format(Date, "%Y-%m-%d")), SampleName = SubjectID,
            Relapse = coalesce(Relapse, E), TimeToRelapse = coalesce(TimeToRelapse, E),
            ClinicalStateSkin = if_else(Study == "RELAD2" & ClinicalStateSkin == E,      # like manifest v4, incl. its spellings
                                        recode(Relapse, relapse = "remisison", `non-relapse` = "remisison", `active AD` = "activeAD",
                                               healthy = "helthy", .default = E), ClinicalStateSkin),
            Sex = if_else(Study == "MicroAD", c("female", "male")[1 + (as.integer(str_extract(SubjectID, "[0-9]+")) %% 3 == 0)], E),
            plate, column, row, well,
            Note = if_else(Relapse %in% "DROPOUT", "DROPOUT", E),
            NoRELAD2 = if_else(Study == "RELAD2" & Relapse != "DROPOUT", as.character(as.integer(str_extract(SubjectID, "[0-9]+$"))), E))

vol <- manifest |> transmute(SampleID, SubjectID, Visit, SampleType, plate, column, row, well,
                             SampleVolume = if_else(SampleType == "dISF", "20", "40"))
vol$SampleVolume[vol$SampleID %in% c("I016", "I019", "I027")] <- c("16", "15", "10.5")
info_sheet <- rbind(matrix(NA, 3, ncol(vol)), names(vol), as.matrix(vol)) |> as.data.frame()
info_sheet[2, 1] <- "Sample list (synthetic)"

layout <- map(1:4, \(p) {
  w <- bind_rows(manifest |> filter(plate == paste("Plate", p)) |> transmute(well, content = SampleID), ctrl)
  grid <- map(LETTERS[1:8], \(r) c(r, w$content[match(paste0(r, 1:12), w$well)])) |> do.call(what = rbind)
  rbind(c(sprintf("PLATE %d", p), rep(NA, 12)), c(NA, 1:12), grid)
}) |> do.call(what = rbind) |> as.data.frame()

manifest_sheet <- rbind(names(manifest), as.matrix(manifest)) |> as.data.frame()
writexl::write_xlsx(list(`Plate Layout` = layout, `Sample Info` = info_sheet, manifest = manifest_sheet),
                    "data_sim/manifest.xlsx", col_names = FALSE)

# ---- LEIP clinical parameters (sheet Key_parameters, same columns as the real file) --------------
lp <- manifest |> filter(Study == "LEIP") |> arrange(SubjectID)
n <- nrow(lp)
clin <- tibble(Olink_SampleID = lp$SampleID, SubjectID = lp$SubjectID, SORB_barcode = 4556000 + seq_len(n),
               age = round(runif(n, 18, 80), 1), sex_MF = sample(c("M", "F"), n, TRUE),
               BMI = round(rnorm(n, 24.5, 3.5), 1), WHR = round(rnorm(n, 0.84, 0.07), 2),
               c_fett = round(rnorm(n, 18, 6), 1), HOMA_IR = round(rlnorm(n, 0.1, 0.6), 2),
               c_CRP = round(rlnorm(n, 0, 0.7), 2), C_CHOL = round(rnorm(n, 5, 0.9), 2),
               C_HDL = round(rnorm(n, 1.6, 0.4), 2), C_LDL = round(rnorm(n, 3, 0.8), 2),
               C_TRIGLY = round(rlnorm(n, 0, 0.4), 2), c_apo = round(rnorm(n, 1.7, 0.3), 2),
               MDRD_kurz = round(rnorm(n, 98, 12), 1), Gluc0_mg_dl = round(rnorm(n, 92, 7), 1))
# LEIP_35 has no clinical data; the simulation still gives it known values (step 08b test checks they are recovered)
clin[n, c("sex_MF", "age", "BMI", "C_HDL")] <- list("M", 72, 31, 1.1)
write.csv(clin[n, c("Olink_SampleID", "SubjectID", "sex_MF", "age", "BMI", "C_HDL")], "data_sim/leip_case_truth.csv", row.names = FALSE)
clin[n, -(1:2)] <- NA                      # like LEIP_35: no clinical data
writexl::write_xlsx(list(Key_parameters = clin), "data_sim/LEIP_clinical.xlsx")
message("Synthetic manifest: ", nrow(manifest), " samples -> data_sim/manifest.xlsx, data_sim/LEIP_clinical.xlsx")
