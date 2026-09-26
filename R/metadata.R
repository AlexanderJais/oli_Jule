# Sample metadata: read the Olink submission sheet, harmonise codes, derive analysis
# variables and flag inconsistencies. Nothing is dropped here; flags are carried along.

excel_date <- function(x) {
  num <- suppressWarnings(as.numeric(x))
  out <- as.Date(num, origin = "1899-12-30")
  txt <- is.na(num) & !is.na(x)
  # each text date is tried against each format separately (mixed formats in one column)
  out[txt] <- dplyr::coalesce(as.Date(x[txt], format = "%Y-%m-%d"), as.Date(x[txt], format = "%d.%m.%Y"),
                              as.Date(x[txt], format = "%d/%m/%Y"))
  out
}

read_manifest_sheet <- function(path) {
  readxl::read_excel(path, sheet = "manifest", col_types = "text", na = c("", "EMPTY", "NA")) |>
    mutate(across(everything(), str_trim))
}

#' Sample volumes from the "Sample Info" sheet (header row is found automatically).
read_sample_volumes <- function(path) {
  raw <- suppressMessages(readxl::read_excel(path, sheet = "Sample Info", col_names = FALSE, col_types = "text"))
  hdr <- which(raw[[1]] == "SampleID")[1]
  if (is.na(hdr)) return(tibble(SampleID = character(), volume_ul = numeric()))
  h <- as.character(unlist(raw[hdr, ]))
  tibble(SampleID = raw[[1]][-seq_len(hdr)],
         volume_ul = suppressWarnings(as.numeric(raw[[which(h == "SampleVolume")[1]]][-seq_len(hdr)]))) |>
    filter(!is.na(SampleID))
}

#' Wells of all plates from the "Plate Layout" sheet: plate, well, content.
read_plate_layout <- function(path) {
  raw <- suppressMessages(readxl::read_excel(path, sheet = "Plate Layout", col_names = FALSE, col_types = "text"))
  out <- list(); plate <- NA_character_
  for (i in seq_len(nrow(raw))) {
    first <- raw[[1]][i]
    if (!is.na(first) && str_detect(first, "^PLATE\\s*\\d+")) {
      plate <- paste("Plate", str_extract(first, "\\d+")); next
    }
    if (!is.na(plate) && !is.na(first) && first %in% LETTERS[1:8]) {
      out[[length(out) + 1]] <- tibble(plate = plate, well = paste0(first, 1:12),
                                       content = as.character(unlist(raw[i, 2:13])))
    }
  }
  bind_rows(out)
}

harmonise_relapse <- function(x) {
  x <- str_to_lower(x)
  case_when(
    x %in% c("relapse") ~ "relapse",
    x %in% c("non-relapse", "no relapse", "nonrelapse") ~ "non-relapse",
    x == "dropout" ~ "dropout",
    x == "healthy" ~ "healthy",
    x %in% c("lesional", "active ad") ~ "active",
    TRUE ~ NA_character_
  )
}

harmonise_time_to_relapse <- function(x) {
  x <- str_to_lower(x)
  case_when(
    str_detect(x, "<\\s*1\\s*w") ~ "<1w",
    str_detect(x, ">\\s*1\\s*w") ~ ">1w",
    str_detect(x, "no relapse") ~ "no relapse",
    TRUE ~ x
  )
}

#' Build the clean sample table from the manifest workbook.
build_metadata <- function(manifest_path, leip_path = NULL) {
  m <- read_manifest_sheet(manifest_path)
  req <- c("SampleID", "SubjectID", "Visit", "Skin", "SampleType", "Group", "Study",
           "Date", "SampleName", "Relapse", "TimeToRelapse", "ClinicalStateSkin",
           "plate", "well", "Note")
  miss <- setdiff(req, names(m))
  if (length(miss)) stop("Manifest is missing columns: ", paste(miss, collapse = ", "))

  meta <- m |>
    transmute(
      SampleID, SubjectID, SampleName,
      matrix = recode(SampleType, dISF = "ISF", Serum = "Serum"),
      cohort = Study,
      group_raw = Group,
      group = case_when(
        Group %in% c("CPOU", "CPUO") ~ "CPUO",
        Group %in% c("GS", "H") ~ "HC",
        Study == "LEIP" ~ "Biobank",
        TRUE ~ Group
      ),
      visit = Visit,
      visit_num = suppressWarnings(as.integer(str_remove(Visit, "^V"))),
      date = excel_date(Date),
      site = if_else(SampleType == "dISF", Skin, NA_character_),
      state = if_else(SampleType == "dISF", ClinicalStateSkin, NA_character_),
      relapse_raw = Relapse,
      relapse = harmonise_relapse(Relapse),
      time_to_relapse = harmonise_time_to_relapse(TimeToRelapse),
      dropout = coalesce(Note == "DROPOUT", FALSE) | coalesce(relapse == "dropout", FALSE),
      plate, well
    ) |>
    mutate(across(c(matrix, cohort, group), as.character))

  # --- timing (MicroAD) ---------------------------------------------------------
  meta <- meta |>
    group_by(SubjectID) |>
    mutate(days_since_v1 = if (all(is.na(date))) NA_real_ else
      as.numeric(date - min(date[visit_num == min(visit_num, na.rm = TRUE)], na.rm = TRUE))) |>
    ungroup()

  # --- relapse visit: first return of 'lesional' at the L site after it had cleared
  l_site <- meta |>
    filter(matrix == "ISF", site == "L", group == "AD", !is.na(visit_num)) |>
    arrange(SubjectID, visit_num)
  relapse_visits <- l_site |>
    group_by(SubjectID) |>
    summarise(
      first_cleared = suppressWarnings(min(visit_num[state %in% "ex-lesional"], na.rm = TRUE)),
      relapse_visit = suppressWarnings(min(visit_num[state %in% "lesional" & visit_num > first_cleared], na.rm = TRUE)),
      .groups = "drop") |>
    mutate(relapse_visit = as.integer(if_else(is.finite(relapse_visit), relapse_visit, NA_real_))) |>
    select(SubjectID, relapse_visit)
  lesion_at_visit <- l_site |>
    transmute(SubjectID, visit_num,
              lesion_state = recode(state, lesional = "active", `ex-lesional` = "cleared",
                                    `non-lesional` = "none"))

  meta <- meta |>
    left_join(relapse_visits, by = "SubjectID") |>
    left_join(lesion_at_visit, by = c("SubjectID", "visit_num")) |>
    mutate(
      visits_to_relapse = relapse_visit - visit_num,
      pre_relapse = !is.na(relapse_visit) & visit_num < relapse_visit
    )

  # --- sample volume -------------------------------------------------------------
  meta <- meta |> left_join(read_sample_volumes(manifest_path), by = "SampleID")

  # --- LEIP clinical data (age, sex, BMI ...) ------------------------------------
  if (!is.null(leip_path) && file.exists(leip_path)) {
    clin <- readxl::read_excel(leip_path, sheet = "Key_parameters") |>
      rename(SampleID = Olink_SampleID) |>
      select(-any_of(c("SubjectID"))) |>
      rename(sex = sex_MF)
    meta <- meta |> left_join(clin, by = "SampleID")
  }

  meta |> arrange(matrix, SampleID)
}

#' Rule-based consistency checks. Returns one row per (sample, issue).
flag_metadata <- function(meta, low_volume = list(ISF = 20, Serum = 40)) {
  f <- list()
  add <- function(ids, issue) if (length(ids)) f[[length(f) + 1]] <<- tibble(SampleID = ids, issue = issue)

  micro <- meta |> filter(cohort == "MicroAD")
  dm <- micro |> group_by(SubjectID, visit) |> filter(n_distinct(date, na.rm = TRUE) > 1) |> ungroup()
  add(dm$SampleID, "different dates within one subject-visit")

  nm <- micro |>
    group_by(SubjectID, visit_num) |> summarise(d = min(date), .groups = "drop") |>
    arrange(SubjectID, visit_num) |> group_by(SubjectID) |>
    filter(any(diff(as.numeric(d)) < 0, na.rm = TRUE)) |> ungroup()
  add(micro$SampleID[micro$SubjectID %in% nm$SubjectID], "visit dates not in chronological order")

  add(meta$SampleID[meta$matrix == "ISF" & meta$site == "L" & meta$visit_num == 1 &
                      meta$state == "ex-lesional" & meta$group == "AD"],
      "L site ex-lesional at baseline (V1)")

  add(meta$SampleID[meta$group == "AD" & coalesce(meta$relapse == "healthy", FALSE)],
      "Group AD but Relapse = healthy")
  add(meta$SampleID[meta$group == "HC" & str_detect(coalesce(meta$SampleName, ""), "^AD")],
      "Group H but SampleName starts with AD")
  add(meta$SampleID[coalesce(meta$relapse == "relapse" & meta$time_to_relapse == "no relapse", FALSE)],
      "Relapse = relapse but TimeToRelapse = no relapse")
  add(meta$SampleID[coalesce(meta$relapse == "non-relapse" & str_detect(meta$time_to_relapse, "^[<>]"), FALSE)],
      "Relapse = non-relapse but TimeToRelapse gives a relapse time")
  add(meta$SampleID[meta$cohort %in% c("RELAD", "RELAD2") & coalesce(meta$relapse == "relapse", FALSE) &
                      is.na(meta$time_to_relapse)],
      "relapse without TimeToRelapse")

  lv <- unlist(low_volume)
  add(meta$SampleID[!is.na(meta$volume_ul) & meta$volume_ul < lv[meta$matrix]], "low sample volume")
  if ("age" %in% names(meta)) add(meta$SampleID[meta$cohort == "LEIP" & is.na(meta$age)], "LEIP sample without clinical data")

  if (!length(f)) return(tibble(SampleID = character(), issue = character()))
  bind_rows(f) |> distinct()
}

#' Hard checks: stop if IDs or plate positions are inconsistent.
validate_metadata <- function(meta, layout) {
  stopifnot(
    "duplicated SampleID in manifest" = !anyDuplicated(meta$SampleID),
    "duplicated plate/well in manifest" = !anyDuplicated(meta[c("plate", "well")]),
    "unknown matrix (SampleType must be dISF or Serum)" = all(meta$matrix %in% c("ISF", "Serum"))
  )
  if (nrow(layout)) {
    chk <- meta |> select(SampleID, plate, well) |>
      left_join(layout, by = c("plate", "well")) |>
      filter(is.na(content) | content != SampleID)
    if (nrow(chk)) {
      print(chk)
      stop(nrow(chk), " sample(s) sit in a different well in 'Plate Layout' than in 'manifest'.")
    }
  }
  dup_isf <- meta |> filter(matrix == "ISF") |> count(SubjectID, visit, site) |> filter(n > 1)
  if (nrow(dup_isf)) {
    print(dup_isf)
    stop("More than one ISF sample for the same subject, visit and site.")
  }
  invisible(TRUE)
}
