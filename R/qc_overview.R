# QC overview (step 02b): proteins above / below LOD per matrix, overall (ring charts) and by protein class
# (Human Protein Atlas "Protein class"). Used by scripts/02b_qc_overview.R and the executive summary (step 18).

library(grid)

# The 14 protein classes shown, and the Human Protein Atlas "Protein class" values each one combines.
qc_protein_classes <- list(
  "Enzymes (incl. metabolic)"                        = c("Enzymes", "Metabolic proteins"),
  "Transcription factors"                            = "Transcription factors",
  "Nuclear receptors"                                = "Nuclear receptors",
  "GPCRs"                                            = "G-protein coupled receptors",
  "Ion channels (voltage-gated)"                     = "Voltage-gated ion channels",
  "Transporters"                                     = "Transporters",
  "Drug related (FDA-approved + potential targets)"  = c("FDA approved drug targets", "Potential drug targets"),
  "Disease related"                                  = c("Disease related genes", "Human disease related genes"),
  "Cancer related"                                   = "Cancer-related genes",
  "Immune related (CD markers, Ig, TCR)"             = c("CD markers", "Immunoglobulin genes", "T-cell receptor genes"),
  "Essential (DepMap)"                               = "Essential proteins",
  "Intracellular"                                    = "Predicted intracellular proteins",
  "Membrane"                                         = "Predicted membrane proteins",
  "Extracellular / secreted (incl. plasma proteins)" = c("Predicted secreted proteins", "Plasma proteins"))
qc_class_block <- setNames(rep(c("function", "relevance", "localisation"), c(6, 5, 3)), names(qc_protein_classes))

#' Protein classes per Olink assay (long: OlinkID, class). Matching per component of combined assays
#' (e.g. "IL12A_IL12B"): UniProt first, then gene symbol, then gene synonym; an assay gets the union of the
#' classes of its components. Returns list(classes = long table, matched = OlinkIDs found in the HPA).
assay_protein_classes <- function(assays, hpa) {
  h <- hpa |> mutate(.row = row_number())
  col <- \(nm) if (nm %in% names(h)) h[[nm]] else rep(NA_character_, nrow(h))
  by_up  <- tibble(.row = h$.row, key = col("Uniprot")) |> separate_rows(key, sep = ",\\s*") |> filter(!is.na(key), key != "")
  by_sym <- tibble(.row = h$.row, key = col("Gene")) |> filter(!is.na(key))
  by_syn <- tibble(.row = h$.row, key = col("Gene synonym")) |> separate_rows(key, sep = ",\\s*") |> filter(!is.na(key), key != "")
  comp <- assays |> distinct(OlinkID, Assay, UniProt) |>
    mutate(up = str_split(coalesce(UniProt, ""), "_"), g = str_split(coalesce(Assay, ""), "_")) |>
    mutate(k = map2(up, g, \(a, b) seq_len(max(length(a), length(b))))) |>
    select(OlinkID, up, g, k) |>
    unnest(c(k)) |> mutate(up = map2_chr(up, k, \(a, i) if (i <= length(a)) a[i] else NA_character_),
                           g = map2_chr(g, k, \(a, i) if (i <= length(a)) a[i] else NA_character_))
  m1 <- comp |> inner_join(by_up, by = c(up = "key"), relationship = "many-to-many")
  rest <- comp |> anti_join(m1, by = c("OlinkID", "k"))
  m2 <- rest |> inner_join(by_sym, by = c(g = "key"), relationship = "many-to-many")
  rest <- rest |> anti_join(m2, by = c("OlinkID", "k"))
  m3 <- rest |> inner_join(by_syn, by = c(g = "key"), relationship = "many-to-many")
  hit <- bind_rows(m1, m2, m3) |> distinct(OlinkID, .row)
  pc <- tibble(.row = h$.row, pc = col("Protein class")) |> separate_rows(pc, sep = ",\\s*") |> filter(!is.na(pc), pc != "")
  apc <- hit |> inner_join(pc, by = ".row", relationship = "many-to-many") |> distinct(OlinkID, pc)
  cls <- imap(qc_protein_classes, \(v, nm) apc |> filter(pc %in% v) |> distinct(OlinkID) |> mutate(class = nm)) |> bind_rows()
  list(classes = cls, matched = unique(hit$OlinkID))
}

#' Per assay and matrix: fraction of samples above LOD and the above / below LOD call.
#' rule "all_samples": above LOD in >= min_frac of all samples of the matrix;
#' rule "analysis_filter": the detection filter of step 02 (>= min_frac in at least one group; column keep).
lod_status <- function(clean, det, rule = "all_samples", min_frac = 0.5, serum_cohorts = "MicroAD") {
  d <- clean |> filter(matrix == "ISF" | (matrix == "Serum" & (length(serum_cohorts) == 0 | cohort %in% serum_cohorts)))
  st <- d |> group_by(matrix, OlinkID, Assay, UniProt) |>
    summarise(n_samples = n_distinct(SampleID), n_with_lod = sum(!is.na(below_lod)),
              frac_above_LOD = if (all(is.na(below_lod))) NA_real_ else mean(!below_lod, na.rm = TRUE), .groups = "drop") |>
    left_join(det |> select(matrix, OlinkID, analysis_filter = keep), by = c("matrix", "OlinkID")) |>
    mutate(status = case_when(is.na(frac_above_LOD) ~ "no LOD",
                              rule == "analysis_filter" & analysis_filter %in% TRUE ~ "above LOD",
                              rule == "analysis_filter" ~ "below LOD",
                              frac_above_LOD >= min_frac ~ "above LOD",
                              TRUE ~ "below LOD"),
           matrix = recode(matrix, ISF = "dISF"))
  st
}

#' Font handling. Nimbus Sans ships with the project (fonts/, URW base35, AGPL-3 with font exception) and is drawn by
#' the package showtext as vector outlines: identical on every computer, but the PDF lists no font and the figure text
#' cannot be selected. Other fonts: the installed font via cairo, else standard Helvetica (a PDF base font).
font_installed <- function(family) {
  if (requireNamespace("systemfonts", quietly = TRUE))
    return(any(tolower(systemfonts::system_fonts()$family) == tolower(family)))
  if (nzchar(Sys.which("fc-list"))) {
    fams <- tryCatch(system2("fc-list", c(":", "family"), stdout = TRUE, stderr = FALSE), error = \(e) character())
    return(any(tolower(trimws(unlist(strsplit(fams, ",")))) == tolower(family)))
  }
  lad <- Sys.getenv("LOCALAPPDATA")
  dirs <- c("~/Library/Fonts", "/Library/Fonts", "/System/Library/Fonts", "C:/Windows/Fonts",
            if (nzchar(lad)) file.path(lad, "Microsoft", "Windows", "Fonts"))
  files <- unlist(lapply(dirs[dir.exists(dirs)], list.files, recursive = TRUE))
  any(grepl(gsub("\\s+", "", family), gsub("[\\s_-]+", "", files, perl = TRUE), ignore.case = TRUE))
}

#' How the figure font is provided: list(mode = "showtext" | "cairo" | "builtin", family, message).
qc_font_setup <- function(font = "Nimbus Sans", font_dir = "fonts") {
  files <- file.path(font_dir, c("NimbusSans-Regular.otf", "NimbusSans-Bold.otf"))
  if (grepl("^nimbus ?sans", font, ignore.case = TRUE) && all(file.exists(files)) &&
      requireNamespace("showtext", quietly = TRUE) && requireNamespace("sysfonts", quietly = TRUE)) {
    if (!"NimbusSansQC" %in% sysfonts::font_families())
      sysfonts::font_add("NimbusSansQC", regular = files[1], bold = files[2])
    return(list(mode = "showtext", family = "NimbusSansQC", message = NULL))
  }
  if (font_installed(font) && capabilities("cairo")) return(list(mode = "cairo", family = font, message = NULL))
  list(mode = "builtin", family = "",
       message = sprintf("Font '%s' is not in fonts/ and not installed: the figure uses Helvetica instead.", font))
}

#' Open a PDF / PNG device for the figure with the font set up; close it with close_device().
open_device <- function(file, type = c("pdf", "png"), fi = qc_font_setup(), width = 11.69, height = 8.27) {
  type <- match.arg(type)
  if (type == "pdf") {
    if (fi$mode == "cairo") grDevices::cairo_pdf(file, width = width, height = height, family = fi$family)
    else grDevices::pdf(file, width = width, height = height, family = "Helvetica")
  } else {
    grDevices::png(file, width = width, height = height, units = "in", res = 300,
                   type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"),
                   family = if (fi$mode == "cairo") fi$family else "sans")
  }
  if (fi$mode == "showtext") { showtext::showtext_opts(dpi = if (type == "png") 300 else 96); showtext::showtext_begin() }
  invisible(fi)
}
close_device <- function(fi) {
  if (fi$mode == "showtext") showtext::showtext_end()
  invisible(grDevices::dev.off())
}

#' Use font family `fam` for all text of a ggplot (theme and text layers).
set_family <- function(p, fam) {
  if (is.null(p) || !nzchar(fam)) return(p)
  p <- p + theme(text = element_text(family = fam))
  for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomText")) p$layers[[i]]$aes_params$family <- fam
  p
}

#' Ring charts (one per matrix) and class bars; returns list(rings, classes) of ggplot objects.
qc_overview_plots <- function(st, cls, colours) {
  fills <- c(`above LOD` = colours$above %||% "#69005F", `below LOD` = colours$below %||% "#FF506E", `no LOD` = colours$no_lod %||% "#BDBDBD")
  mx <- c("Serum", "dISF")
  ring <- st |> count(matrix, status) |>
    mutate(matrix = factor(matrix, mx), status = factor(status, c("above LOD", "below LOD", "no LOD"))) |>
    arrange(matrix, status) |>
    # shares, not counts: each ring is a full circle even if serum and dISF have different totals
    group_by(matrix) |> mutate(total = sum(n), pct = n / total, ymax = cumsum(n) / total, ymin = ymax - n / total) |>
    # label just outside the ring, aligned away from it (justification follows the angle smoothly, so the gap
    # to the ring stays even); two labels in the top band (e.g. a tiny "no LOD" slice) are pushed apart sideways
    mutate(ang = 2 * pi * (ymin + ymax) / 2,
           hjust = pmin(1, pmax(0, 0.5 - 2 * sin(ang))),
           vjust = pmin(1, pmax(0, 0.5 - cos(ang))),
           top = cos(ang) > 0.5,
           pushed = sum(top) > 1 & top,
           hjust = if (sum(top) > 1) if_else(top, if_else(sin(ang) == min(sin(ang)[top]), 1, 0), hjust) else hjust) |>
    ungroup() |>
    mutate(label = sprintf("%s\n%.0f%%", format(n, big.mark = ",", trim = TRUE), 100 * pct),
           # a small gap between two labels pushed apart at the top
           label = case_when(pushed & hjust == 1 ~ gsub("(\n|$)", "   \\1", label),
                             pushed & hjust == 0 ~ gsub("(^|\n)", "\\1   ", label),
                             TRUE ~ label))
  rings <- ggplot(ring) +
    geom_rect(aes(xmin = 2, xmax = 3.4, ymin = ymin, ymax = ymax, fill = status), colour = "white", linewidth = 0.8) +
    geom_text(aes(x = 3.75, y = (ymin + ymax) / 2, hjust = hjust, vjust = vjust, label = label),
              size = 4, lineheight = 0.9, colour = "grey15") +
    geom_text(data = distinct(ring, matrix, total), aes(x = 0, y = 0, label = paste0(matrix, "\n", format(total, big.mark = ",", trim = TRUE))),
              size = 5.2, fontface = "bold", lineheight = 0.95, colour = "grey10") +
    coord_polar(theta = "y", clip = "off") + xlim(0, 4.5) + facet_wrap(~matrix) +
    scale_fill_manual(values = fills, name = NULL, drop = TRUE) + theme_void(base_size = 12) +
    theme(strip.text = element_blank(), legend.position = "top", legend.text = element_text(size = 12),
          panel.spacing = unit(4.5, "cm"))   # rings sit above the serum and dISF bar panels
  if (is.null(cls) || !nrow(cls)) return(list(rings = rings, classes = NULL))
  lv <- names(qc_protein_classes)
  cb <- st |> filter(status != "no LOD") |> inner_join(cls, by = "OlinkID", relationship = "many-to-many") |>
    count(matrix, class, status) |>
    complete(matrix = mx, class = lv, status = c("above LOD", "below LOD"), fill = list(n = 0)) |>
    group_by(matrix, class) |> mutate(total = sum(n)) |> ungroup() |>
    mutate(matrix = factor(matrix, mx), block = factor(qc_class_block[class], c("function", "relevance", "localisation")),
           class = factor(class, rev(lv)), status = factor(status, c("below LOD", "above LOD")))
  pct <- cb |> filter(status == "above LOD") |> mutate(lab = if_else(total > 0, sprintf("%.0f%%", 100 * n / total), ""))
  classes <- ggplot(cb, aes(n, class, fill = status)) +
    geom_col(width = 0.72) +
    geom_text(data = pct, aes(x = total, y = class, label = lab), inherit.aes = FALSE, hjust = -0.15, size = 3.2, colour = "grey20") +
    facet_grid(block ~ matrix, scales = "free_y", space = "free_y") +
    scale_fill_manual(values = fills, breaks = c("above LOD", "below LOD"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12)), labels = \(x) format(x, big.mark = ",", trim = TRUE)) +
    labs(x = "Number of proteins", y = NULL, title = "Distribution across protein classes") +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.text.y = element_blank(), strip.background = element_blank(),
          strip.text.x = element_text(face = "bold", size = 12), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), panel.spacing.y = unit(0.35, "lines"), plot.title.position = "plot",
          plot.title = element_text(size = 13))
  list(rings = rings, classes = classes)
}

#' Draw the one-page overview on the open device.
draw_qc_overview <- function(plots, title = "Olink Explore HT: proteins above and below LOD", no_class_note = NULL, family = "") {
  plots <- lapply(plots, set_family, fam = family)
  fg <- if (nzchar(family)) list(fontfamily = family) else list()
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(0.06, 0.36, 0.58), "npc"))))
  grid.text(title, x = 0.02, just = "left", gp = do.call(gpar, c(list(fontsize = 16, fontface = "bold"), fg)), vp = viewport(layout.pos.row = 1))
  # layout cell first, then the shifted viewport (x / width are ignored when given together with layout.pos.row)
  pushViewport(viewport(layout.pos.row = 2))
  print(plots$rings, vp = viewport(x = if (is.null(plots$classes)) 0.5 else 0.62, width = 0.76))
  popViewport()
  if (!is.null(plots$classes)) print(plots$classes, vp = viewport(layout.pos.row = 3))
  else if (!is.null(no_class_note)) grid.text(no_class_note, vp = viewport(layout.pos.row = 3), gp = do.call(gpar, c(list(fontsize = 11, col = "grey30"), fg)))
  popViewport()
}
