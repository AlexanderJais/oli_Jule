# Helpers for the PDF executive summary (base grid + ggplot only; no extra packages).

library(grid)

pdf_open <- function(file) grDevices::pdf(file, width = 11.69, height = 8.27, onefile = TRUE)   # A4 landscape

page_header <- function(title, subtitle = NULL) {
  grid.text(title, x = 0.04, y = 0.95, just = "left", gp = gpar(fontsize = 18, fontface = "bold"))
  if (!is.null(subtitle)) grid.text(subtitle, x = 0.04, y = 0.91, just = "left", gp = gpar(fontsize = 10, col = "grey35"))
  grid.lines(x = c(0.04, 0.96), y = c(0.89, 0.89), gp = gpar(col = "grey70"))
}

#' Text page: items is a character vector; each item becomes a wrapped bullet
#' (items starting with "## " are sub-headings).
page_text <- function(title, items, subtitle = NULL, width = 125, size = 10.5) {
  grid.newpage(); page_header(title, subtitle)
  y <- 0.86; line_h <- size / 72 / 8.27 * 1.35
  for (it in items) {
    head <- startsWith(it, "## ")
    txt <- if (head) sub("^## ", "", it) else strwrap(it, width = width, exdent = 2, prefix = "", initial = "- ")
    if (y - length(txt) * line_h < 0.04) {                  # continue on a new page
      grid.newpage(); page_header(paste(title, "(continued)"), subtitle); y <- 0.86
    }
    if (head) y <- y - 0.4 * line_h
    for (l in txt) {
      grid.text(l, x = 0.05, y = y, just = c("left", "top"),
                gp = gpar(fontsize = if (head) size + 1.5 else size, fontface = if (head) "bold" else "plain"))
      y <- y - line_h
    }
    y <- y - 0.35 * line_h
  }
}

#' Table page: a data frame drawn in a monospace font (long tables continue on further pages).
page_table <- function(title, df, subtitle = NULL, note = NULL, rows_per_page = 30, digits = 3) {
  if (is.null(df) || !nrow(df)) { page_text(title, "No results available.", subtitle); return(invisible()) }
  df <- as.data.frame(df)
  for (j in seq_along(df)) {
    if (is.numeric(df[[j]])) {
      x <- df[[j]]
      df[[j]] <- if (all(is.na(x) | x == round(x))) ifelse(is.na(x), "", format(round(x), big.mark = "", trim = TRUE))
                 else ifelse(is.na(x), "", vapply(x, \(v) format(signif(v, digits), scientific = abs(v) < 1e-3 && v != 0), ""))
    }
    df[[j]] <- str_trunc(as.character(df[[j]] %||% ""), 38)
    df[[j]][is.na(df[[j]])] <- ""
  }
  w <- pmax(nchar(names(df)), vapply(df, \(x) max(nchar(x), 0), numeric(1)))
  fmt <- \(v) paste(mapply(\(x, n) formatC(x, width = -n), v, w), collapse = "  ")
  chunks <- split(seq_len(nrow(df)), ceiling(seq_len(nrow(df)) / rows_per_page))
  total_w <- sum(w) + 2 * length(w)
  fs <- max(5.5, min(9.5, 9.5 * 150 / max(total_w, 150)))
  for (k in seq_along(chunks)) {
    grid.newpage(); page_header(if (k == 1) title else paste(title, "(continued)"), subtitle)
    lines <- c(fmt(names(df)), strrep("-", total_w), apply(df[chunks[[k]], , drop = FALSE], 1, fmt))
    for (i in seq_along(lines))
      grid.text(lines[i], x = 0.04, y = 0.86 - (i - 1) * 0.026, just = c("left", "top"),
                gp = gpar(fontfamily = "mono", fontsize = fs, fontface = if (i == 1) "bold" else "plain"))
    if (!is.null(note) && k == length(chunks))
      grid.text(note, x = 0.04, y = 0.03, just = "left", gp = gpar(fontsize = 8.5, col = "grey35"))
  }
}

page_plot <- function(p, title = NULL) {
  if (is.null(p)) return(invisible())
  if (!is.null(title)) p <- p + labs(title = title)
  print(p)
}

#' Run one section; if it fails, write a note page instead of stopping the report.
section <- function(name, expr) {
  tryCatch(expr, error = \(e) page_text(name, paste("This section could not be produced:", conditionMessage(e))))
}

#' Read an output CSV if it exists (NULL otherwise).
out_csv <- function(cfg, ...) {
  p <- file.path(cfg$paths$output, ...)
  if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL
}

top_names <- function(df, n = 10) if (is.null(df) || !nrow(df)) "none" else paste(head(unique(df$Assay), n), collapse = ", ")
