# 03 - Exploratory overview per matrix
# PCA coloured by biology and plate; variance explained by subject, skin state, plate, cohort.
# Out: output/explore/*

source("R/utils.R")
source("R/qc.R")
cfg  <- load_config()
meta <- read_step(cfg, "metadata", "sample_metadata.rds", step = "scripts/01_metadata.R")
wide <- read_step(cfg, "data", "npx_wide.rds", step = "scripts/02_import_qc.R")

for (mx in names(wide)) {
  m <- wide[[mx]]
  if (is.null(m) || !nrow(m)) next
  sc <- pca_scores(m)
  ve <- attr(sc, "var_explained")
  sc <- sc |> left_join(meta, by = "SampleID")
  colour_by <- if (mx == "ISF") c("state", "group", "visit", "plate") else c("group", "cohort", "plate", "relapse")
  for (v in colour_by) {
    p <- ggplot(sc, aes(PC1, PC2, colour = .data[[v]])) + geom_point(size = 1.8, alpha = 0.8) +
      labs(title = sprintf("%s - PCA coloured by %s", mx, v),
           x = sprintf("PC1 (%.0f%%)", ve[1]), y = sprintf("PC2 (%.0f%%)", ve[2]))
    save_plot(p, cfg, "explore", sprintf("pca_%s_%s.png", mx, v), width = 7, height = 5)
  }

  # variance partitioning (random effects for all terms, as recommended by variancePartition)
  info <- meta |> filter(SampleID %in% colnames(m)) |> as.data.frame()
  rownames(info) <- info$SampleID
  if (mx == "ISF") {
    info$state <- coalesce(info$state, "none")
    form <- ~ (1 | SubjectID) + (1 | state) + (1 | plate)
  } else {
    form <- ~ (1 | group) + (1 | cohort) + (1 | plate)
  }
  info <- info[colnames(m), ]
  vp <- suppressMessages(suppressWarnings(
    variancePartition::fitExtractVarPartModel(m, form, info, BPPARAM = bpparam_cores(), quiet = TRUE)))
  vp_df <- as.data.frame(vp) |> rownames_to_column("OlinkID")
  save_csv(vp_df, cfg, "explore", sprintf("variance_partition_%s.csv", mx))
  p <- vp_df |> pivot_longer(-OlinkID, names_to = "term", values_to = "fraction") |>
    ggplot(aes(reorder(term, -fraction, median), fraction)) + geom_violin(fill = "grey85") +
    geom_boxplot(width = 0.1, outlier.size = 0.3) +
    labs(title = sprintf("%s - share of variance per protein", mx), x = NULL, y = "fraction of variance")
  save_plot(p, cfg, "explore", sprintf("variance_partition_%s.png", mx), width = 7, height = 5)
  msg("%s: median variance share - %s", mx,
      paste(sprintf("%s %.2f", names(vp), apply(vp, 2, median)), collapse = ", "))
}
