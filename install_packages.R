# Install the R packages used by the pipeline (R >= 4.3).
cran <- c("OlinkAnalyze", "arrow", "dplyr", "tidyr", "readr", "readxl", "stringr", "purrr", "tibble",
          "ggplot2", "yaml", "writexl", "lme4", "lmerTest", "rmcorr", "msigdbr", "BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install(c("limma", "variancePartition", "fgsea", "BiocParallel"), update = FALSE, ask = FALSE)
