# Install the R packages used by the LEIP galanin study (R >= 4.3).
cran <- c("OlinkAnalyze", "arrow", "dplyr", "tidyr", "readr", "readxl", "stringr", "stringi", "purrr", "tibble",
          "ggplot2", "yaml", "writexl", "msigdbr", "BiocManager")
install.packages(setdiff(cran, rownames(installed.packages())))
BiocManager::install("fgsea", update = FALSE, ask = FALSE)   # gene-set enrichment (question 2)
