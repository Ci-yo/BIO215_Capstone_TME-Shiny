#!/usr/bin/env Rscript
# Install the packages loaded by app.R. Run once from the repository root.
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("shiny", "shinythemes", "DT", "dplyr", "ggplot2",
              "randomForest", "shinydashboard", "httr", "jsonlite")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
message("TME Explorer dependencies are ready.")
