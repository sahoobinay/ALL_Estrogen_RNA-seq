# install_packages.R
# Run this script once to install all required packages before running pipeline.R

message("Installing CRAN packages...")
cran_packages <- c(
  "tidyverse",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "patchwork",
  "reshape2",
  "scales",
  "UpSetR"
)

install_if_missing <- function(pkgs, installer) {
  to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    installer(to_install)
  } else {
    message("All packages already installed.")
  }
}

install_if_missing(cran_packages,
  function(pkgs) install.packages(pkgs, repos = "https://cloud.r-project.org"))

message("Installing Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

bioc_packages <- c(
  "edgeR",
  "limma",
  "fgsea",
  "msigdbr",
  "viper",
  "dorothea",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

install_if_missing(bioc_packages,
  function(pkgs) BiocManager::install(pkgs, ask = FALSE))

message("Installing igraph (CRAN)...")
install_if_missing("igraph",
  function(pkgs) install.packages(pkgs, repos = "https://cloud.r-project.org"))

message("All packages installed. Run pipeline.R to start the analysis.")
