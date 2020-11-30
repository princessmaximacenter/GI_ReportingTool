# R script containing the setup information for the running the genetic
# interaction reporting tool. 
# All packages are installed by Docker. This setup will check if the 
# installation is correctly done and otherwise missing packages will be installed.
#
# Author(s): Denise Kersjes
# Date of creation:  11  September 2020
# Date of last edit: 30  November 2020

# Specify the CRAN packages
packages = c("shiny", "cgdsr", "mclust", "stringr", "XML", "seqinr", 
             "dplyr", "tidyr")

# Avoid error message about missing mirror site 
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

# Load the package when available, otherwise install the package
lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE, quietly = TRUE)) {
    install.packages(x, dependencies = TRUE, quitely = TRUE)
    library(x, character.only = TRUE, quietly = TRUE)
  }
})

# Install and load BiocManager packages with BiocManager
biocmanager.packages = c("data.table", "maftools")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

lapply(biocmanager.packages, FUN = function(x) {
  if (!require(x, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(x)
    library(x, character.only = TRUE, quietly = TRUE)
  }
})
