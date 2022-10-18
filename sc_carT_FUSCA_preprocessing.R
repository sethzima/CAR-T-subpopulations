list.of.packages <- c('cccd', 'grid', 'tsne', 'Rtsne', 'igraph', 'mclust', 'ggplot2', 'pheatmap', 'reshape', 'reshape2')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages, repos = c("http://cran.rstudio.com/", "https://bioconductor.org/biocLite.R"))

# In case of a missing package after running the steps above, please install them by typing:
#source('http://bioconductor.org/biocLite.R')
#biocLite('package_name')

# The package Vennerable also needs to be installed:
install.packages("Vennerable", repos = "http://R-Forge.R-project.org")

#required to fix error in multiprocessing
library(doParallel)
library(hdf5r)
# Install FUSCA
library(devtools)
#devtools::install_github("edroaldo/fusca")

list.of.packages<-c("cccd", "grid", "Rcpp", "tsne", "uwot", "dplyr", "fusca", "Rtsne", "igraph", "Matrix", "mclust", "tibble", 
                    "cowplot", "ggplot2", "reshape", "reshape2", "pheatmap")
lapply(list.of.packages, require, character.only = TRUE)

set.seed(42)
setwd


#cellrouter.het <- CreateCellRouter(data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)