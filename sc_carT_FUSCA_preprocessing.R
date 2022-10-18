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


library(Seurat)
library(patchwork)
library(dplyr)
library(hdf5r)
#C:\Users\Patron\Documents\MyDocuments\R\CAR-T\CAR-T-subpopulations\GSE145809_RAW

# files <- list.files()
# 
# for (i in 1:length(files)) {
#   fileName <- files
#   file_new <- paste0("Sample", 1:length(files))
#   
#   sample <- Read10x()
#   
# }
#data = "C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/WT/"

#read in the data
set.seed(42)

WTsample.1 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/WT/GSM4333345_RBG9802_filtered_gene_bc_matrices_h5.h5")
WTsample.2 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/WT/GSM4333346_RBG13746_filtered_gene_bc_matrices_h5.h5")
WTsample.3 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/WT/GSM4333347_RBG13747_filtered_gene_bc_matrices_h5.h5")

Stim.1 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/Stimulated/GSM4333351_RBG13751_filtered_gene_bc_matrices_h5.h5")
Stim.2 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/Stimulated/GSM4333352_RBG13752_filtered_gene_bc_matrices_h5.h5")
Stim.3 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/Stimulated/GSM4333353_RBG13753_filtered_gene_bc_matrices_h5.h5")

Ustim.1 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/Unstimulated/GSM4333348_RBG13748_filtered_gene_bc_matrices_h5.h5")
Ustim.2 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/Unstimulated/GSM4333349_RBG13749_filtered_gene_bc_matrices_h5.h5")
Ustim.3 <- Read10X_h5("C:/Users/Patron/Documents/MyDocuments/R/CAR-T/CAR-T-subpopulations/GSE145809_RAW/Unstimulated/GSM4333350_RBG13750_filtered_gene_bc_matrices_h5.h5")


WT1 <- CreateSeuratObject(counts = WTsample.1, assay = "RNA", project = "WT T-cells")
WT2 <- CreateSeuratObject(counts = WTsample.2, assay = "RNA", project = "WT T-cells")
WT3 <- CreateSeuratObject(counts = WTsample.3, assay = "RNA", project = "WT T-cells")

St1 <- CreateSeuratObject(counts = Stim.1, assay = "RNA", project = "Stimulated T-cells")
St2 <- CreateSeuratObject(counts = Stim.2, assay = "RNA", project = "Stimulated T-cells")
St3 <- CreateSeuratObject(counts = Stim.3, assay = "RNA", project = "Stimulated T-cells")

Us1 <- CreateSeuratObject(counts = Ustim.1, assay = "RNA", project = "Unstimulated T-cells")
Us2 <- CreateSeuratObject(counts = Ustim.2, assay = "RNA", project = "Unstimulated T-cells")
Us3 <- CreateSeuratObject(counts = Ustim.3, assay = "RNA", project = "Unstimulated T-cells")

WT.combined <- merge(WT1, y = c(WT2, WT3), add.cell.ids = c("Donor1", "Donor2", "Donor3"), project = "WT T-cells")
St.combined <- merge(St1, y = c(St2, St3), add.cell.ids = c("Donor1", "Donor2", "Donor3"), project = "Stimulated T-cells")
Us.combined <- merge(Us1, y = c(Us2, Us3), add.cell.ids = c("Donor1", "Donor2", "Donor3"), project = "Unstimulated T-cells")

wt.data <- GetAssayData(WT.combined, slot = "counts")
st.data <- GetAssayData(St.combined, slot = "counts")
us.data <- GetAssayData(Us.combined, slot = "counts")


#quality control
#wt
cellrouter.wt <- CreateCellRouter(wt.data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)
cellrouter.st <- CreateCellRouter(st.data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)
cellrouter.us <- CreateCellRouter(us.data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)
mito.genes.wt <- grep(pattern = "^MT-", x = rownames(x = cellrouter.wt@assays$RNA@ndata), value = TRUE)
percent.mito.wt <- Matrix::colSums(cellrouter.wt@assays$RNA@ndata[mito.genes.wt, ]) / Matrix::colSums(cellrouter.wt@assays$RNA@ndata)
ribo.genes.wt <- grep(pattern = "^RP[SL]", x = rownames(x = cellrouter.wt@assays$RNA@ndata), value = TRUE)
percent.ribo.wt <- Matrix::colSums(cellrouter.wt@assays$RNA@ndata[ribo.genes.wt, ]) / Matrix::colSums(cellrouter.wt@assays$RNA@ndata)

cellrouter.wt@assays$RNA@sampTab$percent.mito.wt = percent.mito.wt
cellrouter.wt@assays$RNA@sampTab$percent.ribo.wt = percent.ribo.wt

p1 = ggplot(cellrouter.wt@assays$RNA@sampTab, aes(x = "percent.mito", y = percent.mito.wt)) + 
  geom_violin(fill = "grey80", colour = "#FF0000") + theme(legend.position = "none") + xlab("") 
p2 = ggplot(cellrouter.wt@assays$RNA@sampTab, aes(x = "percent.ribo", y = percent.ribo.wt)) + 
  geom_violin(fill = "grey80", colour = "#FF9900") + theme(legend.position = "none") + xlab("") 
p3 = ggplot(cellrouter.wt@assays$RNA@sampTab, aes(x = "nGene", y = nGene)) + 
  geom_violin(fill = "grey80", colour = "#3366FF") + theme(legend.position = "none") + xlab("") 
plot_grid(p1, p2, p3, nrow = 1) # labels = "")


#st
mito.genes.st <- grep(pattern = "^MT-", x = rownames(x = cellrouter.st@assays$RNA@ndata), value = TRUE)
percent.mito.st <- Matrix::colSums(cellrouter.st@assays$RNA@ndata[mito.genes.st, ]) / Matrix::colSums(cellrouter.st@assays$RNA@ndata)
ribo.genes.st <- grep(pattern = "^RP[SL]", x = rownames(x = cellrouter.st@assays$RNA@ndata), value = TRUE)
percent.ribo.st <- Matrix::colSums(cellrouter.st@assays$RNA@ndata[ribo.genes.st, ]) / Matrix::colSums(cellrouter.st@assays$RNA@ndata)

cellrouter.st@assays$RNA@sampTab$percent.mito.st = percent.mito.st
cellrouter.st@assays$RNA@sampTab$percent.ribo.st = percent.ribo.st

p1 = ggplot(cellrouter.st@assays$RNA@sampTab, aes(x = "percent.mito", y = percent.mito.st)) + 
  geom_violin(fill = "grey80", colour = "#FF0000") + theme(legend.position = "none") + xlab("") 
p2 = ggplot(cellrouter.st@assays$RNA@sampTab, aes(x = "percent.ribo", y = percent.ribo.st)) + 
  geom_violin(fill = "grey80", colour = "#FF9900") + theme(legend.position = "none") + xlab("") 
p3 = ggplot(cellrouter.st@assays$RNA@sampTab, aes(x = "nGene", y = nGene)) + 
  geom_violin(fill = "grey80", colour = "#3366FF") + theme(legend.position = "none") + xlab("") 
plot_grid(p1, p2, p3, nrow = 1) # labels = "")


#us
mito.genes.us <- grep(pattern = "^MT-", x = rownames(x = cellrouter.us@assays$RNA@ndata), value = TRUE)
percent.mito.us <- Matrix::colSums(cellrouter.us@assays$RNA@ndata[mito.genes.us, ]) / Matrix::colSums(cellrouter.us@assays$RNA@ndata)
ribo.genes.us <- grep(pattern = "^RP[SL]", x = rownames(x = cellrouter.us@assays$RNA@ndata), value = TRUE)
percent.ribo.us <- Matrix::colSums(cellrouter.us@assays$RNA@ndata[ribo.genes.us, ]) / Matrix::colSums(cellrouter.us@assays$RNA@ndata)

cellrouter.us@assays$RNA@sampTab$percent.mito.us = percent.mito.us
cellrouter.us@assays$RNA@sampTab$percent.ribo.us = percent.ribo.us

p1 = ggplot(cellrouter.us@assays$RNA@sampTab, aes(x = "percent.mito", y = percent.mito.us)) + 
  geom_violin(fill = "grey80", colour = "#FF0000") + theme(legend.position = "none") + xlab("") 
p2 = ggplot(cellrouter.us@assays$RNA@sampTab, aes(x = "percent.ribo", y = percent.ribo.us)) + 
  geom_violin(fill = "grey80", colour = "#FF9900") + theme(legend.position = "none") + xlab("") 
p3 = ggplot(cellrouter.us@assays$RNA@sampTab, aes(x = "nGene", y = nGene)) + 
  geom_violin(fill = "grey80", colour = "#3366FF") + theme(legend.position = "none") + xlab("") 
plot_grid(p1, p2, p3, nrow = 1) # labels = "")

#filtering
cellrouter.wt <- filterCells(cellrouter.wt, assay.type = "RNA", variables = c("nGene", "percent.mito.wt"),
                             thresholds.low = c(400, -Inf), thresholds.high = c(10000, 0.08))

cellrouter.st <- filterCells(cellrouter.st, assay.type = "RNA", variables = c("nGene", "percent.mito.st"),
                             thresholds.low = c(1000, -Inf), thresholds.high = c(30000, 0.08))

cellrouter.us <- filterCells(cellrouter.us, assay.type = "RNA", variables = c("nGene", "percent.mito.us"),
                                thresholds.low = c(1000, -Inf), thresholds.high = c(30000, 0.08))

