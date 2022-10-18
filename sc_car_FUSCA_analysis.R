install.packages("devtools")
library("devtools")
devtools::install_github("edroaldo/fusca")
library("dplyr")
library("fusca")
library("igraph")
library("scales")
library("ggplot2")
library(Seurat)
library(patchwork)
library(dplyr)
library(hdf5r)

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


cellrouter.wt <- CreateCellRouter(wt.data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)
cellrouter.st <- CreateCellRouter(st.data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)
cellrouter.us <- CreateCellRouter(us.data, assay.type = "RNA", min.genes = 200, min.cells = 3, is.expr = 0)


#preprocessing
var.genes.wt <- FindVariableGenes(cellrouter.wt, assay.type = "RNA", 
                               method = "coefficient_variation", pvalue = 0.05)

var.genes.wt <- var.genes.wt[order(var.genes.wt$adj.pvalue, decreasing = FALSE), ]
cellrouter.wt@var.genes <- rownames(var.genes.wt[1:2000, ])

cellrouter.wt <- Normalize(cellrouter.wt, assay.type = "RNA")

cellrouter.wt <- scaleData(cellrouter.wt, assay.type = "RNA", 
                        genes.use = cellrouter.wt@var.genes, blocksize = nrow(cellrouter.wt@assays$RNA@ndata))





var.genes.st <- FindVariableGenes(cellrouter.st, assay.type = "RNA", 
                                  method = "coefficient_variation", pvalue = 0.05)

var.genes.st <- var.genes.st[order(var.genes.st$adj.pvalue, decreasing = FALSE), ]
cellrouter.st@var.genes <- rownames(var.genes.st[1:2000, ])

cellrouter.st <- Normalize(cellrouter.st, assay.type = "RNA")

cellrouter.st <- scaleData(cellrouter.st, assay.type = "RNA", 
                           genes.use = cellrouter.st@var.genes, blocksize = nrow(cellrouter.st@assays$RNA@ndata))






var.genes.us <- FindVariableGenes(cellrouter.us, assay.type = "RNA", 
                                  method = "coefficient_variation", pvalue = 0.05)

var.genes.us <- var.genes.us[order(var.genes.us$adj.pvalue, decreasing = FALSE), ]
cellrouter.us@var.genes <- rownames(var.genes.us[1:2000, ])

cellrouter.us <- Normalize(cellrouter.us, assay.type = "RNA")

cellrouter.us <- scaleData(cellrouter.us, assay.type = "RNA", 
                           genes.use = cellrouter.us@var.genes, blocksize = nrow(cellrouter.us@assays$RNA@ndata))




