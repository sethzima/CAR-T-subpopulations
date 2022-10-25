#install.packages("devtools")
library(devtools)
#devtools::install_github("edroaldo/fusca")
#library(fusca)
library(igraph)
library(scales)
library(ggplot2)
library(Seurat)
library(patchwork)
library(dplyr)
library(hdf5r)
library(doParallel)

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
options(timeout = max(300, getOption("timeout")))

if (!requireNamespace("renv")) install.packages("renv")
library(renv)
renv::restore()

remotes::install_github("carmonalab/scGate")
remotes::install_github("carmonalab/ProjecTILs")

library(ProjecTILs)

#looking at reference data
ref <- load.reference.map()
refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000",
             "#87f6a5", "#e812dd")
DimPlot(ref, label = T, cols = refCols)

markers <- c("Cd4", "Cd8a", "Ccr7", "Tcf7", "Pdcd1", "Havcr2", "Tox", "Izumo1r",
             "Cxcr6", "Xcl1", "Gzmb", "Gzmk", "Ifng", "Foxp3")
VlnPlot(ref, features = markers, stack = T, flip = T, assay = "RNA")


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


#preprocessing

#MT genes
WT.combined[["percent.mt"]] <- PercentageFeatureSet(WT.combined, pattern = "^MT-")
St.combined[["percent.mt"]] <- PercentageFeatureSet(St.combined, pattern = "^MT-")
Us.combined[["percent.mt"]] <- PercentageFeatureSet(Us.combined, pattern = "^MT-")

#QC metrics
#at least 400 genes for leuka, 1000 for product
#less than 8% mito
#fileter out more than 10000 and 30000

WT.combined <- subset(WT.combined, subset = nFeature_RNA > 400 & nCount_RNA < 10000 & percent.mt < 8)
St.combined <- subset(St.combined, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 8)
Us.combined <- subset(Us.combined, subset = nFeature_RNA > 1000 & nCount_RNA < 30000 & percent.mt < 8)


#log normalizing
WT.combined <- NormalizeData(WT.combined, normalization.method = "LogNormalize", scale.factor = 10000)
St.combined <- NormalizeData(St.combined, normalization.method = "LogNormalize", scale.factor = 10000)
Us.combined <- NormalizeData(Us.combined, normalization.method = "LogNormalize", scale.factor = 10000)

#finding highly variable genes
WT.combined <- FindVariableFeatures(WT.combined, selection.method = "vst", nfeatures = 2000)
St.combined <- FindVariableFeatures(St.combined, selection.method = "vst", nfeatures = 2000)
Us.combined <- FindVariableFeatures(Us.combined, selection.method = "vst", nfeatures = 2000)

#scaling the data
all.genes.wt <- rownames(WT.combined)
all.genes.st <- rownames(St.combined)
all.genes.us <- rownames(Us.combined)

WT.combined <- ScaleData(WT.combined, features = all.genes.wt)
St.combined <- ScaleData(St.combined, features = all.genes.st)
Us.combined <- ScaleData(Us.combined, features = all.genes.us)


#run PCA
WT.combined <- RunPCA(WT.combined)
St.combined <- RunPCA(St.combined)
Us.combined <- RunPCA(Us.combined)


#clustering
WT.combined <- FindNeighbors(WT.combined, dims = 1:50)
WT.combined <- FindClusters(WT.combined, resolution = 0.5)

St.combined <- FindNeighbors(St.combined, dims = 1:50)
St.combined <- FindClusters(St.combined, resolution = 0.5)

Us.combined <- FindNeighbors(Us.combined, dims = 1:50)
Us.combined <- FindClusters(Us.combined, resolution = 0.5)

#UMAP
WT.combined <- RunUMAP(WT.combined, dims = 1:50)
St.combined <- RunUMAP(St.combined, dims = 1:50)
Us.combined <- RunUMAP(Us.combined, dims = 1:50)

DimPlot(WT.combined, reduction = "umap", label = TRUE)
DimPlot(St.combined, reduction = "umap", label = TRUE)
DimPlot(Us.combined, reduction = "umap", label = TRUE)


#project the query datasets onto the reference atlas
wt.projected <- Run.ProjecTILs(WT.combined, ref = ref, filter.cells = FALSE)
plot.projection(ref, wt.projected)
plot.statepred.composition(ref, wt.projected, metric = "Percent")


st.projected <- Run.ProjecTILs(St.combined, ref = ref, filter.cells = FALSE)
plot.projection(ref, st.projected)
plot.statepred.composition(ref, st.projected, metric = "Percent")

us.projected <- Run.ProjecTILs(Us.combined, ref = ref, filter.cells = FALSE)
plot.projection(ref, us.projected)
plot.statepred.composition(ref, us.projected, metric = "Percent")

plot.states.radar(ref, query = wt.projected)
plot.states.radar(ref, query = st.projected)
plot.states.radar(ref, query = us.projected)

#comparing gene programs
remotes::install_github("carmonalab/SignatuR")
library(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)

library(UCell)
ref <- AddModuleScore_UCell(ref, features = programs, assay = "RNA", name = NULL)

#wild type
wt.projected <- AddModuleScore_UCell(wt.projected, features = programs, assay = "RNA",
                                        name = NULL)
plot.states.radar(ref, query = wt.projected, meta4radar = names(programs))

#st
st.projected <- AddModuleScore_UCell(st.projected, features = programs, assay = "RNA",
                                     name = NULL)
plot.states.radar(ref, query = st.projected, meta4radar = names(programs))

#us
us.projected <- AddModuleScore_UCell(us.projected, features = programs, assay = "RNA",
                                     name = NULL)
plot.states.radar(ref, query = us.projected, meta4radar = names(programs))



#compare states across conditions
query.control <- subset(us.projected)
query.perturb <- subset(st.projected)

plot.states.radar(ref, query = list(Control = query.control, Query = query.perturb))


#
