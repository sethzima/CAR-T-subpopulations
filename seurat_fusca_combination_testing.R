#install.packages("devtools")
library(devtools)
devtools::install_github("edroaldo/fusca")
library(dplyr)
library(fusca)
library(igraph)
library(scales)
library(ggplot2)
library(Seurat)
library(patchwork)
library(dplyr)
library(hdf5r)
library(doParallel)


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

# VlnPlot(WT.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(St.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(Us.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

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


ElbowPlot(WT.combined)
ElbowPlot(St.combined)
ElbowPlot(Us.combined)


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

#marker genes from supplemental info
features.wt <- c("FTL", "S100A8", "CD14", "GNLY", "NKG7", "CD19", "CD79A", "CD3D", "CD8A", "CD4")

my.features.wt <- c("CD14", "LYZ", "MS4A7", "CST3", "ACSL1", "AGTRAP", "AIF1", "FCN1", "GNS")

b.features.wt <- c("MS4A1", "CD79A", "CD19")

lymph.features.wt <- c("IL7R", "CCR7", "S100A4", "CD8A", "CD4", "GNLY", "NKG7", "CD48", "CD96")

t.features.wt <- c("CD2", "CD3D", "CD3E", "CD3G", "CD48", "CD8A", "IL7R", "IL16", "IL23A")



#visulatizations
FeaturePlot(WT.combined, features = b.features.wt)

FeaturePlot(WT.combined, features = my.features.wt)

FeaturePlot(WT.combined, features = lymph.features.wt)

FeaturePlot(WT.combined, features = t.features.wt)


#find cluster 5 markers
wt.cluster5.markers <- FindMarkers(WT.combined, ident.1 = 5, min.pct = 0.25)
write.csv(wt.cluster5.markers, file="wt_cluster5_markers_seurat.csv", row.names=T)


#annotated WT clusters
wt.cluster.ids <- c("T cells", "T cells", "B lymphocytes", "T cells", "T cells", "T cells",
                     "NK", "B lymphocytes", "Myeloid/monocyte", "T cells", "T cells", "Myeloid/monocyte")

names(wt.cluster.ids) <- levels(WT.combined)

WT.combined <- RenameIdents(WT.combined, wt.cluster.ids)

DimPlot(WT.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
