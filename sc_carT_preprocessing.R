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

# top10.wt <- head(VariableFeatures(WT.combined), 10)
# top10.st <- head(VariableFeatures(St.combined), 10)
# top10.us <- head(VariableFeatures(Us.combined), 10)
# 
# plot.wt <- VariableFeaturePlot(WT.combined)
# plot.wt2 <- LabelPoints(plot = plot.wt, points = top10.wt, repel = TRUE)
# plot.wt + plot.wt2
# 
# 
# plot.st <- VariableFeaturePlot(St.combined)
# plot.st2 <- LabelPoints(plot = plot.st, points = top10.st, repel = TRUE)
# plot.st + plot.st2
# 
# plot.us <- VariableFeaturePlot(Us.combined)
# plot.us2 <- LabelPoints(plot = plot.us, points = top10.us, repel = TRUE)
# plot.us + plot.us2

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


#differential gene expression
WT.cluster.markers <- FindAllMarkers(WT.combined)
St.cluster.markers <- FindAllMarkers(St.combined)
Us.cluster.markers <- FindAllMarkers(Us.combined)

write.csv(WT.cluster.markers, file="wt_clusters.csv", row.names=F)
write.csv(St.cluster.markers, file="st_clusters.csv", row.names=F)
write.csv(Us.cluster.markers, file="us_clusters.csv", row.names=F)

WT.cluster.narrowed <- WT.cluster.markers %>%
                          group_by(cluster) %>%
                          slice_max(n = 20, order_by = avg_log2FC)

St.cluster.narrowed <- St.cluster.markers %>%
                          group_by(cluster) %>%
                          slice_max(n = 20, order_by = avg_log2FC)

Us.cluster.narrowed <- Us.cluster.markers %>%
                          group_by(cluster) %>%
                          slice_max(n=20, order_by = avg_log2FC)

WT.genes <-  WT.cluster.narrowed[ , 7]

wt.cluster0 <- WT.genes[1:20, ]
wt.cluster1 <- WT.genes[21:40, ]
wt.cluster2 <- WT.genes[41:60, ]
wt.cluster3 <- WT.genes[61:80, ]
wt.cluster4 <- WT.genes[81:100, ]
wt.cluster5 <- WT.genes[101:120, ]
wt.cluster6 <- WT.genes[121:140, ]
wt.cluster7 <- WT.genes[141:160, ]
wt.cluster8 <- WT.genes[161:180, ]
wt.cluster9 <- WT.genes[181:200, ]
wt.cluster10 <- WT.genes[201:220, ]
wt.cluster11 <- WT.genes[221:240, ]


St.genes <- St.cluster.narrowed[ , 7]

St.cluster0 <- St.genes[1:20, ]
St.cluster1 <- St.genes[21:40, ]
St.cluster2 <- St.genes[41:60, ]
St.cluster3 <- St.genes[61:80, ]
St.cluster4 <- St.genes[81:100, ]
St.cluster5 <- St.genes[101:120, ]
St.cluster6 <- St.genes[121:140, ]
St.cluster7 <- St.genes[141:160, ]
St.cluster8 <- St.genes[161:180, ]
St.cluster9 <- St.genes[181:200, ]
St.cluster10 <- St.genes[201:220, ]
St.cluster11 <- St.genes[221:240, ]
St.cluster12 <- St.genes[241:260, ]


Us.genes <- Us.cluster.narrowed[ , 7]
Us.cluster0 <- Us.genes[1:20, ]
Us.cluster1 <- Us.genes[21:40, ]
Us.cluster2 <- Us.genes[41:60, ]
Us.cluster3 <- Us.genes[61:80, ]
Us.cluster4 <- Us.genes[81:100, ]
Us.cluster5 <- Us.genes[101:120, ]
Us.cluster6 <- Us.genes[121:140, ]
Us.cluster7 <- Us.genes[141:160, ]
Us.cluster8 <- Us.genes[161:180, ]
Us.cluster9 <- Us.genes[181:200, ]
Us.cluster10 <- Us.genes[201:220, ]
Us.cluster11 <- Us.genes[221:240, ]


#in order to do more, I'd like to have the clusters annotated with their cell types

#cluster annotation
wt.new.cluster.ids <- c("CD4+ T cell", "CD4+ naive/CD8+ memeory effector T", "CD19+ B Cells", "CD56+ NK cells", 
                        "CD4+ Naive T cells", "CD8+ T cells", "CD56+ NK cells", "CD19+ B cells", "CD33+ Myeloid",
                        "T lymphocyte", "CD56+ NK cells", "Myeloid Dendritic Cell")
names(wt.new.cluster.ids) <- levels(WT.combined)
WT.combined <- RenameIdents(WT.combined, wt.new.cluster.ids)
DimPlot(WT.combined, reduction = "umap", label=TRUE)


st.new.cluster.ids <- c("CD8+ T cells", "CD4+ central memory T cells", "CD8+ Effectory memory T cells", "Dendritic cells",
                        "NK cells", "CD8+ Effectory memory T cells", "Proliferating NK/T cells", "CD4+ Central memory T cells",
                        "CD4+ central memory T cells", "B lymphoblasts", "CD56+ NK cells", "Dendritic cells", "CD4+ T cells")
names(st.new.cluster.ids) <- levels(St.combined)
St.combined <- RenameIdents(St.combined, st.new.cluster.ids)
DimPlot(St.combined, reduction = "umap", label=TRUE)

us.new.cluster.ids <- c("CD4+ memory T", "Naive CD4+ T", "B lymphocytes", "CD56+ NK", "CD8+ T", "Dendritic",
                       "Proliferating NK/T", "Naive regulatory T", "CD4+ proliferating T", "Dendritic", "Proliferating macrophage", "CD4+ T")
names(us.new.cluster.ids) <- levels(Us.combined)
Us.combined <- RenameIdents(Us.combined, us.new.cluster.ids)
DimPlot(Us.combined, reduction = "umap", label=TRUE)