install.packages("devtools")
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



#clustering
#wild type
cellrouter.wt <- computePCA(cellrouter.wt, assay.type = "RNA", 
                         seed = 42, num.pcs = 50, genes.use = cellrouter.wt@var.genes) 
plot(cellrouter.wt@pca$sdev)

# cellrouter.wt <- computeTSNE(cellrouter.wt, 
#                           seed = 1, num.pcs = 8, max_iter = 1000)

cellrouter.wt <- computeUMAP(cellrouter.wt, 
                          seed = 1, num.pcs = 8, metric = "cosine", n_neighbors = 30, spread = 1, min_dist = 0.1)

cellrouter.wt <- findClusters(cellrouter.wt, assay.type = "RNA", 
                           num.pcs = 8, nn.type = "knn", k = 20)

# plotReducedDimension(cellrouter.wt, assay.type = "RNA", reduction.type = "tsne", annotation = "population", annotation.color = "colors",
#                      dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)


plotReducedDimension(cellrouter.wt, assay.type = "RNA", reduction.type = "umap", annotation = "population", annotation.color = 'colors',
                     dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)





#stimulated
cellrouter.st <- computePCA(cellrouter.st, assay.type = "RNA", 
                            seed = 42, num.pcs = 50, genes.use = cellrouter.st@var.genes) 
plot(cellrouter.st@pca$sdev)

# cellrouter.wt <- computeTSNE(cellrouter.wt, 
#                           seed = 1, num.pcs = 8, max_iter = 1000)

cellrouter.st <- computeUMAP(cellrouter.st, 
                             seed = 1, num.pcs = 8, metric = "cosine", n_neighbors = 30, spread = 1, min_dist = 0.1)

cellrouter.st <- findClusters(cellrouter.st, assay.type = "RNA", 
                              num.pcs = 8, nn.type = "knn", k = 20)

# plotReducedDimension(cellrouter.wt, assay.type = "RNA", reduction.type = "tsne", annotation = "population", annotation.color = "colors",
#                      dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)


plotReducedDimension(cellrouter.st, assay.type = "RNA", reduction.type = "umap", annotation = "population", annotation.color = 'colors',
                     dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)


#unstimulated
cellrouter.us <- computePCA(cellrouter.us, assay.type = "RNA", 
                            seed = 42, num.pcs = 50, genes.use = cellrouter.us@var.genes) 
plot(cellrouter.us@pca$sdev)

# cellrouter.wt <- computeTSNE(cellrouter.wt, 
#                           seed = 1, num.pcs = 8, max_iter = 1000)

cellrouter.us <- computeUMAP(cellrouter.us, 
                             seed = 1, num.pcs = 8, metric = "cosine", n_neighbors = 30, spread = 1, min_dist = 0.1)

cellrouter.us <- findClusters(cellrouter.us, assay.type = "RNA", 
                              num.pcs = 8, nn.type = "knn", k = 20)

# plotReducedDimension(cellrouter.wt, assay.type = "RNA", reduction.type = "tsne", annotation = "population", annotation.color = "colors",
#                      dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)


plotReducedDimension(cellrouter.us, assay.type = "RNA", reduction.type = "umap", annotation = "population", annotation.color = 'colors',
                     dotsize = 1.5, showlabels = TRUE, labelsize = 5, convex = FALSE)


#cluster-specific gene signatures and marker genes
#wt clustering
markers.wt <- findSignatures(cellrouter.wt, assay.type = "RNA", column = "population", pos.only = TRUE, 
                             fc.threshold = 0.2, nCores = 10)

markers.st <- findSignatures(cellrouter.st, assay.type = "RNA", column = "population", pos.only = TRUE, 
                             fc.threshold = 0.2, nCores = 10)
  
markers.us <- findSignatures(cellrouter.us, assay.type = "RNA", column = "population", pos.only = TRUE, 
                             fc.threshold = 0.2, nCores = 10)


markers.wt.down <-  findSignatures(cellrouter.wt, assay.type = "RNA", column = "population", pos.only = FALSE, 
                                   fc.threshold = 0.2, nCores = 10)




write.csv(markers.wt, file="wt_clusters_fusca.csv", row.names=F)
write.csv(markers.st, file="st_clusters_fusca.csv", row.names=F)
write.csv(markers.us, file="us_clusters_fusca.csv", row.names=F)



#cluster annotations
#WT
tmp.wt <- recode(as.character(cellrouter.wt@assays$RNA@sampTab$population), 
              "1"="CD4+ T cell", 
              "2"="CD4+ T cell",
              "3"="Dendritic cell",
              "4"="CD33+ Myeloid",
              "5"="Macrophage",
              "6"="CD14+ monocyte",
              "7"="721 B Lymphocyte",
              "8"="CD8+ T cell",
              "9"="721 B Lymphoctye",
              "10"="CD19+ B cell",
              "11"="Double negative T",
              "12"="Intermediate B cells",
              "13"="CD56+ NK cells"
              )
names(tmp.wt) <- rownames(cellrouter.wt@assays$RNA@sampTab)

df <- data.frame(cluster = cellrouter.wt@assays$RNA@sampTab$population, celltype = tmp.wt)
rownames(df) <- rownames(cellrouter.wt@assays$RNA@sampTab)

cellrouter.wt <- addInfo(cellrouter.wt, assay.type = "RNA", sample.name = sample_names[[1]],
                      metadata = df, colname = "celltype",
                      metadata.column = "celltype")

plotReducedDimension(cellrouter.wt, assay.type = 'RNA', 
                     reduction.type = 'umap', annotation="celltype", annotation.color = 'celltype_color',
                     dotsize=1, showlabels = TRUE, labelsize=5, convex=FALSE)



#calculate mean expression of ligands and receptors per cluster
#wt
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))# from NicheNet
head(lr_network)

pairs <- lr_network
pairs$Pair.Name <- paste(pairs$from, pairs$to, sep = "_")

ligands <- unique(lr_network$from)
ligands <- intersect(ligands, rownames(cellrouter.wt@assays$RNA@ndata))

receptors <- unique(lr_network$to)
receptors <- intersect(receptors, rownames(cellrouter.wt@assays$RNA@ndata))

ligands.receptors <- unique(c(ligands, receptors))

mean.expr <- computeValue(cellrouter.wt, assay.type = "RNA", 
                          genelist = ligands.receptors, "celltype", fun = "mean"); gc();

interactions <- population.pairing(mean.expr = mean.expr, pairs = pairs, ligands = ligands, receptors = receptors, threshold = 0.25)

interactions <- calculateObservedMean(mean.expr = mean.expr, interactions = interactions)
head(interactions)

markers <- findSignatures(cellrouter.wt, assay.type = "RNA", 
                          column = "celltype", pos.only = TRUE, fc.threshold = 0.2, nCores = 10)


#calculate null distribution of intracluster gene expression means
genelist <- unique(c(interactions$ligand, interactions$receptor))

p <- clusterPermutation(cellrouter.wt, assay.type = "RNA", 
                        genelist = genelist, interactions = interactions, cluster.label = "celltype", nPerm = 1000, nCores = 10)

interactions.p <- calculatePvalue(p, nPerm = 1000, interactions2 = interactions)

tmp.wt <- interactions.p[which(interactions.p$pvalue < 0.01),]
head(tmp.wt)


#calculate a network based on the number of interactions (ligand/receptors), 
#cluster centroids and distances between these clusters
my_matrix <- interactionmatrix(tmp)
head(matrix)

my_graph <- cellnetwork3(tmp, threshold = 5)


cellrouter.wt <- calculateCentroids(cellrouter.wt, assay.type = "RNA", sample.name = "Sample1", 
                                 cluster.column = "celltype", cluster.type = "Cluster")

cellrouter.wt <- calculateDistanceMatrix(cellrouter.wt, assay.type = "RNA", sample.name = "Sample1", 
                                      cluster.type = "Cluster", spot_distance = 100, normalize = FALSE)

plots <- predictCellInteractions(cellrouter.wt, assay.type = "RNA", sample.name = "Sample1", 
                                 cluster.type = "Cluster", graph = my_graph, distance.threshold = 0.75)
options(repr.plot.width = 22, repr.plot.height = 7)
gridExtra::grid.arrange(grobs = plots, ncol = 1)

