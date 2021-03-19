library(Seurat)
library(cowplot)
library(ggplot2)

## DATA IMPORT LOCATIONS
gex_matrix <- "H:\\shared_data\\IFN_scSeq\\Run_179-BP_IFN-06-08\\Sample_IFN_pool_08\\outs\\filtered_feature_bc_matrix"
hto_matrix <- "H:/shared_data/SM_Pilot/combined_ADT-HTO/umi_count"

# **************************************************

# Load in UMIs
pbmc.umis <- Read10X(gex_matrix)
pbmc.umis

# Load in HTOs
pbmc.htos <- Read10X(hto_matrix, gene.column=1)
pbmc.htos
rownames(pbmc.htos)

# Remove unmapped reads
pbmc.htos <- pbmc.htos[setdiff(rownames(x = pbmc.htos), "unmapped"), ]
pbmc.htos
rownames(pbmc.htos)

# Select cell barcodes detected by both RNA and HTO
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

# ***************************************************************

## INITIALIZE SEURAT OBJECT

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)

# Find and scale variable featuresot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)

# ***************************************************************

## PROCESS SEURAT OBJECT

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using the
# default settings
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.85) ## Adjust this to control stringecy
pbmc.hashtag

# Global classification results
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot", verbose = FALSE)
table(pbmc.hashtag$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:3], ncol = 1)

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
FeatureScatter(pbmc.hashtag, feature1 = "HTO1", feature2 = "HTO3")

# Compare number of UMIs for singlets, doublets and negative cells
Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.hashtag[["percent.mt"]] <- PercentageFeatureSet(pbmc.hashtag, pattern = "^mt-")
VlnPlot(pbmc.hashtag, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0.1, ncol = 3)

# Filter based on cell quality metrics
pbmc.hashtag
pbmc.hashtag <- subset(pbmc.hashtag, subset = nFeature_RNA > 500 & nFeature_RNA < 7500 & nCount_RNA > 100 & nCount_RNA < 60000 & percent.mt < 15) ## Adjust these to set dead cell parameters
pbmc.hashtag

# Run PCA/UMAP on most variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))
pbmc.hashtag <- RunPCA(pbmc.hashtag, verbose = FALSE)
pbmc.hashtag <- RunTSNE(pbmc.hashtag, dims = 1:5, perplexity = 100)

Idents(pbmc.hashtag) <- "HTO_classification" # More granular classification
DimPlot(pbmc.hashtag, ncol = 3, split.by = "HTO_classification")
Idents(pbmc.hashtag) <- "HTO_classification.global" # Summary classification
DimPlot(pbmc.hashtag, ncol = 2, split.by = "HTO_classification.global")

FeaturePlot(pbmc.hashtag, features = c("CD3D", "CD33", "MT-ATP6", "GNLY"), min.cutoff = "q5", max.cutoff = "q99", ncol = 2, split.by = "HTO_classification.global") ## Substitute genes of interest as QC step
FeaturePlot(pbmc.hashtag, features = c("CD3D", "CD33", "MT-ATP6", "GNLY"), min.cutoff = "q5", max.cutoff = "q99", ncol = 2, split.by = "HTO_maxID")

# Remove negative cells from the object
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)
pbmc.hashtag.subset

# Run PCA on most variable features
pbmc.hashtag.subset <- FindVariableFeatures(pbmc.hashtag.subset, selection.method = "mean.var.plot")
pbmc.hashtag.subset <- ScaleData(pbmc.hashtag.subset, features = VariableFeatures(pbmc.hashtag.subset))
pbmc.hashtag.subset <- RunPCA(pbmc.hashtag.subset, verbose = FALSE)
pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, dims = 1:5, perplexity = 100)

Idents(pbmc.hashtag.subset) <- "HTO_classification" # More granular classification
DimPlot(pbmc.hashtag.subset, ncol = 3, split.by = "HTO_classification")
Idents(pbmc.hashtag.subset) <- "HTO_classification.global" # Summary classification
DimPlot(pbmc.hashtag.subset, ncol = 2, split.by = "HTO_classification.global")

# View HTO heatmap (to increase plotting efficiency, use num.cells argument)
HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000)

# Extract the singlets
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")
pbmc.singlet

# Select the top 1000 most variable features
pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")
# Scaling RNA data (only scale the variable features here for efficiency)
pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
# Run PCA
pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))
ElbowPlot(pbmc.singlet)
# Use top 10 PCs for clustering and tSNE based on PCElbowPlot
pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = "pca", dims = 1:13)
pbmc.singlet <- FindClusters(pbmc.singlet, resolution = 0.8)
pbmc.singlet <- RunUMAP(pbmc.singlet, reduction = "pca", dims = 1:13)

# Project singlet identities on TSNE visualization
DimPlot(pbmc.singlet, group.by = "HTO_classification")
DimPlot(pbmc.singlet, split.by = "HTO_classification", ncol = 2)

RidgePlot(pbmc.singlet, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:4], ncol = 2)

# Subset into individual HTOS
Idents(pbmc.singlet) <- "HTO_classification"

pbmc.hto1 <- subset(pbmc.singlet, idents = "HTO1")
pbmc.hto1

pbmc.hto2 <- subset(pbmc.singlet, idents = "HTO2")
pbmc.hto2

pbmc.hto3 <- subset(pbmc.singlet, idents = "HTO3")
pbmc.hto3

pbmc.hto4 <- subset(pbmc.singlet, idents = "HTO4")
pbmc.hto4

# Save out objects for each HTO
out_folder <- "H:\\shared_data\\IFN_scSeq\\IFN Seurat Data\\"
saveRDS(pbmc.hto1, file = (paste(out_folder,"XX.rds", sep="")))
saveRDS(pbmc.hto2, file = (paste(out_folder,"XX.rds", sep="")))
saveRDS(pbmc.hto3, file = (paste(out_folder,"XX.rds", sep="")))
saveRDS(pbmc.hto4, file = (paste(out_folder,"XX.rds", sep="")))

# Save global singlet object as archive
saveRDS(pbmc.singlet, file = (paste(out_folder,"XX.rds", sep="")))