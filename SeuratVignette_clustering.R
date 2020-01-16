library(Seurat)
library(dplyr)

pbmc.data <- Read10X(data.dir = 'data/seurat_vignette/filtered_gene_bc_matrices/hg19/')
pbmc <- CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, min.features = 200)
?CreateSeuratObject

pbmc

# looking into this object

class(pbmc)
class(pbmc.data)
rownames(pbmc.data)

'CHD7' %in% rownames(pbmc.data)

pbmc.data[c(1,2,3),1:20]

#QC Metrics

head(pbmc@meta.data,5)

# adding percentage mitochondrial RNA

pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

head(pbmc@meta.data,5)
dim(pbmc@meta.data)

library(ggplot2)
ggplot(pbmc@meta.data, aes(x = nFeature_RNA)) + geom_density()
ggplot(pbmc@meta.data, aes(x = nFeature_RNA)) + geom_density() + 
        xlim(xmin = 1500, xmax = 3500)

ggplot(pbmc@meta.data, aes(x = nCount_RNA)) + geom_density()
ggplot(pbmc@meta.data, aes(x = nCount_RNA)) + geom_density() +
        xlim(xmin = 5000, xmax = 16000)

ggplot(pbmc@meta.data, aes(x = percent.mt)) + geom_density()
ggplot(pbmc@meta.data, aes(x = percent.mt)) + geom_density() + 
        xlim(xmin = 5, xmax = max(pbmc@meta.data$percent.mt))

# xlim isn't the best way to do this since it only looks at hte density of that region, rather
# I would just like to look at that region

# This is the subset that they use
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalizing the data

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

?NormalizeData

# Finding variable features, they return 2000 by default

?FindVariableFeatures

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
top10

# Plot to visual these points

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

# Next scaling the data set before doing PCs

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
?ScaleData

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
?RunPCA
print(pbmc[['pca']], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')
DimPlot(pbmc, reduction = 'pca')
pbmc[['pca']]

DimHeatmap(pbmc, dims = 1:9, cells = 500, balanced = T)

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:20)

?JackStraw

ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
?FindNeighbors



