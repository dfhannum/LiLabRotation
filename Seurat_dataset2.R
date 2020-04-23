# Looking at the second set of experiments

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot)

exp2 <- Read10X(data.dir= 'data/Experiment2/filtered_feature_bc_matrix/')

# Init the seurat object with the raw (non-normalized data)

exp2 <- CreateSeuratObject(counts = exp2, project = 'exp2', min.cells = 3, min.features = 200)

head(exp2@meta.data)

exp2[['percent.mt']] <- PercentageFeatureSet(exp2, pattern = '^mt')
head(exp2@meta.data)

VlnPlot(exp2, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

ggplot(data = exp2@meta.data, aes(x = percent.mt)) +
        geom_density() + theme_bw()

summary(exp2$percent.mt > 10)

ggplot(exp2@meta.data, aes(x = nFeature_RNA)) + geom_density() + theme_bw()

summary(exp2$nFeature_RNA > 5000)

ggplot(exp2@meta.data, aes(x = nCount_RNA)) + geom_density() + theme_bw()

# filtering

dim(exp2)
exp2 <- subset(exp2, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
dim(exp2)

# going from 5,947 cells to 4,584

# other things to temper the data

exp2 <- NormalizeData(exp2, normalization.method = 'LogNormalize', scale.factor = 10000)

exp2 <- FindVariableFeatures(exp2, selection.method = 'vst', nFeatures = 2000)

top10 <- head(VariableFeatures(exp2),10)
top10

plt1 <- VariableFeaturePlot(exp2)
LabelPoints(plot = plt1, points = top10, xnudge = 0, ynudge = 0, repel = T)

all.genes <- rownames(exp2)
exp2 <- ScaleData(exp2, features = all.genes)
exp2 <- RunPCA(exp2, features = VariableFeatures(object = exp2))

DimPlot(exp2, reduction = 'pca')

ElbowPlot(exp2)

exp2 <- FindNeighbors(exp2, dims = 1:10)

# looking at different resolutions

lst_of_res <- seq(0,1,by = 0.05)

resul = c()
cnt = 1
for (i in lst_of_res){
        x <- FindClusters(exp2, resolution = i)
        resul[cnt] <- length(levels(x$seurat_clusters))
        cnt = cnt + 1
}

plot(resul)

lst_of_res[10]

exp2 <- FindClusters(exp2, resolution = 0.45)

exp2 <- RunUMAP(exp2, dims = 1:10)

DimPlot(exp2, reduction = 'umap', label = T)
ggsave('./figures/EXP2_umap_projection.png', device = 'png',
       units = 'in', height = 5, width = 5, dpi = 400)

saveRDS(exp2, file = './data/EXP2_clustered_filtered.rds')
