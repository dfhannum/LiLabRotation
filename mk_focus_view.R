library(Seurat)

wbm <- readRDS('data/wbm_clustered_filtered_named.rds')

wbm.list <- SplitObject(wbm, split.by = 'cluster_IDs')

mk <- wbm.list$MK
mk@meta.data
dim(mk@meta.data)
dim(mk)

summary(as.factor(mk@meta.data$state))
head(mk@meta.data)


mk<- FindVariableFeatures(mk, selection.method = 'vst', nfeatures = 1000)
top10 <- head(VariableFeatures(mk), 10)
top10

all.genes <- row.names(mk)

mk <- ScaleData(mk, features = all.genes)

mk <- RunPCA(mk, features = VariableFeatures(object = mk))

DimPlot(mk, reduction = 'pca', group.by = 'state')

ElbowPlot(mk)

mk <-FindNeighbors(mk, dims = 1:10)

lst_of_resolutions <- seq(0,1, by = .05)

resul = c()
cnt = 1
for (i in lst_of_resolutions){
        #print(i)
        x <- FindClusters(mk, resolution = i)
        resul[cnt] <- length(levels(x$seurat_clusters))
        cnt = cnt + 1
}

plot(resul)
lst_of_resolutions[12]

mk <- FindClusters(mk, resolution = 0.55)
mk <- RunUMAP(mk, dims = 1:10)

DimPlot(mk, reduction = 'umap', label = T) + NoLegend()
DimPlot(mk, reduction = 'umap', group.by = 'state')
DimPlot(mk, reduction = 'umap', label = T, split.by = 'state')
