library(Seurat)

wbm.int <- readRDS('./data/wbm.integrated.RDS')

DimPlot(wbm.int, reduction = 'umap', split.by = 'orig.ident', label = T,
        repel = T) + NoLegend()


head(wbm.int@meta.data)

wbm.list <- SplitObject(wbm.int, split.by = 'integrated_snn_res.0.4')

iMK <- wbm.list$`8`

summary(iMK@meta.data)

summary(as.factor(iMK$cluster_name))

iMK <- FindVariableFeatures(iMK, selection.method = 'vst', nfeatures = 100)
head(VariableFeatures(iMK),10)

all.genes <- row.names(iMK)

iMK <- ScaleData(iMK, features = all.genes)
iMK <- RunPCA(iMK, features = VariableFeatures(object = iMK))

ElbowPlot(iMK)
iMK <- FindNeighbors(iMK, dims = 1:20)

lst_ = seq(0,1,.05)
resul = c()
cnt = 1
for (i in lst_){
        x <- FindClusters(iMK, resolution = i)
        resul[cnt] <- length(levels(x$seurat_clusters))
        cnt = cnt + 1
}
plot(resul)
lst_[15]

iMK <- FindClusters(iMK, resolution = 0.9)
iMK <- RunUMAP(iMK, dims = 1:20)

DimPlot(iMK, reduction = 'umap', label = T) + NoLegend()
DimPlot(iMK, reduction = 'umap', split.by = 'orig.ident')
DimPlot(iMK, reduction = 'umap', split.by = 'state')
