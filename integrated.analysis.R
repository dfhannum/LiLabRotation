# combining the two experiments
library(data.table)
library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
# creating a single Seurat object

wbm1 <- readRDS('./data/wbm_clustered_filtered_named.rds')
wbm2 <- readRDS('./data/EXP2_clustered_filtered_SingleR_labels.rds')

wbm.comb <- merge(wbm1,wbm2, add.cell.ids = c('exp1','exp2'), project = 'WBM_comb',
                  merge.data = T)

wbm.comb
head(colnames(wbm.comb), 20)

summary(as.factor(tstrsplit(colnames(wbm.comb),'_', keep = 1)[[1]]))

head(wbm.comb@meta.data)
tail(wbm.comb@meta.data)

wbm.comb@meta.data$orig.ident <- ifelse(wbm.comb@meta.data$orig.ident == 'exp2','exp2','exp1')

summary(as.factor(wbm.comb@meta.data$orig.ident))

wbm.list <- SplitObject(wbm.comb, split.by = 'orig.ident')
wbm.list <- wbm.list[c('exp1','exp2')]

for (i in 1:length(wbm.list)){
        wbm.list[[i]] <- NormalizeData(wbm.list[[i]], verbose = F)
        wbm.list[[i]] <- FindVariableFeatures(wbm.list[[i]], selection.method = 'vst',
                                              nfeatures = 2000, verbose = F)
}

ref.list <- wbm.list
wbm.anchors <- FindIntegrationAnchors(object.list = ref.list, dims = 1:30)
wbm.int <- IntegrateData(anchorset = wbm.anchors, dims = 1:30)

DefaultAssay(wbm.int) <- 'integrated'

wbm.int <- ScaleData(wbm.int, verbose = F)
wbm.int <- RunPCA(wbm.int, npcs = 30, verbose = F)
ElbowPlot(wbm.int)
wbm.int <- RunUMAP(wbm.int, reduction = 'pca', dims = 1:30)

tail(wbm.int@meta.data)

wbm.int@meta.data$cluster_name <- ifelse(is.na(wbm.int@meta.data$cluster_IDs), wbm.int@meta.data$cluster,
                                    wbm.int@meta.data$cluster_IDs)
                                   
summary(as.factor(wbm.int@meta.data$cluster_name))
length(unique(wbm.int@meta.data$cluster_name))

p1 <- DimPlot(wbm.int, reduction = 'umap', group.by = 'orig.ident')
p2 <- DimPlot(wbm.int, reduction = 'umap', group.by = 'cluster_name', label = T, repel = T,
        cols = DiscretePalette(29)) + 
        NoLegend()
p1 + p2
ggsave('./figures/INT_umap.png', device = 'png', units = 'in', height = 4, width = 9.5, dpi = 400)


#++++++
# Getting the new clusters

head(wbm.int@meta.data)

DimPlot(wbm.int, reduction = 'umap', group.by = 'seurat_clusters')

wbm.int <- FindNeighbors(wbm.int, dims = 1:30)

lst_ <- seq(0,1, by = 0.05)

resul = c()
cnt = 1

for (i in lst_){
        x <- FindClusters(wbm.int, resolution = i)
        resul[cnt] <- length(levels(x$seurat_clusters))
        cnt = cnt + 1
}

plot(resul)
lst_[9]
wbm.int <- FindClusters(wbm.int, resolution = .4)

DimPlot(wbm.int, reduction = 'umap', label = T)

saveRDS(wbm.int,'./data/wbm.integrated.RDS')

