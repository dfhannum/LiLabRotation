# EDA on experiment one

library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)


wbm <- Read10X(data.dir = './data/filtered_feature_bc_matrix/')

wbm <- CreateSeuratObject(counts = wbm, project = 'wbm_1', min.cells = 3, min.features = 200)

wbm

wbm[['percent.mt']] <- PercentageFeatureSet(wbm, pattern = '^mt')
head(wbm@meta.data)

htos <- read.table('./data/HTOs.csv', sep = ',', header = T)
head(htos)

htos$Barcode <- tstrsplit(htos$Barcode,'-', keep = 1)[[1]]

summary(tstrsplit(htos$Barcode,'-', keep = 1)[[1]] %in% rownames(wbm@meta.data))
# already filtered out some when creating the seurat object

summary(htos$HTOs)

htos$HTOs <- as.factor(ifelse(htos$HTOs == '', 'unknown', as.character(htos$HTOs)))
summary(htos)

htos <- htos[htos$Barcode %in% rownames(wbm@meta.data),]

summary(htos$Barcode == rownames(wbm@meta.data))

HTOS <- ifelse(htos$HTOs == 'HTO-1', 'mut',
       ifelse(htos$HTOs == 'HTO-2', 'enr_mut',
              ifelse(htos$HTOs == 'HTO-3', 'cnt',
                     ifelse(htos$HTOs == 'HTO-4', 'enr_cnt', 'unknown'))))

wbm[['HTO']] <- HTOS

head(wbm@meta.data)

wbm <- subset(wbm, subset = HTO != 'unknown')

wbm <- NormalizeData(wbm, normalization.method = 'LogNormalize', scale.factor = 10000)

wbm <- FindVariableFeatures(wbm, selection.method = 'vst', nfeatures = 2000)
all.genes <- rownames(wbm)
wbm <- ScaleData(wbm, features = all.genes)

wbm <- RunPCA(wbm, features = VariableFeatures(object = wbm))
ElbowPlot(wbm)

wbm <- FindNeighbors(wbm, dims = 1:10)
wbm <- FindClusters(wbm, resolution = 0.45)
wbm <- RunUMAP(wbm, dims = 1:10)

DimPlot(wbm, reduction = 'umap', cols = DiscretePalette(14))
ggsave('./figures/eda_plots/FILTERING_umap.png', device = 'png',
       units = 'in', height = 4, width = 7)
library(patchwork)

## Want to get some labels on these cells first so it is more informative

library(SingleR)
library(scRNAseq)

m.ref <- ImmGenData()
m.ref2 <- MouseRNAseqData()

SCwbm <- as.SingleCellExperiment(wbm)

pred_cluster <- SingleR(test = SCwbm, ref = list(m.ref,m.ref2),
                        labels = list(m.ref$label.main, m.ref2$label.main),
                        method = 'cluster',
                        clusters = SCwbm$seurat_clusters)

cluster_ids <- pred_cluster[,4]

names(cluster_ids) <- levels(wbm)
cluster_ids

wbm <- RenameIdents(wbm, cluster_ids)

wbm@meta.data$cluster_ids <- mapvalues(wbm@meta.data$seurat_clusters, 
                                       from = levels(wbm@meta.data$seurat_clusters),
                                       to = cluster_ids)


DimPlot(wbm, reduction = 'umap')
ggsave('./figures/eda_plots/FILTERING_umap_NAMED.png', device = 'png',
       units = 'in', height = 4, width = 7)

head(wbm@meta.data)

## Back to look at the comparisons

#percent.mt
ggplot(data = wbm@meta.data, aes(x = seurat_clusters, y = nFeature_RNA)) + geom_violin() +
        theme_bw() + geom_hline(yintercept = c(500,5000), colour = 'red', linetype = 2)
ggsave('./figures/eda_plots/FILTERING_nFeature_violin_plot.png', device = 'png',
       units = 'in', height = 4, width = 7)

ggplot(data = wbm@meta.data, aes(x = seurat_clusters, y = percent.mt)) + geom_violin() +
        theme_bw() + geom_hline(yintercept = 10, colour = 'red', linetype = 2)
ggsave('./figures/eda_plots/FILTERING_percent.mt_violin_plot.png', device = 'png',
       units = 'in', height = 4, width = 7)

#nFeatures
ggplot(data = wbm@meta.data, aes(x = cluster_ids, y = nFeature_RNA)) + geom_violin() +
        theme_bw() + geom_hline(yintercept = c(500,5000), colour = 'red', linetype = 2)

ggsave('./figures/eda_plots/FILTERING_nFeature_violin_plot_NAMED_CLUSTERS.png', device = 'png',
       units = 'in', height = 4, width = 7)

ggplot(data = wbm@meta.data, aes(x = cluster_ids, y = percent.mt)) + geom_violin() +
        theme_bw() + geom_hline(yintercept = 10, colour = 'red', linetype = 2)
ggsave('./figures/eda_plots/FILTERING_percent.mt_violin_plot_NAMED_CLUSTERS.png', device = 'png',
       units = 'in', height = 4, width = 7)


wbm[['filtered_out_by_nFeature_RNA']] <- as.factor(ifelse(wbm@meta.data$nFeature_RNA < 500, 'too low',
                                                ifelse(wbm@meta.data$nFeature_RNA > 5000, 'too high', 'keep')))

DimPlot(wbm, reduction = 'umap', group.by = 'filtered_out_by_nFeature_RNA',
        cols = c('grey','pink','darkred'))

wbm[['mt_filter']] <- as.factor(ifelse(wbm@meta.data$percent.mt > 10, 'toss','keep'))

DimPlot(wbm, reduction = 'umap', group.by = 'mt_filter',
        cols = c('grey','red'))

wbm[['toss_status']] <- ifelse(wbm@meta.data$filtered_out_by_nFeature_RNA != 'keep' &
                                                 wbm@meta.data$mt_filter =='toss', 'double toss (155)',
                                         ifelse(wbm@meta.data$mt_filter == 'toss', 'mt.toss (228)',
                                                ifelse (wbm@meta.data$filtered_out_by_nFeature_RNA == 'keep',
                                                        'keep (6138)', 'nFeature.toss (707)')))

summary(as.factor(wbm$toss_status))



DimPlot(wbm, reduction = 'umap', group.by = 'toss_status',
        cols = c('black','grey','purple','red'))
ggsave('./figures/eda_plots/FILTERING_umap_toss.png',  device = 'png',
       units = 'in', height = 4, width = 7)
DimPlot(wbm, reduction = 'umap', split.by = 'toss_status')
ggsave('./figures/eda_plots/FILTERING_umap_toss_split.png', device = 'png',
       units = 'in', height = 5, width = 10)
