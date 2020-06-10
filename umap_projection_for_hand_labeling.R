# BIG IMAGE OF THE SEURAT CLUSTERS FOR ANNOTATION

library(Seurat)
library(ggplot2)

wbm <- readRDS('./data/wbm_clustered_filtered_named.rds')

head(wbm@meta.data)

new_cluster_ids <- 0:12
names(new_cluster_ids) <- levels(wbm)
new_cluster_ids
wbm <- RenameIdents(wbm, new_cluster_ids)

DimPlot(wbm, reduction = 'umap', label = T)
ggsave('./figures/large_umap_projection_for_annotation.pdf', device = 'pdf',
       units = 'in', width = 10, height = 7.5)
