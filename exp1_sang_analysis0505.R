# Updated analysis for exp1

# have a meeting with the Sang lab today and wanted to create some potential
# images for the manuscript. Looking specifically at Hay 2018 and recreating
# some of their figures and looking at their marker genes in our exp1 dataset.
# Ideally run this analysis again for exp2, and for the integrated analysis, but
# still waiting on HTO tags.

library(Seurat)
library(ggplot2)

wbm <- readRDS('./data/wbm_clustered_filtered_named.rds')

head(wbm@meta.data)

DimPlot(wbm, reduction = 'umap')

#features from the paper


features <- c('Avp','Spink2','Negr1','Cygb','Mdp-2','Fcer1a')
features <- features[features %in% rownames(wbm)]
features

VlnPlot(wbm, features = features)

FeaturePlot(wbm, features = features)

DotPlot(wbm, features = features) + RotatedAxis()

DoHeatmap(subset(wbm, downsample = 100), features = features, size = 3)

features <- c('Fcer1a','Itga2b','Ank1','Csf3r','Sox4','Gata3','Top2a','Mki67')

FeaturePlot(wbm, features = features)

VlnPlot(wbm, features = features[1:3], pt.size = 0)

DotPlot(wbm, features = features) + RotatedAxis()

DoHeatmap(subset(wbm, downsample = 100), features = features, size = 3)

feat2 <- c('nFeature_RNA','nCount_RNA','percent.mt')

VlnPlot(wbm, features = feat2, ncol =3, pt.size = 0 )
FeaturePlot(wbm, features = feat2)

?VlnPlot


# Figure generating 
