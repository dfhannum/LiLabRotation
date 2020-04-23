library(Seurat)
library(SingleR)
library(ggplot2)

wbm_exp <- readRDS('data/wbm_exp.rds')

head(wbm_exp@meta.data)

# creating a state variable combining celltype and condition for later analysis
wbm_exp@meta.data$state <- paste0(wbm_exp@meta.data$celltype, '-',
                                  wbm_exp@meta.data$condition)

## Splitting of the unknown state cells

wbm.list <- SplitObject(wbm_exp, split.by = 'assigned_celltype')

wbm <- wbm.list['Yes']$Yes
wbm

head(wbm@meta.data)

## Looking at how to filter the dataset
# dashed lines represent the cutoffs being used

summary(wbm@meta.data)

ggplot(wbm@meta.data, aes(x = nFeature_RNA, colour = state)) + geom_density() +
        geom_vline(xintercept = c(500,5000), colour = 'red', linetype = 'dashed')
ggsave('./figures/nFeature_RNA.png', device = 'png', units = 'in',
       width = 6, height = 4, dpi = 300)

ggplot(wbm@meta.data, aes(x = percent.mt, colour = state)) + geom_density() +
        geom_vline(xintercept = 10, colour = 'red', linetype = 'dashed')
ggsave('./figures/percent_mt.png', device = 'png', units = 'in',
       width = 6, height = 4, dpi = 300)

ggplot(wbm@meta.data, aes(x = nCount_RNA, colour = state)) + geom_density()
ggsave('./figures/nCount_RNA.png', device = 'png', units = 'in',
       width = 6, height = 4, dpi = 300)

# counts of the data/state before subsetting
summary(as.factor(wbm@meta.data$state))

wbm <- subset(wbm, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

# counts of the data per state after subsetting
summary(as.factor(wbm@meta.data$state))

## Running analysis to get the clusters and visualizaiton for the data

wbm <- NormalizeData(wbm, normalization.method = 'LogNormalize', scale.factor = 10000)

wbm <- FindVariableFeatures(wbm, selection.method = 'vst', nfeatures = 2000)

top10 <- head(VariableFeatures(wbm), 10)

plot1 <- VariableFeaturePlot(wbm)
LabelPoints(plot = plot1, points = top10, xnudge = 0, ynudge = 0,repel = T)

all.genes <- row.names(wbm)
wbm <- ScaleData(wbm, features = all.genes)

wbm <- RunPCA(wbm, features = VariableFeatures(object = wbm))

DimPlot(wbm, reduction = 'pca')
DimPlot(wbm, reduction = 'pca', group.by = 'state')
# seems like the second PCA is distinguishing the enr_WBM-mut from the others

ElbowPlot(wbm)
# going to use the first 10 PCs

wbm <- FindNeighbors(wbm, dims = 1:10)

## Looking to see how the resolution affects the number of clusters
lst_of_resolutions <- seq(0,1, by = .05)

resul = c()
cnt = 1
for (i in lst_of_resolutions){
        #print(i)
        x <- FindClusters(wbm, resolution = i)
        resul[cnt] <- length(levels(x$seurat_clusters))
        cnt = cnt + 1
}

plot(resul)

# going with resolution 10
lst_of_resolutions[10]

wbm <- FindClusters(wbm, resolution = 0.45)

wbm <- RunUMAP(wbm, dims = 1:10)
DimPlot(wbm, reduction = 'umap', label = T)
ggsave('./figures/umap_projection.png', device = 'png', units = 'in',
       height = 5, width = 6, dpi = 400)
DimPlot(wbm, reduction = 'umap', label = T, split.by = 'state')
ggsave('./figures/umap_projection_split_by_state.png', device = 'png',
       units = 'in', height = 5, width = 12, dpi = 400)

## Getting tables to show the distribution of cells between states and clusters

 
# Now normalizing over the rows
tbl_rowPerc <- tbl

for (i in 1:dim(tbl)[1]){
        for (j in 1:dim(tbl)[2]){
                tbl_rowPerc[i,j] <- round(tbl[i,j] / rowSums(tbl)[i],2) * 100
        }
}

tbl1 <- as.matrix(tbl_rowPerc)

heatmap(tbl1)
cluster_order <- paste0('Cluster ',c(4,3,0,1,2,9,8,7,12,11,10,6,5))

tbl_rowPerc <- reshape(tbl_rowPerc, idvar = 'state', ids = row.names(tbl_rowPerc), 
                       times = names(tbl_rowPerc),  varying = list(names(tbl_rowPerc)),
                       direction = 'long')

rownames(tbl_rowPerc) <- NULL


tbl_rowPerc$time <- factor(tbl_rowPerc$time, levels = cluster_order)

library(plyr)

tbl_rowPerc$state <- as.factor(tbl_rowPerc$state)
new_states <- mapvalues(tbl_rowPerc$state, from = levels(tbl_rowPerc$state), 
                        to = paste0(levels(tbl_rowPerc$state),' (n=', rowSums(tbl),')'))
tbl_rowPerc$state <- new_states

tbl_rowPerc$time2 <- rep(c('Granulocyte1','Granulocyte2','Granulocyte3','B-cell','MK','Granulocyte4',
      'Stem Cell', 'Monocyte','Macrophage','Erythroid','B-cell pro.','T-cell/NK','MEP'),each = 4)

ggplot(data = tbl_rowPerc, aes(x = time2, y = state, fill = `Cluster 0`)) + geom_tile() +
        scale_fill_gradient2(high = 'darkred', low = 'white', mid = 'red', midpoint = 50,
                             limit = c(0,100), space ="Lab",
                             name = 'Percentage of Cells\nper Cluster ') +
        theme_minimal() + ylab('Experimental States') + xlab('') +
        ggtitle ("Heatmap of Cell Distributions within Experimental States") + 
        theme (plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_text(aes(time2, state, label = `Cluster 0`), color = 'black', size = 4)+
        coord_fixed()

ggsave('./figures/NAMED_heatmap_cells_within_states.png', device = 'png', units = 'in',
       height = 5, width = 8, dpi = 400)

# normalizing over columns
tbl_colPerc <- tbl

for (i in 1:dim(tbl)[1]){
        for (j in 1:dim(tbl)[2]){
                tbl_colPerc[i,j] <- round(tbl[i,j] / colSums(tbl)[j],2) *100
        }
}

tbl_colPerc <- reshape(tbl_colPerc, idvar = 'state', ids = row.names(tbl_colPerc), 
                       times = names(tbl_colPerc),  varying = list(names(tbl_colPerc)),
                       direction = 'long')

rownames(tbl_rowPerc) <- NULL

tbl_colPerc$time <- factor(tbl_colPerc$time, levels = cluster_order)

library(plyr)

tbl_colPerc$state <- as.factor(tbl_colPerc$state)
new_states <- mapvalues(tbl_colPerc$state, from = levels(tbl_colPerc$state), 
                        to = paste0(levels(tbl_colPerc$state),' (n=', rowSums(tbl),')'))
tbl_colPerc$state <- new_states
tbl_colPerc$time2 <- rep(c('Granulocyte1','Granulocyte2','Granulocyte3','B-cell','MK','Granulocyte4',
                           'Stem Cell', 'Monocyte','Macrophage','Erythroid','B-cell pro.','T-cell/NK','MEP'),each = 4)
ggplot(data = tbl_colPerc, aes(x = time2, y = state, fill = `Cluster 0`)) + geom_tile() +
        scale_fill_gradient2(high = 'darkred', low = 'white', mid = 'red', midpoint = 50,
                             limit = c(0,100), space ="Lab",
                             name = 'Percentage of Cells\n') +
        theme_minimal() + ylab('Experimental States') + xlab('') +
        ggtitle ("Heatmap of Cell Distributions within Clusters") + 
        theme (plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_text(aes(time2, state, label = `Cluster 0`), color = 'black', size = 4)+
        coord_fixed()

ggsave('./figures/NAMED_heatmap_cells_within_clusters.png', device = 'png', units = 'in',
       height = 5, width = 8, dpi = 400)


saveRDS(wbm, file = './data/wbm_clustered_filtered.rds')

# Normalizing to get differential cell counts

tbl <- as.data.frame(dcast(setDT(wbm@meta.data), state~seurat_clusters, length))
tbl
rownames(tbl) <- tbl$state
tbl <- tbl [,-1]
colnames(tbl) <- paste0('Cluster ',0:12)
tbl

for (i in 1:dim(tbl)[1]){
        for (j in 1:dim(tbl)[2]){
                tbl[i,j] <- round(tbl[i,j] / rowSums(tbl)[i],3)*1000
        }
}
tbl

# Looking at the UMAP in different ways

wbm <- readRDS('./data/wbm_clustered_filtered.rds')
DimPlot(wbm, reduction = 'umap', label = T, split.by = 'state')
DimPlot(wbm, reduction = 'umap', label = T, repel = T, cols = DiscretePalette(13))
ggsave('./figures/umap_projection.png', device = 'png', units = 'in',
       height = 5, width = 6, dpi = 400)
DimPlot(wbm, reduction = 'umap', split.by = 'state', cols = DiscretePalette(13))
ggsave('./figures/umap_projection_split_by_state.png', device = 'png',
       units = 'in', height = 5, width = 12, dpi = 400)

?DimPlot
head(wbm@meta.data)

x <- DimPlot(wbm, reduction = 'umap')
x + theme_bw()

# Genes for MK
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Itga2b > 1)),
        cols.highlight = 'red') + ggtitle('Itga2b') + theme(plot.title = element_text(hjust = 0.5))
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Pf4 > 2)),
        cols.highlight = 'red') + ggtitle('Pf4') + theme(plot.title = element_text(hjust = 0.5))
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Selp > 1)),
        cols.highlight = 'red') + ggtitle('Selp') + theme(plot.title = element_text(hjust = 0.5))
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cd47 > 3)),
        cols.highlight = 'red') + ggtitle('Cd47') + theme(plot.title = element_text(hjust = 0.5))
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Gata2 > 1)),
        cols.highlight = 'red') + ggtitle('Gata2') + theme(plot.title = element_text(hjust = 0.5))
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Plk3 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Tspan9 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Treml1 > 2)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Runx1 > 3)),
        cols.highlight = 'red')

# Ery
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Klf1 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Tmod1 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Ank1 > 2)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cfp > 2)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Grsf1 > 1)),
        cols.highlight = 'red')

# Ery Prog.
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = `Hba-a1` > 3)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = `Hbb-bt` > 3)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = `Hba-a2` > 3)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Snca > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cd59a > 1)),
        cols.highlight = 'red')

# Mye
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Elane > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Mpo > 2)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Csf3r > 1)),
        cols.highlight = 'red')



# Lym
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Ighd > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Sox4 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Lef1 > 1)),
        cols.highlight = 'red')


# T-cells (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4384382/)
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Gata3 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Tbx21 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Rorc > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Foxp3 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cd5 > 1)),
        cols.highlight = 'red')

# Plasma
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Sik1 > 1)),
        cols.highlight = 'red')

# Proliferating Cells/Stem Cells
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Top2a > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Mki67 > 1)),
        cols.highlight = 'red')

# Macrophages (?monocytes) https://panglaodb.se/markers.html?cell_type=%27Macrophages%27
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cd14 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Ccr5 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cd5 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Slamf9 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Lilra5 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Mgl2 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Ccl12 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Clec4a2 > 2)),
        cols.highlight = 'red')

# Monocytes
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Clec12a > 2)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Cxcl10 > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Psap > 3)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Ifitm3 > 4)),
        cols.highlight = 'red')

# NK
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Il2rb > 1)),
        cols.highlight = 'red')
DimPlot(wbm, reduction = 'umap', cells.highlight = colnames(subset(wbm, subset = Zfp683 > 1)),
        cols.highlight = 'red')

# Looking into the genes that define cluster 9 and 12, I think that 9 is erythroid/mk/mep,
# and I am uncertain about 12 (singleR said neutro/granulo/myeloid)

cluster12.markers <- FindMarkers(wbm, ident.1 = 12, min.pct = .25)
head(cluster12.markers,20)
cluster9.markers <- FindMarkers(wbm, ident.1 = 9, min.pct = .25)
head(cluster9.markers,20)
# alas2: an erythroid-specific mitchondrially located enzyme
# bpgm found at high concentrations in rbc
# Rhag associated with human erythrocyte membrane preoteins


# Trying to pull out the UMAP coordinates
head(wbm@meta.data)
summary(wbm@meta.data)

wbm@meta.data$nFeature_RNA_quart <- ifelse(wbm@meta.data$nFeature_RNA < 1475, 1,
                                           ifelse(wbm@meta.data$nFeature_RNA < 2094, 2,
                                           ifelse(wbm@meta.data$nFeature_RNA < 3168, 3,4)))
wbm@meta.data$nCount_RNA_quart <- ifelse(wbm@meta.data$nCount_RNA < 4946, 1,
                                           ifelse(wbm@meta.data$nCount_RNA < 9847, 2,
                                                  ifelse(wbm@meta.data$nCount_RNA < 17314, 3,4)))
wbm@meta.data$percent.mt_quart <- ifelse(wbm@meta.data$percent.mt < 1.0033, 1,
                                         ifelse(wbm@meta.data$percent.mt < 1.794, 2,
                                                ifelse(wbm@meta.data$percent.mt < 3.6237, 3,4)))

DimPlot(wbm, reduction = 'umap', group.by = 'nFeature_RNA_quart', cols = c('yellow','orange','red','darkred'))
DimPlot(wbm, reduction = 'umap', group.by = 'nCount_RNA_quart', cols = c('yellow','orange','red','darkred'))
DimPlot(wbm, reduction = 'umap', group.by = 'percent.mt_quart', cols = c('yellow','orange','red','darkred'))

# Giving the clusters names!!!

new.cluster.ids <- c('Granulocyte1','Granulocyte2','Granulocyte3','B cell','MK','Granulocyte4','Stem Cell',
                     'Monocyte','Macrophage','Erythroid','B-cell pro','T-cell/NK','MEP')
names(new.cluster.ids) <- levels(wbm)

wbm <- RenameIdents(wbm, new.cluster.ids)

DimPlot(wbm, reduction = 'umap', label = T, pt.size = 0.5, repel = T, split.by = 'state',
        cols = DiscretePalette(13)) + NoLegend()

DimPlot(wbm, reduction = 'umap', label = T, repel = T, cols = DiscretePalette(13))
# ggsave('./figures/NAMED_umap_projection_split_by_state.png', device = 'png',
#        units = 'in', height = 5, width = 12, dpi = 400)
DimPlot(wbm, reduction = 'umap', label = T, repel = T, cols = DiscretePalette(13))

head(wbm@meta.data)
new.cluster.ids
wbm
?RenameIdents

x <- Idents(wbm)
class(x)

wbm@meta.data$cluster_IDs <- Idents(wbm)

# differential expression between wbm and enr_wbm 
wbm <- readRDS('data/wbm_clustered_filtered_named.rds')
library(DESeq2)
mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'WBM-control', ident.2 = 'WBM-mut',
                         verbose = T, group.by = 'state',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'DESeq2')
summary(mk.DE.wbm$p_val_adj < 0.05)
DESeq_genes <- rownames(mk.DE.wbm[mk.DE.wbm$p_val_adj < 0.05,])
DESeq_genes <- DESeq_genes[-grep('NA', DESeq_genes)]

mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'WBM-control', ident.2 = 'WBM-mut',
                         verbose = T, group.by = 'state',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'wilcox')
summary(mk.DE.wbm$p_val_adj < 0.05)
wilcox_genes <- rownames(mk.DE.wbm[mk.DE.wbm$p_val_adj < 0.05,])

mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'WBM-control', ident.2 = 'WBM-mut',
                         verbose = T, group.by = 'state',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'bimod')
summary(mk.DE.wbm$p_val_adj < 0.05)
bimod_genes <- rownames(mk.DE.wbm[mk.DE.wbm$p_val_adj < 0.05,])

mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'WBM-control', ident.2 = 'WBM-mut',
                         verbose = T, group.by = 'state',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'roc')
summary(mk.DE.wbm$power > 0.8)
roc_genes <- rownames(mk.DE.wbm[mk.DE.wbm$power > 0.8,])

ggplot(mk.DE.wbm, aes(x = avg_logFC, y = -log10(p_val))) + geom_point() + 
        geom_hline(yintercept = -log10(0.005), colour = 'red', linetype = 'dashed')
        
write.csv(mk.DE.wbm, './data/DEgenes.MK.wbm.csv', quote = F)

saveRDS(wbm, file = './data/wbm_clustered_filtered_named.rds')


Reduce(intersect, list(DESeq_genes,wilcox_genes,bimod_genes,roc_genes))
Reduce(intersect, list(DESeq_genes,bimod_genes,roc_genes))
Reduce(intersect, list(wilcox_genes,roc_genes))
Reduce(intersect, list(wilcox_genes,bimod_genes))
Reduce(intersect, list(DESeq_genes,wilcox_genes))
Reduce(intersect, list(DESeq_genes,bimod_genes))
Reduce(intersect, list(DESeq_genes,roc_genes))
Reduce(intersect, list(bimod_genes,roc_genes))

all_mk_genes <- unique(c(DESeq_genes, wilcox_genes, bimod_genes, roc_genes))
tb <- matrix(ncol = 4, nrow = length(all_mk_genes), FALSE)
colnames(tb) <- c('DESeq2','Wilcox','Bimod','Roc')
rownames(tb) <- all_mk_genes
for (i in 1:dim(tb)[1]){
        if (rownames(tb)[i] %in% DESeq_genes){
                tb[i,1] = T
        }
        if (rownames(tb)[i] %in% wilcox_genes){
                tb[i,2] = T
        }
        if (rownames(tb)[i] %in% bimod_genes){
                tb[i,3] = T
        }
        if (rownames(tb)[i] %in% roc_genes){
                tb[i,4] = T
        }
        
}

tb[order(tb[,1],tb[,2],tb[,3],tb[,4], decreasing = T),]
#______________________________________________________________________________
md <- wbm@meta.data

head(md)
md <- md[,-10]

library(data.table)
tbl <- as.data.frame(dcast(setDT(wbm@meta.data), state~seurat_clusters, length))
tbl
rownames(tbl) <- tbl$state
tbl <- tbl [,-1]
colnames(tbl) <- paste0('Cluster ',0:12)
colnames(tbl) <- c('Gran1','Gran2','Gran3','B-Cell','MK','Gran4','Stem','Mono','Macro','Eryth','B-Cell pro.','T-Cell','MEP')
tbl


# Getting markers for macrophage and stem cell clusters to help clarify cell labeling

macro <- FindMarkers(wbm, ident.1 = 'Macrophage')
write.csv(macro,'./data/macrophage.cluster.markers.csv', quote = F)

sc <- FindMarkers(wbm, ident.1 = 'Stem Cell')
write.csv(sc,'./data/stemcell.cluster.markers.csv', quote = F)

write.csv(wbm@meta.data,'./data/metadata0330.csv', quote = F)
