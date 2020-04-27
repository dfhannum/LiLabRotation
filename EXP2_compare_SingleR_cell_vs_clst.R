# singleR for experiment 2

library(SingleR)
library(ggplot2)
library(scRNAseq)

exp2 <- readRDS('data/EXP2_clustered_filtered.rds')

head(exp2@meta.data)

m.ref <- ImmGenData()
m.ref2 <- MouseRNAseqData()

m.ref

SCexp2 <- as.SingleCellExperiment(exp2)

# For the clusters 

pred_cluster <- SingleR(test = SCexp2, ref = list(m.ref, m.ref2),
                        labels = list(m.ref$label.main, m.ref2$label.main),
                        method = 'cluster',
                        clusters = SCexp2$seurat_clusters)
pred_cell <- SingleR(test = SCexp2, ref = list(m.ref, m.ref2),
                     labels = list(m.ref$label.main, m.ref2$label.main),
                     method = 'single')

md <- exp2@meta.data

summary(rownames(md)== rownames(pred_cell))

table(pred_cluster$labels, pred_cluster$pruned.labels)

head(pred_cell[,1:4])

head(md)

md <- cbind(md, pred_cell[,4])
colnames(md)[7] <- 'Cell ID'


# testing 
library(plyr)

mapvalues(pred_cluster$cluster, from = pred_cluster$cluster, to = pred_cluster$pruned.labels)

x <- as.factor(c(0:15,1,12,13,14))
x
mapvalues(x, from = levels(x), to = pred_cluster$pruned.labels )

# this works
pr.lab <- pred_cluster$pruned.labels
summary(as.factor(pr.lab))
pr.lab[c(1:5,11)] <- paste0('Granuloyctes', 1:6)
pr.lab[c(8,9,12,13)] <- paste0('B cells', 1:4)
pr.lab[7] <- 'unknown'

md$`Cluster ID` <- mapvalues(md$seurat_clusters, from = levels(md$seurat_clusters), to = pr.lab)

# plotting

tbl <- table(md$`Cluster ID`, md$`Cell ID`)
tbl
cl_vs_cel <- as.data.frame(tbl)

head(cl_vs_cel,20)

clst_counts <- c(summary(as.factor(md$`Cluster ID`)))
cl_vs_cel$cluster_counts <- rep(clst_counts, 15)
cl_vs_cel$cluster_perc <- round(cl_vs_cel$Freq / cl_vs_cel$cluster_counts,2)*100

head(cl_vs_cel)
colnames(cl_vs_cel)[1:3] <- c('Cluster ID','Cell ID', 'Count')

ggplot(data = cl_vs_cel, aes(x = `Cell ID`, y = `Cluster ID`, fill = cluster_perc)) + geom_tile() +
        scale_fill_gradient2(high = 'darkred', low = 'white', mid = 'red', midpoint = 50,
                             limit = c(0,100), space ="Lab",
                             name = 'Percentage of Cells\n') +
        ylab('Cluster Labels') + xlab('Cell Labels') +
        ggtitle ("Heatmap of Individual Cell Labels within Cluster Labels") + 
        theme (plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
        coord_fixed()

ggsave('./figures/EXP2clusterIDS_vs_cellIDS_heatmap.png', device = 'png', units = 'in',
       height = 5, width = 8, dpi = 400)


# adding the cluster ids to the seurat object

head(exp2@meta.data)

exp2$cluster <- mapvalues(exp2$seurat_clusters, from = levels(exp2$seurat_clusters), to = pr.lab)

summary(as.factor(paste0(exp2$seurat_clusters, exp2$cluster)))



names(pr.lab) <- levels(exp2)
pr.lab
exp2 <- RenameIdents(exp2, pr.lab)

DimPlot(exp2, reduction = 'umap', cols = DiscretePalette(16))
DimPlot(exp2, reduction = 'umap', label = TRUE, repel = T,cols = DiscretePalette(16)) + NoLegend()
ggsave('./figures/EXP2_umap_projection.png', device = 'png', units = 'in',
       height = 5, width = 12, dpi = 400)


saveRDS(exp2,'data/EXP2_clustered_filtered_SingleR_labels.rds')

