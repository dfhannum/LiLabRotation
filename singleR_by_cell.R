# singleR run looking at correlation for individual cell

library(Seurat)
library(SingleR)
library(ggplot2)
library(scRNAseq)

wbm <- readRDS('data/wbm_clustered_filtered_named.rds')
head(wbm@meta.data)
summary(as.factor(paste0(wbm@meta.data$seurat_clusters,wbm@meta.data$cluster_IDs)))

m.ref <- ImmGenData()
m.ref2 <- MouseRNAseqData()

SCwbm <- as.SingleCellExperiment(wbm)

pred <- SingleR(test = SCwbm,
                ref = list(igd = m.ref, mrd = m.ref2),
                labels = list(m.ref$label.main, m.ref2$label.main),
                method = 'single')
table(pred$labels)
pred[1:10,1:4]

#BiocManager::install('pheatmap')
library(pheatmap)

plotScoreHeatmap(pred)

# Looks interesting and now time to compare!

summary(rownames(pred) == rownames(wbm@meta.data))

wbm@meta.data$cell_id <- pred$pruned.labels

table(wbm@meta.data$cluster_IDs)
table(wbm@meta.data$cell_id)
tbl <- table(wbm@meta.data$cluster_IDs,wbm@meta.data$cell_id)
class(tbl)
head(as.data.frame(tbl))

clst_vs_cell_id <- as.data.frame(tbl)

colnames(clst_vs_cell_id) <- c('Cluster ID', 'Cell ID', 'Count')

cluster_counts <- c(summary(as.factor(wbm@meta.data$cluster_IDs)))
clst_vs_cell_id$cluster_counts <- rep(cluster_counts, 14)
clst_vs_cell_id$cluster_perc <- round(clst_vs_cell_id$Count / clst_vs_cell_id$cluster_counts,2) *100

ggplot(data = clst_vs_cell_id, aes(x = `Cell ID`, y = `Cluster ID`, fill = cluster_perc)) + geom_tile() +
        scale_fill_gradient2(high = 'darkred', low = 'white', mid = 'red', midpoint = 50,
                             limit = c(0,100), space ="Lab",
                             name = 'Percentage of Cells\n') +
        ylab('Cluster Labels') + xlab('Cell Labels') +
        #theme_minimal() +
        ggtitle ("Heatmap of Individual Cell Labels within Cluster Labels") + 
        theme (plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
        #geom_text(aes(time2, state, label = `Cluster 0`), color = 'black', size = 4)+
        coord_fixed()

ggsave('./figures/clusterIDS_vs_cellIDS_heatmap.png', device = 'png', units = 'in',
       height = 5, width = 8, dpi = 400)
        