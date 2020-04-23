# Running anlaysis with Unassigned cells

library(Seurat)
library(ggplot2)

wbm_exp <- readRDS('data/wbm_exp.rds')

head(wbm_exp@meta.data)

VlnPlot(wbm_exp, features = c('nFeature_RNA','nCount_RNA','percent.mt'), ncol = 3)

ggplot(wbm_exp@meta.data, aes(x = percent.mt, colour = celltype)) + geom_density() +
        geom_vline(xintercept = 10, colour = 'black', linetype = 'dashed')

ggplot(wbm_exp@meta.data, aes( x = nFeature_RNA, colour = celltype)) + 
        geom_vline(xintercept = c(500,5000), colour = 'red', linetype = 'dashed') +
        geom_density()

ggplot(wbm_exp@meta.data, aes( x= nCount_RNA, colour = celltype)) + 
        geom_density() 
        

wbm_exp <- subset(wbm_exp, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

summary(as.factor(wbm_exp@meta.data$celltype))

wbm_exp <- NormalizeData(wbm_exp, normalization.method = 'LogNormalize', scale.factor = 10000)

wbm_exp <- FindVariableFeatures(wbm_exp, selection.method = 'vst', nfeatures = 2000)

top10 <- head(VariableFeatures(wbm_exp), 10)

plot1 <- VariableFeaturePlot(wbm_exp)
LabelPoints(plot = plot1, points = top10, repel = T)

all.genes <- rownames(wbm_exp)
wbm_exp <- ScaleData(wbm_exp, features = all.genes)

wbm_exp <- RunPCA(wbm_exp, features = VariableFeatures(object = wbm_exp))

DimPlot(wbm_exp, reduction = 'pca')
DimPlot(wbm_exp, reduction = 'pca', split.by = 'celltype')
# Not using
wbm_exp <- JackStraw(wbm_exp, num.replicate = 100)

wbm_exp <- ScoreJackStraw(wbm_exp, dims = 1:20)

JackStrawPlot(wbm_exp, dims = 1:20)

ElbowPlot(wbm_exp)

#using 15

wbm_exp <- FindNeighbors(wbm_exp, dims = 1:15)


# x <- FindClusters(wbm_exp, resolution = 0.4)
# length(levels(x$seurat_clusters))

lst_of_resolutions <- seq(0,1, by = .05)

resul = c()
cnt = 1
for (i in lst_of_resolutions){
        #print(i)
        x <- FindClusters(wbm_exp, resolution = i)
        resul[cnt] <- length(levels(x$seurat_clusters))
        cnt = cnt + 1
}

plot(resul)
lst_of_resolutions[10]

wbm_exp <- FindClusters(wbm_exp, resolution = 0.45)

wbm_exp <- RunUMAP(wbm_exp, dims = 1:15)
DimPlot(wbm_exp, reduction = 'umap', label = T, split.by = 'celltype')

# Getting the counts of cell type per cluster

md <- wbm_exp@meta.data
head(md)
md <- md[,c(1:8, 10)]
md$state <- as.factor(paste0(md$celltype,'-',md$condition))

library(dplyr)
tbl <- as.data.frame(summarise(group_by(md, state, seurat_clusters), count = n()))
summary(md$state)

tbl <- reshape(tbl, idvar = 'state', timevar = 'seurat_clusters', direction = 'wide')

rownames(tbl) <- tbl$state
tbl <- tbl[,-1]
colnames(tbl) <- paste0('Cluster ', 0:14)
tbl[2,12] <- 0
tbl

# Percentage of each state in a cluster

tbl_rowPerc <- tbl

for (i in 1:dim(tbl)[1]){
        for (j in 1:dim(tbl)[2]){
                tbl_rowPerc[i,j] <- round(tbl[i,j] / rowSums(tbl)[i],2) * 100
        }
}

tbl1 <- as.matrix(tbl_rowPerc)

heatmap(tbl1)
cluster_order <- paste0('Cluster ',c(3,4,2,1,0,5,9,7,11,10,14,13,12,6,8))

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

ggplot(data = tbl_rowPerc, aes(x = time, y = state, fill = `Cluster 0`)) + geom_tile() +
        scale_fill_gradient2(high = 'darkred', low = 'white', mid = 'red', midpoint = 50,
                             limit = c(0,100), space ="Lab",
                             name = 'Percentage of Cells\nper Cluster ') +
        theme_minimal() + ylab('Experimental States') + xlab('') +
        ggtitle ("Heatmap of Cell Distributions within Experimental States") + 
        theme (plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_text(aes(time, state, label = `Cluster 0`), color = 'black', size = 4)+
        coord_fixed()


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

ggplot(data = tbl_colPerc, aes(x = time, y = state, fill = `Cluster 0`)) + geom_tile() +
        scale_fill_gradient2(high = 'darkred', low = 'white', mid = 'red', midpoint = 50,
                             limit = c(0,100), space ="Lab",
                             name = 'Percentage of Cells\n') +
        theme_minimal() + ylab('Experimental States') + xlab('') +
        ggtitle ("Heatmap of Cell Distributions within Clusters") + 
        theme (plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_text(aes(time, state, label = `Cluster 0`), color = 'black', size = 4)+
        coord_fixed()
