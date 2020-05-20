# adding a script to put together potential figures for the paper

library(Seurat)
library(ggplot2)

wbm <- readRDS('./data/wbm_clustered_filtered_named.rds')

head(wbm@meta.data)

DimPlot(wbm, reduction = 'umap', cols= DiscretePalette(13))
ggsave('./figures/exp1_umap.tiff', device = 'tiff', height = 5, width = 8, units = 'in', dpi =300)
?ggsave

DimPlot(wbm, reduction = 'umap', cols = DiscretePalette(13), split.by = 'state') + NoLegend()
ggsave('./figures/exp1_umap_split_by_state.tiff', device = 'tiff', 
       height = 5, width = 14, units = 'in', dpi = 300)

# these aren't translating properly to adobe, the sizes are all out of wack

# DE Genes for MKs

?FindMarkers

mk_markers <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                          verbose = T, group.by = 'condition',
                          subset.ident = 'MK', min.diff.pct = .2,
                          logfc.threshold = log(2),
                          test.use = 'wilcox')
dim(mk_markers[mk_markers$p_val_adj < 0.01,])

DE_genes <- list()
idents <- levels(wbm@active.ident)
for (i in idents){
        print(i)
        
        x <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                         verbose = T, group.by = 'condition',
                         subset.ident = i, min.diff.pct = .2,
                         logfc.threshold = log(2),
                         test.use = 'wilcox')
        DE_genes[[i]] <- x
        
}
DE_genes

DE_genes$MK <- DE_genes$MK[DE_genes$MK$p_val_adj < 0.01,]
mk_markers <- DE_genes$MK

up_mk_mrks <- mk_markers[mk_markers$avg_logFC > 0,]
down_mk_mrks <- mk_markers[mk_markers$avg_logFC < 0,]
dim(mk_markers)
dim(up_mk_mrks)
dim(down_mk_mrks)

up_genes_mk <- rownames(up_mk_mrks)
down_genes_mk <- rownames(down_mk_mrks)

write.table(up_genes_mk,'./data/mk_up_genes.txt', quote = F, row.names = F)
write.table(down_genes_mk, './data/mk_down_genes.txt', quote = F, row.names = F)
write.table(mk_markers, './data/mk_DE_table.txt', quote = F)
dim(wbm)

# down regulated genes go term analysis

go_down_mk <- read.table('./data/go_term_MK_down_genes_analysis', header = T, sep = '\t')
head(go_down_mk)

colnames(go_down_mk)[2:3] <- c('Biological Process', '-log10 FDR')
go_down_mk$`-log10 FDR` <- -log10(go_down_mk$`-log10 FDR`)
go_down_mk <- go_down_mk[order(go_down_mk$`-log10 FDR`, decreasing = F),]
go_down_mk_order <- go_down_mk$`Biological Process`

error_line <- -log10(0.05)

ggplot(data = go_down_mk, aes(x = `Biological Process`, y = `-log10 FDR`)) +
        geom_bar(stat = 'identity', fill = 'blue') + coord_flip() +
        scale_x_discrete(limits = go_down_mk_order) +
        #ggtitle('Down-Regulated Genes') + 
        ylab ('-log10 Corrected P-Value') + xlab('Biological Process') +
        geom_hline(yintercept = error_line, linetype = 2,
                   color = 'red', show.legend = T, size = 1.1) +
        theme_bw() +
        theme(text=element_text(size = 10, family = 'sans'),
              axis.text = element_text(size = 6.5))

ggsave('./figures//go_term_plot_MK_down_reg.tiff', device = 'tiff', height = 3,
       width = 3.75, units = 'in')

# heatmap or dot plot

features <- c('Fcer1a','Itga2b','Ank1','Csf3r','Sox4','Gata3','Top2a','Mki67', 'Fcgr1', 'Dysf')

DotPlot(wbm, features = features) + RotatedAxis()

wbm.markers <- FindAllMarkers(wbm, only.pos = T, min.pct = 0.25, logfc.threshold = log(2))
library(dplyr)
top5 <- wbm.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(wbm, features = top5$gene) + NoLegend()
ggsave('./figures/top5_heatmap.tiff', device = 'tiff', height = 4,
       width = 8, units = 'in')

