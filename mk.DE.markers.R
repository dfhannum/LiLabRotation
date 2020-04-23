# MK Differential Expression

library(Seurat)
library(DESeq2)

wbm <- readRDS('data/wbm_clustered_filtered_named.rds')
library(DESeq2)

head(wbm@meta.data)
levels(Idents(wbm))
t
for (t in tests){
        print(t)
        if (t %in% tests[1:3]){
                print ('not roc')
        }
}
tests[1:3]

DE_list <- list()
for (i in levels(Idents(wbm))){
        DE_list[[i]] <- list()
        
        for (t in tests){
                DE_list[[i]][[t]] <- c(1,2,3,33,4,5)
        }
}

for (i in levels(Idents(wbm))){
        mtx <- as.data.frame(matrix(ncol = length(tests)))
        colnames(mtx) <- tests
        
        DE_list[[i]] <- list()
        for (t in tests){
                DE_list[[i]]
        }
}

DE_list <- list()
for (i in levels(Idents(wbm))){
        
        tests <- c('DESeq2','wilcox','bimod','roc')
        
        DE_list[[i]] <- list()
        
        for (t in tests){
                print(paste0(i,'-',t))
                mrk <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                                         verbose = T, group.by = 'condition',
                                         subset.ident = i, logfc.threshold = log(2),
                                         test.use = t)
                
                if (t %in% tests[1:3]){
                        mrks <- mrk[mrk$p_val_adj < 0.05,]
                }
                else{
                        mrks <- mrk[mrk$power > 0.80,]
                }
                DE_list[[i]][[t]] <- mrks
        }
}

saveRDS(DE_list,'./data/DEgenes_all_list')

# going through the list

DE_list$Granulocyte1$DESeq2[c('Rpsa','Thbs1','BC100530'),]

g1 <- c('Rpsa','Thbs1','BC100530')

DE_list$Granulocyte2$DESeq2[c('Wfdc17','Lrg1','Steap4','Ifitm1','Apoe'),]


DE_list$Granulocyte4$DESeq2[c('H1f0','S100a6','Stfa1','BC100530'),]

DE_list$`B cell`$DESeq2[c(rownames(DE_list$`B cell`$wilcox),rownames(DE_list$`B cell`$bimod)),]

DE_list$`B-cell pro`

dg <- rownames(DE_list$MK$DESeq2)
wg <- rownames(DE_list$MK$wilcox)
bg <- rownames(DE_list$MK$bimod)
rg <- rownames(DE_list$MK$roc)
dw <- intersect(dg,wg)
db <- intersect(dg,bg)
dr <- intersect(dg,rg)

# rg genes are in the other three lists, so those 20 are in all 4!

de_MK <- DE_list$MK$wilcox[wg %in% rg,]

write.csv(rg, './data/DEgenes.MK.wbm.csv', quote = F, row.names = F)
# no go terms popped up

DE_list$`Stem Cell`$DESeq2[rownames(DE_list$`Stem Cell`$wilcox),]

dg <- rownames(DE_list$Monocyte$DESeq2)
wg <- rownames(DE_list$Monocyte$wilcox)
bg <- rownames(DE_list$Monocyte$bimod)
de_mono <- DE_list$Monocyte$wilcox[intersect(dg,intersect(wg,bg)),]
write.csv(rownames(de_mono), './data/DEgenes.Mono.wbm.csv', quote = F, row.names = F)
# a bunch of 'response' processes

wg <- rownames(DE_list$Macrophage$wilcox)
dg <- rownames(DE_list$Macrophage$DESeq2)
bg <- rownames(DE_list$Macrophage$bimod)
intersect(wg,intersect(dg,bg))

DE_list$Erythroid
wg <- rownames(DE_list$Erythroid$wilcox)
dg <- rownames(DE_list$Erythroid$DESeq2)
bg <- rownames(DE_list$Erythroid$bimod)
DE_list$Erythroid$wilcox[intersect(wg,bg),]

BWDE_list$MEP



dim(wbm)

mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                         verbose = T, group.by = 'condition',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'wilcox')


wilcox.mk.markers <- mk.DE.wbm[mk.DE.wbm$p_val_adj < 0.05,]

mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                         verbose = T, group.by = 'condition',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'bimod')
bimod.mk.markers <- mk.DE.wbm[mk.DE.wbm$p_val_adj < 0.05,]

mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                         verbose = T, group.by = 'condition',
                         subset.ident = 'MK', logfc.threshold = log(2),
                         test.use = 'roc')
roc.mk.markers <- mk.DE.wbm[mk.DE.wbm$power > 0.80,]

DESeq_genes <- rownames(DEseq.mk.markers)
wilcox_genes <- rownames(wilcox.mk.markers)
bimod_genes <- rownames(bimod.mk.markers)
roc_genes <- rownames(roc.mk.markers)

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

tb <- tb[order(tb[,1],tb[,2],tb[,3],tb[,4], decreasing = T),]
all.markers <- tb[1:20,]
all.markers
marker.genes <- rownames(all.markers)

DEseq.mk.markers[rownames(DEseq.mk.markers) %in% marker.genes,][1:5,]
wilcox.mk.markers[rownames(wilcox.mk.markers) %in% marker.genes,]
bimod.mk.markers[rownames(bimod.mk.markers) %in% 'Slpi',]
roc.mk.markers[rownames(roc.mk.markers) %in% 'Slpi',]

# going to send over the wilcox

write.csv(wilcox.mk.markers[rownames(wilcox.mk.markers) %in% marker.genes,],'./data/wilcox_DE_inMK_control_vs_mutant.csv', quote = F, row.names = T)
x <- read.csv('./data/wilcox_DE_inMK_control_vs_mutant.csv')

rownames(wilcox.mk.markers[rownames(wilcox.mk.markers) %in% marker.genes,])
write.csv(rownames(wilcox.mk.markers[rownames(wilcox.mk.markers) %in% marker.genes,]), 
          './data/wilcox_DE_inMK_control_vs_mutant_GENES.csv', row.names = F,
          quote = F)


# not sure if I need this chunk
mk.DE.wbm <- FindMarkers(wbm, ident.1 = 'control', ident.2 = 'mut',
                         verbose = T, group.by = 'condition',
                         subset.ident = i, logfc.threshold = log(2),
                         test.use = t)
print(i,t)
DEseq.mk.markers <- mk.DE.wbm[mk.DE.wbm$p_val_adj < 0.05,]
head(DEseq.mk.markers)
