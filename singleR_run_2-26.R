# Using SingleR full data (including Unassigned)

library(Seurat)
library(SingleR)
library(ggplot2)
library(scRNAseq)
library(data.table)

wbm <- readRDS('data/wbm_clustered_filtered.rds')
head(wbm@meta.data)

?SingleR

m.ref <- ImmGenData()
m.ref2 <- MouseRNAseqData()
m.ref$label.main[grepl('mega', ignore.case = T, m.ref$label.main)]

# the reference with the SingleR package doesn't have meagkaryocytes but
# still going to run to see what results I get
# I do have reference data but would need to figure out how to format it
m.ref

# example to figure out the data type
test <- LaMannoBrainData('human-es')


SCwbm <- as.SingleCellExperiment(wbm)

pred <- SingleR(test = SCwbm, ref = m.ref, labels = m.ref$label.main,
                method = 'cluster',
                clusters = SCwbm$seurat_clusters)
pred$labels
table(pred$labels)

lab <- tstrsplit(pred$pruned.labels, ' \\(',keep = 1)
lab <- lab[[1]]
summary(as.factor(lab))

wbm[['pred_labl']] <- lab


tbl <- as.data.frame(dcast(setDT(wbm@meta.data), seurat_clusters~pred_labl, length))
tbl
rownames(tbl) <- tbl$seurat_clusters
tbl <- tbl [,-1]
tbl

heatmap(as.matrix(tbl))

## poor results, now doing it with both data references

pred2 <- SingleR(test = SCwbm,
                 ref = list(igd = m.ref, mrd = m.ref2),
                 labels = list(m.ref$label.main, m.ref2$label.main))
table(pred2$labels)
wbm[['pred_labl']] <- pred2$labels
tbl <- as.data.frame(dcast(setDT(wbm@meta.data), seurat_clusters~pred_labl, length))
tbl
rownames(tbl) <- tbl$seurat_clusters
tbl <- tbl [,-1]
tbl

heatmap(as.matrix(tbl))

#using the cluster method
pred3 <- SingleR(test = SCwbm,
                 ref = m.ref2,
                 labels = m.ref2$label.main,
                 method = 'cluster',
                 clusters = SCwbm$seurat_clusters)
pred2c <- SingleR(test = SCwbm,
                 ref = list(igd = m.ref, mrd = m.ref2),
                 labels = list(m.ref$label.main, m.ref2$label.main),
                 method = 'cluster',
                 clusters = SCwbm$seurat_clusters)
table(pred2c$labels)
pred2c[1:13,1:4]
# got results I am not completely in love with, also missing key lables

## Getting a dataframe 

library(data.table)
ref <- read.table('./data/reference_qianyi/BlueprintEncode_GeneExp.txt', row.names = 1)
ref[1:10,1:10]
dim(ref)

head(ref[,grepl('B.cell', ignore.case = T, x = colnames(ref))])
head(ref[,grepl('T.cell', ignore.case = T, x = colnames(ref))])
head(ref[,grepl('kill', ignore.case = T, x = colnames(ref))])
head(ref[,grepl('neutro', ignore.case = T, x = colnames(ref))])



sref <- ref[,grepl('mega',ignore.case = T, x = colnames(ref)) |
                    grepl('lympho',ignore.case = T, x = colnames(ref)) |
                    grepl('erthy',ignore.case = T, x = colnames(ref)) |
                    grepl('myelo',ignore.case = T, x = colnames(ref)) |
                    grepl('B.cell',ignore.case = T, x = colnames(ref)) |
                    grepl('T.cell',ignore.case = T, x = colnames(ref)) |
                    grepl('kill',ignore.case = T, x = colnames(ref)) |
                    grepl('plasma',ignore.case = T, x = colnames(ref)) |
                    grepl('macro',ignore.case = T, x = colnames(ref)) |
                    grepl('stem',ignore.case = T, x = colnames(ref)) |
                    grepl('mono',ignore.case = T, x = colnames(ref)) |
                    grepl('neutro',ignore.case = T, x = colnames(ref))]
dim(sref)

colnames(ref)

celltypes <- read.table('./data/reference_qianyi/BlueprintEncode_celltype.txt', sep = '\t')
colnames(celltypes) <- c('label.fine','label.main')
rownames(celltypes) <- rownames
library(scater)

rownames(ref) <- paste0(str_sub(rownames(ref),1,1),tolower(str_sub(rownames(ref),2))) 
compl_ref <- ref[complete.cases(ref),]

compl_ref <- compl_ref[rownames(compl_ref) %in% rownames(SCwbm),]
ref2 <- SummarizedExperiment(list(logcounts = compl_ref),
                             colData= celltypes)

library(stringr)

pred2c <- SingleR(test = SCwbm,
                  ref = m.ref,
                  labels = m.ref$label.main,
                  method = 'cluster',
                  clusters = SCwbm$seurat_clusters)

pred2 <- SingleR(test = SCwbm,
                 ref = ref2, 
                 labels = ref2$label.main,
                 method = 'cluster',
                 clusters = SCwbm$seurat_clusters)

# Need to fix all the freaking row names because they don't match up

sc <- rownames(SCwbm)
rn <- rownames(ref)

summary(sc %in% rn)
t1_rn <- rn
class(rn[1])

rn[1]
library(stringr)

summary(sc %in% paste0(str_sub(rn,1,1),tolower(str_sub(rn,2))))

rownames(ref) <- paste0(str_sub(rn,1,1),tolower(str_sub(rn,2))) 

wbm_count_data <- GetAssayData(object = wbm[['RNA']])
dim(wbm_count_data)

as.matrix(wbm_count_data)[1:10,1:10]
wbm_count_data <- as.matrix(wbm_count_data)
VariableFeatures(wbm)

#++++++++++++++++++++++++++++++++++++++++++++++++
# Go #2
library(stringr)
library(Seurat)
library(SingleR)
library(ggplot2)

wbm <- readRDS('data/wbm_clustered_filtered.rds')
head(wbm@meta.data)
ref <- read.table('./data/reference_qianyi/BlueprintEncode_GeneExp.txt', row.names = 1)
celltypes <- read.table('./data/reference_qianyi/BlueprintEncode_celltype.txt', sep = '\t')

rn <- rownames(ref)
summary(VariableFeatures(wbm) %in% paste0(str_sub(rn,1,1),tolower(str_sub(rn,2))) )
VariableFeatures(wbm)[!(VariableFeatures(wbm) %in% paste0(str_sub(rn,1,1),tolower(str_sub(rn,2))))]

rn[grepl('-', ignore.case = T, rn)]
rownames(ref) <- paste0(str_sub(rn,1,1),tolower(str_sub(rn,2)))
         
ref2 <- ref[rownames(ref) %in% VariableFeatures(wbm),]
wbm2 <- GetAssayData(object = wbm[['RNA']])
wbm2 <- wbm2[rownames(wbm2) %in% VariableFeatures(wbm),]
wbm2 <- wbm2[rownames(wbm2) %in% rownames(ref2),]
wbm2 <- as.matrix(wbm2)

pred2 <- SingleR(test = wbm2, ref = ref, labels = colnames(ref))

# following the example in ?SingleR
which(is.na(ref2), arr.ind = T)
ref2 <- ref2[!rownames(ref2) == 'Ero1l',]

ref <- SummarizedExperiment(
        list(logcounts = ref2),
        colData = DataFrame(label = celltypes)
)
rownames(ref) <- rownames(ref2)

#need this step and only get errors
ref <- scater::logNormCounts(rf)
trained <- trainSingleR(ref, ref$label.V2)

which(is.na(ref2), arr.ind = T)
ref2 <- ref2[!rownames(ref2) == 'Ero1l',]
# getting clusters

wbm2 <- wbm2[!rownames(wbm2) == 'Ero1l',]
cn <- colnames(wbm2)

test <- SummarizedExperiment(
        list(counts = wbm2),
        colData = DataFrame(cluster = wbm@meta.data$seurat_clusters)
)

rownames(test) <- rownames(wbm2)

test <- scater::logNormCounts(test)
pred2 <- SingleR(test,ref, labels = ref$label.V2, method = 'cluster',
                 clusters = test$cluster)
