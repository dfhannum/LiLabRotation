# Getting marker genes from Qianyi dataset with the combined run

wbm_exp <- readRDS('data/wbm_exp.rds')

head(wbm_exp@meta.data)

library(Seurat)
library(SingleR)
library(data.table)

?SingleR

hpca.se <- hpca
hpca.se$main_types
class(hpca.se)

sessionInfo()

ref <- as.data.frame(fread('./data/reference_qianyi/BlueprintEncode_GeneExp.txt'))
ref[1:10,1:4]
rownames(ref) <- ref$V1
ref <- ref[,-1]
ref[1:10,1:4]
summary(ref[1:5])

?makeSummarizedExperimentFromDataFrame

sref <- makeSummarizedExperimentFromDataFrame(ref)

rownames(ref)
wbm_exp

rownames(data) <- toupper(rownames(data))
rownames(ref)

# recreating the seurat object

common_genes [rownames(data) %in% rownames(ref))

rownames(data)[grepl(pattern = 'Cd4', ignore.case = T,x = rownames(data))]

data <- Read10X('./data/filtered_feature_bc_matrix/')

gene_in_list <- function(gene_name, list = ct){
        return(list[grepl(pattern = gene_name, ignore.case = T, x = list)])
}

gene_in_data('g6b')
# alt
gene_in_data('CLEC1B')
gene_in_data('itga2b')


ct <- colnames(ref)

mk_cells <- ref[,grepl(pattern = 'positive.mega', ignore.case = T, x = colnames(ref))]

head(mk_cells)

not_mk_cells <- ref[,!grepl(pattern = 'positive.mega', ignore.case = T, x = colnames(ref))]

results <- as.data.frame(matrix(ncol = 6, nrow = dim(ref)[1], 0))
colnames(results) <- c('cell_avg','otro_cell_avg','fc','p.value', 'pos_cell_perc', 'neg_cell_perc')
rownames(results) <- rownames(ref)


colnames(ref)


#subsetting the reference down
gene_in_list('mega', ct)

gene_in_list('eryth')

sref <- ref[,grepl('mega',ignore.case = T, x = colnames(ref)) |
                    grepl('lympho',ignore.case = T, x = colnames(ref)) |
                    grepl('erthy',ignore.case = T, x = colnames(ref)) |
                    grepl('myelo',ignore.case = T, x = colnames(ref))]

sref <- sref[rowSums(sref)>0,]

mk_cells <- sref[,grepl('positive.megak', ignore.case = T, x = colnames(sref))]
not_mk_cells <- sref[,!grepl('positive.megak', ignore.case = T, x = colnames(sref))]

write.table(sref, './data/')

for (i in 1:dim(sref)[1]){
        test <- t.test(mk_cells[i,], not_mk_cells[i,])
        results[i,1] <- as.numeric(test$estimate[1])
        results[i,2] <- as.numeric(test$estimate[2])
        results[i,3] <- results[i,1] / results[i,2]
        results[i,4] <- test$p.value
        results[i,5] <- sum(mk_cells[i,] > 0) / ncol(mk_cells)
        results[i,6] <- sum(not_mk_cells[i,] > 0) / ncol(not_mk_cells)
        
}

resul1 <- results[results$p.value < 0.0001 & results$fc >5 & results$neg_cell_perc > .5,]

# now makign it into a function to run the others

get_sig_diff_genes <- function(data = sref, cell_name, p.val.thresh = 0.05){
        pos_type <- data[,grepl(pattern = cell_name, ignore.case = T, x = colnames(data))]
        neg_type <- data[,!grepl(pattern = cell_name, ignore.case = T, x = colnames(data))]
        
        results <- as.data.frame(matrix(ncol = 6, nrow = dim(data)[1], 0))
        colnames(results) <- c('cell_avg','otro_cell_avg','fc','p.value', 'pos_cell_perc', 'neg_cell_perc')
        rownames(results) <- rownames(data)
        
        for (i in 1:dim(data)[1]){
                test <- t.test(mk_cells[i,], not_mk_cells[i,])
                results[i,1] <- as.numeric(test$estimate[1])
                results[i,2] <- as.numeric(test$estimate[2])
                results[i,3] <- results[i,1] / results[i,2]
                results[i,4] <- test$p.value
                results[i,5] <- sum(mk_cells[i,] > 0) / ncol(mk_cells)
                results[i,6] <- sum(not_mk_cells[i,] > 0) / ncol(not_mk_cells)
                
        }
        
        return(results[results$p.value < p.val.thresh,])
}


data = sref
cell_name = 'lymph'

pos_type <- data[,grepl(pattern = cell_name, ignore.case = T, x = colnames(data))]
neg_type <- data[,!grepl(pattern = cell_name, ignore.case = T, x = colnames(data))]

results <- as.data.frame(matrix(ncol = 6, nrow = dim(data)[1], 0))
colnames(results) <- c('cell_avg','neg_cell_avg','fc','p.value', 'pos_cell_perc', 'neg_cell_perc')
rownames(results) <- rownames(data)

for (i in 1:dim(data)[1]){
        test <- t.test(pos_type[i,], neg_type[i,])
        results[i,1] <- as.numeric(test$estimate[1])
        results[i,2] <- as.numeric(test$estimate[2])
        results[i,3] <- results[i,1] / results[i,2]
        results[i,4] <- test$p.value
        results[i,5] <- sum(mk_cells[i,] > 0) / ncol(mk_cells)
        results[i,6] <- sum(not_mk_cells[i,] > 0) / ncol(not_mk_cells)
        
}

# res_myeloid <- results[results$p.value < 0.0001 & results$fc >5 & results$neg_cell_perc > .5,]
results <- results[rowSums(results) != 0,]
head(results[order(results$p.value, decreasing = F),],100)
results <- results[order(results$p.value, decreasing = F),]
results[results$fc > 1,][1:20,]


