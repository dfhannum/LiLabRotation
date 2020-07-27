# Running experiment one in Monocle

# http://cole-trapnell-lab.github.io/monocle-release/docs/#getting-started-with-monocle

library(monocle)

mtx <- readMM('./data/filtered_feature_bc_matrix/matrix.mtx')

mtx[1:10,1:10]

dim(mtx)

features <- read.table('./data/filtered_feature_bc_matrix/features.tsv')
colnames(features)[2] <- 'gene_short_name'
head(features)
dim(features)

barcodes <- read.table('./data/filtered_feature_bc_matrix/barcodes.tsv')
head(barcodes)

pd <- new("AnnotatedDataFrame", data = barcodes)
fd <- new("AnnotatedDataFrame", data = features)

library(Matrix)
cds <- newCellDataSet(as(mtx, "sparseMatrix"),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersionspDa(cds)
cds <- detectGenes()