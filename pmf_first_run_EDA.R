# First look and EDA on PMF data

library(dplyr)
library(Seurat)

pmf.data <- Read10X(data.dir = "data/filtered_feature_bc_matrix/")

pmf <- CreateSeuratObject(counts = pmf.data, project = 'pmf11k', min.cells = 3, min.features = 200)
pmf

pmf.data
