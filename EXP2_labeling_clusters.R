exp2 <- readRDS('./data/EXP2_clustered_filtered_SingleR_labels.rds')


# MK Markers
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Itga2 >.5)),
        cols.highlight = 'red') + ggtitle('Itga2b') + theme(plot.title = element_text(hjust = 0.5)) +
        NoLegend()
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Pf4 > 3)),
        cols.highlight = 'red') + ggtitle('Pf4') + theme(plot.title = element_text(hjust = 0.5)) +
        NoLegend()
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd47 > 2)),
        cols.highlight = 'red') + ggtitle('Cd47') + theme(plot.title = element_text(hjust = 0.5)) + 
        NoLegend()
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd47 > 3)),
        cols.highlight = 'red') + ggtitle('Cd47') + theme(plot.title = element_text(hjust = 0.5)) +
        NoLegend()
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd47 > 2)),
        cols.highlight = 'red') + ggtitle('Cd47') + theme(plot.title = element_text(hjust = 0.5)) +
        NoLegend()

DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Gata2 > 1)),
        cols.highlight = 'red') + ggtitle('Gata2') + theme(plot.title = element_text(hjust = 0.5)) +
        NoLegend()

# Ery
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Klf1 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Tmod1 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Ank1 > 2)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cfp > 2)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Grsf1 > 1)),
        cols.highlight = 'red')

# Ery Prog.
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = `Hba-a1` > 3)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = `Hbb-bt` > 3)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = `Hba-a2` > 3)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Snca > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd59a > 1)),
        cols.highlight = 'red')

# Mye
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Elane > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Mpo > 2)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Csf3r > 1)),
        cols.highlight = 'red')



# Lym
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Ighd > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Sox4 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Lef1 > 1)),
        cols.highlight = 'red')


# T-cells (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4384382/)
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Gata3 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Tbx21 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Rorc > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Foxp3 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd5 > 1)),
        cols.highlight = 'red')

# Plasma
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Sik1 > 1)),
        cols.highlight = 'red')

# Proliferating Cells/Stem Cells
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Top2a > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Mki67 > 1)),
        cols.highlight = 'red')

# Macrophages (?monocytes) https://panglaodb.se/markers.html?cell_type=%27Macrophages%27
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd14 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Ccr5 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cd5 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Slamf9 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Lilra5 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Mgl2 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Ccl12 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Clec4a2 > 2)),
        cols.highlight = 'red')

# Monocytes
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Clec12a > 2)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Cxcl10 > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Psap > 3)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Ifitm3 > 4)),
        cols.highlight = 'red')

# NK
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Il2rb > 1)),
        cols.highlight = 'red')
DimPlot(exp2, reduction = 'umap', cells.highlight = colnames(subset(exp2, subset = Zfp683 > 1)),
        cols.highlight = 'red')

