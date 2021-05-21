library(Seurat)
library(dplyr)
library(Matrix)

filenames<-'marques'
data <- readRDS(paste0(filenames,'.rds'))@assays$data$counts
# data <- read.csv("data.csv", row.names = 1)
pbmc <- CreateSeuratObject(counts = data)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
write.csv(Idents(pbmc), file=paste0(filenames,'_data_Seurat.csv'))