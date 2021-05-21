library(SingleCellExperiment)
library(SC3)
library(scater)

filenames<-'zeisel'
data <- readRDS(paste0(filenames,'.rds'))@assays$data$counts

sce <- SingleCellExperiment(
    assays = list(
    counts = as.matrix(data),
    logcounts = log2(as.matrix(data) + 1)))

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

sce <- sc3(sce, ks = 47, biology = FALSE, n_cores = 10)
label<-as.integer(sce@colData@listData[[1]])
write.csv(label, file=paste0(filenames,'_data_SC3.csv'))