library(SIMLR)
library(SingleCellExperiment)

file.list<-dir('/home/qmliu/projects/DMCluster/rdata')
ks.list<-c(3,7,3,3,39,3,11,49,5,5,11,47)

# large.data<-c('Data_Macosko','Data_Zeisel','Data_Tasic')
large.data<-c('Data_Zeisel','Data_Tasic')
# for (i in c(1:(length(file.list)-1))) {
ncell = 200
for (i in c(5:5)) {
  data.names<-unlist(strsplit(x = file.list[i], split = '.', fixed = TRUE))[1]
  data <- readRDS(paste0('/home/qmliu/projects/DMCluster/rdata/',file.list[i]))
  true.label<- data@colData@listData[["clust_id"]][1:ncell]
  if (is.null(true.label)){true.label<- data@colData@listData[["V1"]]}
  data.counts<-data@assays@data@listData[["counts"]][1:ncell,]
  # data.counts<-data@assays[["data"]]@listData[["counts"]]
  
  if (isFALSE(data.names %in% large.data)){
    t2<-proc.time()
    res<-SIMLR(X = data.counts, c = ks.list[i], cores.ratio = 0)
  }else{
    t2<-proc.time()
    res<-SIMLR_Large_Scale(X = data.counts, c= ks.list[i], k = 10, kk = 100, if.impute = FALSE, normalize = FALSE)
  }
  print(paste('the whole costs time:', (proc.time() -t2)[3]))
  label<-res[["y"]][["cluster"]]
  
  if (isFALSE(file.exists('/home/qmliu/projects/DMCluster/rdata/res'))) {
    dir.create('./home/qmliu/projects/DMCluster/rdata/res')
  }
  
  label.df<-data.frame('prediction'=label, 'label'=true.label)
  write.csv(label.df, file=file.path('/home/qmliu/projects/DMCluster/rdata/res',paste0(data.names,'_SIMLR.csv')))
}