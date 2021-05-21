
# Note ----------------------------------------------------------------------
# we simulate single-cell RNA sequencing count data with a varying number of cells and cell groups, 
# with different degree of cluster separability and varying rate of dropouts.

library(splatter)


# simulation step 1 -------------------------------------------------------
# a varying number of cells and cell groups
params <- newSplatParams()
nCells.list<-c(500,1000,3000)
nGroups.list<-list(c(0.25,0.25,0.25,0.25),
                   c(rep(0.1,10)),
                   c(rep(0.05,20)))
sim <-list()
params.list <- list(mean.rate = 0.5,
                    seed = 644049,
                    dropout.type = 'experiment',
                    nGenes = 1000,
                    dropout.mid = 0,
                    de.prob = 0.5,
                    # nBatch=1,
                    dropout.shape = -1)
for (j in 1:length(nGroups.list)) {
  for (i in 1:length(nCells.list)) {
    params <- setParams(params, params.list)
    sim[[i]] <- splatSimulate(params = params,
                              batchCells = nCells.list[i],
                              # nCells = nCells.list[i],
                              group.prob = Groups.list[[j]],
                              method = 'group',
                              verbose = FALSE)
    sim.groups <- sim[[i]]
    reference.data <- sim.groups@assays@data@listData$TrueCounts
    dropout.data <- sim.groups@assays@data@listData$Dropout
    cell.type <- sim.groups@colData@listData$Group
    
    obversed.data <- reference.data
    obversed.data[dropout.data == TRUE] = 0
    # after.zeros.num = length(obversed.data[obversed.data ==0])
    # zeros.num = length(reference.data[reference.data ==0])
    # 
    # drop_rate = zeros.num/after.zeros.num
    # i
    # drop_rate
    print(paste('cell',i,sep = '_'))
    write.csv(obversed.data, file = paste('cell',nCells.list[i],'group',length(nGroups.list[[j]]),'obsevered_data.csv',sep = '_'))
    write.csv(reference.data, file = paste('cell',nCells.list[i],'group',length(nGroups.list[[j]]), 'reference_data.csv', sep = '_'))
    write.csv(cell.type, file = paste('cell',nCells.list[i],'group',length(nGroups.list[[j]]),'cell_type.csv', sep = '_'))
    
  }
  saveRDS(sim, file = paste('group',length(nGroups.list[[j]]),'sp_1.rds',sep = '_'))
}


# simulation step 2 -------------------------------------------------------
# different degree of cluster separability
params <- newSplatParams()
sim <-list()
de.prob.list <- c(0.3,0.6,0.9)
params.list <- list(mean.rate = 0.5,
                    seed = 644049,
                    dropout.type = 'experiment',
                    nGenes = 1000,
                    batchCells = 500,
                    group.prob = c(0.25,0.25,0.25,0.25),
                    dropout.mid = 0,
                    dropout.shape = -1)

for (i in c(1:length(de.prob.list))) {
  params <- setParams(params, params.list)
  sim[[i]] <- splatSimulate(params = params,
                            # batchCells = nCells.list[i],
                            # nCells = nCells.list[i],
                            # group.prob = Groups.list[[j]],
                            method = 'group',
                            de.prob = de.prob.list[i],
                            verbose = FALSE)
  sim.groups <- sim[[i]]
  reference.data <- sim.groups@assays@data@listData$TrueCounts
  dropout.data <- sim.groups@assays@data@listData$Dropout
  cell.type <- sim.groups@colData@listData$Group
  
  obversed.data <- reference.data
  obversed.data[dropout.data == TRUE] = 0
  # after.zeros.num = length(obversed.data[obversed.data ==0])
  # zeros.num = length(reference.data[reference.data ==0])
  # 
  # drop_rate = zeros.num/after.zeros.num
  # i
  # drop_rate
  print(paste('cell',i,sep = '_'))
  write.csv(obversed.data, file = paste('de_prob',de.prob.list[i],'obsevered_data.csv',sep = '_'))
  write.csv(reference.data, file = paste('de_prob',de.prob.list[i], 'reference_data.csv', sep = '_'))
  write.csv(cell.type, file = paste('de_prob',de.prob.list[i],'cell_type.csv', sep = '_'))
  
}
saveRDS(sim, file = paste('nde_prob',length(de.prob.list),'sp_2.rds',sep = '_'))


# simulation step3 --------------------------------------------------------
# different degree of dropout rate

params <- newSplatParams()
sim <-list()
dropout.mid.list <- c(2,4,6)
params.list <- list(mean.rate = 0.5,
                    seed = 644049,
                    dropout.type = 'experiment',
                    nGenes = 1000,
                    batchCells = 500,
                    de.prob = 0.5,
                    group.prob = c(0.25,0.25,0.25,0.25),
                    dropout.shape = -1)

for (i in c(1:length(de.prob.list))) {
  params <- setParams(params, params.list)
  sim[[i]] <- splatSimulate(params = params,
                            # batchCells = nCells.list[i],
                            # nCells = nCells.list[i],
                            # group.prob = Groups.list[[j]],
                            method = 'group',
                            dropout.mid = dropout.mid.list[i],
                            verbose = FALSE)
  sim.groups <- sim[[i]]
  reference.data <- sim.groups@assays@data@listData$TrueCounts
  dropout.data <- sim.groups@assays@data@listData$Dropout
  cell.type <- sim.groups@colData@listData$Group
  
  obversed.data <- reference.data
  obversed.data[dropout.data == TRUE] = 0
  # after.zeros.num = length(obversed.data[obversed.data ==0])
  # zeros.num = length(reference.data[reference.data ==0])
  # 
  # drop_rate = zeros.num/after.zeros.num
  # i
  # drop_rate
  print(paste('cell',i,sep = '_'))
  write.csv(obversed.data, file = paste('dropout_prob',dropout.mid.list[i],'obsevered_data.csv',sep = '_'))
  write.csv(reference.data, file = paste('dropout_prob',dropout.mid.list[i], 'reference_data.csv', sep = '_'))
  write.csv(cell.type, file = paste('dropout_prob',dropout.mid.list[i],'cell_type.csv', sep = '_'))
  
}
saveRDS(sim, file = paste('ndropout_mid',length(dropout.mid.list),'sp_2.rds',sep = '_'))