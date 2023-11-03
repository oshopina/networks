library(SpiecEasi)

env = read.csv('Data/env_16S.csv')
rownames(env) = env$Sample
bacteria = read.csv('Data/otu_16S_4500.RA.csv')
fungi = read.csv('Data/otu_ITS_7150.RA.csv')
rownames(bacteria) = bacteria$X
rownames(fungi) = fungi$X

bacteria = bacteria[,-1]
fungi = fungi[,-1]

fungi = fungi[,colnames(bacteria)]

rownames(bacteria) = paste0("B.", rownames(bacteria))
rownames(fungi) = paste0("F.", rownames(fungi))

############################ Group below 4.5 ##################################

env_45 = env[env$X16S_clustering == 2 & env$ITS_clustering == 2,]

bacteria_45 = bacteria[, colnames(bacteria) %in% env_45$Sample]
bacteria_45 = bacteria_45[which(apply(bacteria_45, 1, max)>22.5),] |> t()

fungi_45 = fungi[, colnames(fungi) %in% env_45$Sample]
fungi_45 = fungi_45[which(apply(fungi_45, 1, max)>35.75),] |> t()

spieceasi <- spiec.easi(list(bacteria_45, fungi_45), 
                        method = 'mb',lambda.min.ratio = 1e-2,nlambda = 19, 
                        icov.select.params = list(rep.num = 50, ncores = 10))

saveRDS(spieceasi, 'network_45.rds')

############################ Group below 6.5 ##################################

env_65 = env[env$X16S_clustering == 3 & env$ITS_clustering == 1,]

bacteria_65 = bacteria[, colnames(bacteria) %in% env_65$Sample]
bacteria_65 = bacteria_65[which(apply(bacteria_65, 1, max)>22.5),] |> t()

fungi_65 = fungi[, colnames(fungi) %in% env_65$Sample]
fungi_65 = fungi_65[which(apply(fungi_65, 1, max)>35.75),] |> t()

spieceasi <- spiec.easi(list(bacteria_65, fungi_65), 
                        method = 'mb',lambda.min.ratio = 1e-2,nlambda = 19, 
                        icov.select.params = list(rep.num = 50, ncores = 10))

saveRDS(spieceasi, 'network_65.rds')

############################ Group above 7 ##################################

env_7 = env[env$X16S_clustering == 1 & env$ITS_clustering == 3,]

bacteria_7 = bacteria[, colnames(bacteria) %in% env_7$Sample]
bacteria_7 = bacteria_7[which(apply(bacteria_7, 1, max)>22.5),] |> t()

fungi_7 = fungi[, colnames(fungi) %in% env_7$Sample]
fungi_7 = fungi_7[which(apply(fungi_7, 1, max)>35.75),] |> t()

spieceasi <- spiec.easi(list(bacteria_7, fungi_7), 
                        method = 'mb',lambda.min.ratio = 1e-2,nlambda = 19, 
                        icov.select.params = list(rep.num = 50, ncores = 10))

saveRDS(spieceasi, 'network_7.rds')
