load("C:/Users/uqoshopi/Dropbox/Hoosfield/6 Amplicon/Networks/envs_ABS/MA_postSpiecEasi_g7.ABS.Rdata")

## List RDS files in the directory
files <- list.files(path = "Data/SpiecEasi_results_16S_ITS/", pattern = "\\.rds$", full.names = TRUE)
data_list <- list()
for (file in files) {
  data <- readRDS(file)
  data_list[[file]] <- data
}

#extract the adjacency matrix from the spiec.easi object
spieceasi.matrix.g4 <- symBeta(getOptBeta(data_list$`Data/SpiecEasi_results_16S_ITS/spieceasi_multinet_g7.rds`), 
                               mode='maxabs')
spieceasi.matrix.g4 <- as.matrix(spieceasi.matrix.g4)

#add correct row names from the otu table
otu.names.g4 <- c(colnames(otu.table.16S.g7.slim), colnames(otu.table.ITS.g7.slim))
rownames(spieceasi.matrix.g4) <- otu.names.g4
colnames(spieceasi.matrix.g4) <- otu.names.g4

#build a weighted network from the adjacency matrix
net.g4 <- graph.adjacency(spieceasi.matrix.g4,mode = "undirected", weighted = TRUE, diag = FALSE)
V(net.g4)$name <- otu.names.g4

#convert the weighted network to a separate weighted network
net.g4.dist <- net.g4
max(abs(E(net.g4.dist)$weight))
weights.g4.dist <-  1 - abs(E(net.g4.dist)$weight)
E(net.g4.dist)$weight <- weights.g4.dist

#convert the weighted network to a separate absolute network
net.g4.abs <- net.g4
E(net.g4.abs)$weight <- abs(E(net.g4.abs)$weight)

saveRDS(net.g4, 'Data/SpiecEasi_results_16S_ITS/net_g7.rds')
saveRDS(net.g4.abs, 'Data/SpiecEasi_results_16S_ITS/net_g7_abs.rds')
saveRDS(net.g4.dist, 'Data/SpiecEasi_results_16S_ITS/net_g7_dist.rds')

