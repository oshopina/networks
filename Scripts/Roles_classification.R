## Classify all the notes in a network according to roles
library(igraph)

load('Data/Shared_data/Basic_metrics.RData')
clustering_results <- readRDS('Data/Shared_data/clustering_results.rds')

Zs <- list()
Ps <- list()
nodes_class <- list()

for (j in c(4, 45, 5, 55, 6, 65, 7)) {
  net_g_var <- paste0("net_g", j, "_dist")
  net_g <- data_list[[net_g_var]]
  community <- clustering_results[['infomap']][[as.character(j)]]
  
  # Calculate within-module connectivity (z)
  membership <- membership(community)
  z <- numeric(vcount(net_g))
  for (i in 1:vcount(net_g)) {
    module <- membership[i]
    module_nodes <- which(membership == module)
    within_module_degree <- sum(degree(net_g)[module_nodes])
    average_degree <- mean(degree(net_g)[module_nodes])
    std_degree <- sd(degree(net_g)[module_nodes])
    z[i] <- (within_module_degree - average_degree) / std_degree
  }
  Zs[[as.character(j)]] <- z
  
  # Calculate participation coefficient (P)
  P <- numeric(vcount(net_g))
  for (i in 1:vcount(net_g)) {
    neighbors <- neighbors(net_g, i)
    module_neighbors <- unique(membership[neighbors])
    P[i] <- 1 - sum((degree(net_g)[neighbors] / sum(degree(net_g)[membership == module_neighbors]))^2)
  }
  Ps[[as.character(j)]] <- P
  
  # Classification thresholds
  z_threshold <- 2.5
  P_threshold <- 0.62
  
  # Create a vector to store the node classifications
  node_class <- rep(NA, vcount(net_g))
  
  # Classify nodes
  for (i in 1:vcount(net_g)) {
    if (!is.na(z[i]) && !is.na(P[i])) {
      if (z[i] <= z_threshold && P[i] <= P_threshold) {
        node_class[i] <- "Peripheral"
      } else if (z[i] <= z_threshold && P[i] > P_threshold) {
        node_class[i] <- "Connector"
      } else if (z[i] > z_threshold && P[i] <= P_threshold) {
        node_class[i] <- "Module Hub"
      } else {
        node_class[i] <- "Network Hub"
      }
    } else {
      node_class[i] <- NA
    }
  }
  
  nodes_class[[as.character(j)]] <- node_class
  
  # Print the number of nodes in each class for the current pH sample
  print(table(node_class))
}


