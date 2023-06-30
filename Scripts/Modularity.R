##Metrics with one number per sample, distribution, modularity
load('Data/Shared_data/Basic_metrics.RData')

library(igraph)
library(ggplot2)
library(gridExtra)


#############Test different modularity methods #####################

## Define the cluster methods you want to use
##Calculate different cluster methods

cluster_methods <- c("louvain", "infomap", "walktrap", "betweenness", "propagation",
                     "fastgreedy", "leading.eigenvector", "spinglass", "optimal", "leiden")

clustering_results <- list()

for (method in cluster_methods) {
  cat("Running clustering method:", method, "\n")
  
  modularity_c <- c()
  n_modules_c <- c()
  clustering_result <- list()
  

for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  cat("  Running clustering for i =", i, "\n")
  
  net_g_var <- paste0("net_g", i, "_dist")
  net_g <- data_list[[net_g_var]]
  
  cluster_result <- switch(method,
                           louvain = cluster_louvain(net_g),
                           infomap = cluster_infomap(net_g),
                           walktrap = cluster_walktrap(net_g),
                           betweenness = cluster_edge_betweenness(net_g),
                           propagation = cluster_label_prop(net_g),
                           fastgreedy = cluster_fast_greedy(net_g),
                           leading.eigenvector = cluster_leading_eigen(net_g),
                           spinglass = cluster_spinglass(net_g),
                           optimal = cluster_optimal(net_g),
                           leiden = cluster_leiden(net_g))
  
  # Store clustering result in a list
  clustering_result[[as.character(i)]] <- cluster_result
}
  clustering_results[[method]] <- clustering_result
}

saveRDS(clustering_results, "Data/Shared_data/clustering_results.rds")

# Create an empty list to store the plots and clustering results
modularity_plots <- list()

# Iterate over each cluster method
for (method in cluster_methods) {
  cluster = clustering_results[[method]]
  cat("  Running clustering for method", method, "\n")
  
  for (i in c(4, 45, 5, 55, 6, 65, 7)) {
    cat("  Running clustering for i =", i, "\n")
    cluster_result = cluster[[as.character(i)]]
    
    modularity <- modularity(cluster_result)
    n_modules <- length(cluster_result)
    
    modularity_c[as.character(i)] <- modularity
    n_modules_c[as.character(i)] <- n_modules
}
  
  modularity_test <- data.frame(modularity_c, n_modules_c)
  modularity_test$sample <- rownames(modularity_test)
  
  plot_mod <- ggplot(modularity_test, aes(x = sample)) +
    geom_bar(aes(y = modularity_c), stat = "identity") +
    labs(title = paste("Modularity -", method)) +
    scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))
  
  plot_n <- ggplot(modularity_test, aes(x = sample)) +
    geom_bar(aes(y = n_modules_c), stat = "identity") +
    labs(title = paste("Number of Modules -", method)) +
    scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))
  
  combined_plot <- grid.arrange(plot_mod, plot_n, nrow = 1)
  
  # Save the plot and clustering results to the respective lists
  modularity_plots[[method]] <- combined_plot
}

pdf("Figures/modularity_plots.pdf")
for (plot in modularity_plots) {
  print(plot)
}
dev.off()

modularity_plots = list()
for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  net_dist_var <- paste0("net_g", i, "_dist")
  net_abs_var = paste0("net_g", i, "_abs")
  net_g_dist <- data_list[[net_dist_var]]
  net_g_abs <- data_list[[net_abs_var]]
  
  cluster_greedy = cluster_fast_greedy(net_g_abs)
  modularity = modularity(cluster_greedy)
  n_modules = length(cluster_greedy)
  l = layout.fruchterman.reingold(net_g_abs)
  community_colors = rainbow(n_modules)
  
  plot(net_g_abs, layout = l, vertex.size = 5, 
       vertex.color = community_colors[membership(cluster_greedy)],
       vertex.label = NA, edge.color = "gray", main = paste0('pH ', i))
  
  # Add text for modularity and number of clusters at the middle bottom
  mtext(side = 1, line = 2, paste("Modularity:", round(modularity, 2)), col = "black")
  mtext(side = 1, line = 1, paste("Number of Clusters:", n_modules), col = "black")
  
  plot = recordPlot()
  modularity_plots[[i]] = plot
}

