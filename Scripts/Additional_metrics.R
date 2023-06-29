##Metrics with one number per sample, distribution, modularity
load('Data/Shared_data/Basic_metrics.RData')

library(igraph)
library(ggplot2)
library(gridExtra)


#############Test different modularity methods #####################

# Define the cluster methods you want to use

cluster_methods <- c("louvain", "infomap", "walktrap", "betweenness", "propagation",
                     "fastgreedy", "leading.eigenvector", "spinglass", "optimal", "leiden")

# Create an empty list to store the plots and clustering results
modularity_plots <- list()
clustering_results <- list()

# Iterate over each cluster method
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
    
    modularity <- modularity(cluster_result)
    n_modules <- length(cluster_result)
    
    modularity_c[as.character(i)] <- modularity
    n_modules_c[as.character(i)] <- n_modules
    
    # Store clustering result in a list
    clustering_result[[as.character(i)]] <- cluster_result
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
  clustering_results[[method]] <- clustering_result
}

pdf("Figures/modularity_plots.pdf")
for (plot in modularity_plots) {
  print(plot)
}
dev.off()

saveRDS(clustering_results, "Data/Shared_data/clustering_results.rds")



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

average.path.length(data_list$net_g4_dist)
# Diameter
diameter <- diameter(data_list$net_g4_dist)
# Density
density <- graph.density(data_list$net_g4_dist)

##Important nodes


# Define the roles based on centrality measures and community membership
peripherals <- V(net.g4)$name[degree <= quantile(degree, 0.25)]
connectors <- V(net.g4)$name[betweenness >= quantile(betweenness, 0.75) & degree > quantile(degree, 0.25)]
module_hubs <- V(net.g4)$name[degree > quantile(degree, 0.75) & !V(net.g4)$name %in% connectors]
network_hubs <- V(net.g4)$name[degree > quantile(degree, 0.75)]



# Combine network measures into a data frame
network_measures <- data.frame(
  Degree = degree,
  Betweenness = betweenness,
  Closeness = closeness,
  Eigenvector = eigenvector
)

# Define the thresholds for each network measure
degree_threshold <- 15
betweenness_threshold <- 1500
closeness_threshold <- 0
eigenvector_threshold <- 0.70

# Identify keystone OTUs based on high network measures
keystone_otus <- row.names(network_measures[
  network_measures$Degree >= degree_threshold &
    network_measures$Betweenness >= betweenness_threshold &
    network_measures$Closeness >= closeness_threshold &
    network_measures$Eigenvector >= eigenvector_threshold,
])

# Create a color palette for community visualization
community_colors <- rainbow(15)

# Create a layout using the Fruchterman-Reingold algorithm
layout <- layout_with_kk(net.g4.dist)

# Plot the community network
plot(net.g4, layout = layout, vertex.size = 5, vertex.color = community_colors[membership(c1)],
     vertex.label = NA, edge.color = "gray", main = "Community Network")

# Add a legend for the community colors
legend("topright", legend = 1:num_communities, fill = community_colors, title = "Community")





# Degree distribution
degree_distribution <- degree.distribution(net.g4)
