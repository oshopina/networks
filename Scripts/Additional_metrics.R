##Metrics with one number per sample, distribution, modularity
load('Data/Shared_data/Basic_metrics.RData')


##Modularity

modularity_plots = list()
for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  net_dist_var <- paste0("net_g", i, "_dist")
  net_g_dist <- data_list[[net_dist_var]]
  
  cluster_greedy = cluster_fast_greedy(net_g_dist)
  modularity = modularity(cluster_greedy)
  n_modules = length(cluster_greedy)
  l = layout.fruchterman.reingold(net_g_dist)
  community_colors = rainbow(n_modules)
  
  plot(net_g_dist, layout = l, vertex.size = 5, 
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
