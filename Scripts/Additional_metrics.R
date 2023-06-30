##Metrics with one number per sample, distribution, modularity
load('Data/Shared_data/Basic_metrics.RData')

library(igraph)
library(ggplot2)
library(gridExtra)


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
