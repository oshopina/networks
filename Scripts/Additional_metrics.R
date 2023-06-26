## Script for calculating new metrics for network analysis 
library(SpiecEasi)
library(igraph)

load("Data/SpiecEasi_results/MA_postSpiecEasi_g4.ABS.Rdata")

#extract the adjacency matrix from the spiec.easi object
spieceasi.matrix.g4 <- symBeta(getOptBeta(spieceasi.multinet.g4), mode='maxabs')
spieceasi.matrix.g4 <- as.matrix(spieceasi.matrix.g4)

#add correct row names from the otu table
otu.names.g4 <- c(colnames(otu.table.16S.g4.slim), colnames(otu.table.ITS.g4.slim))
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

# greedy method (hiearchical, fast method)
c1 = cluster_fast_greedy(net.g4.abs)

# modularity measure
modularity(c1)
membership(c1)
length(c1)
sizes(c1)

l = layout_in_circle(net.g4.abs)

plot(c1, net.g4.abs, layout = l)
plot_dendrogram(c1)

hubs = hub_score(net.g4.abs)
sort(hubs$vector)
sort(degree(net.g4.abs))
sort(closeness(net.g4.dist))
sort(betweenness(net.g4.dist,v = V(net.g4.dist)))

a = cluster_spinglass(net.g4)

# Calculate degree centrality
degree <- degree(net.g4)

# Calculate betweenness centrality
betweenness <- betweenness(net.g4.dist)

# Define the roles based on centrality measures and community membership
peripherals <- V(net.g4)$name[degree <= quantile(degree, 0.25)]
connectors <- V(net.g4)$name[betweenness >= quantile(betweenness, 0.75) & degree > quantile(degree, 0.25)]
module_hubs <- V(net.g4)$name[degree > quantile(degree, 0.75) & !V(net.g4)$name %in% connectors]
network_hubs <- V(net.g4)$name[degree > quantile(degree, 0.75)]

# Calculate closeness centrality
closeness <- closeness(net.g4.dist)

# Calculate eigenvector centrality
eigenvector <- eigen_centrality(net.g4)$vector

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

# Calculate the clustering coefficient
clustering_coef <- transitivity(net.g4, type = "local")
average_path_length <- average.path.length(net.g4.dist)

# Degree distribution
degree_distribution <- degree.distribution(net.g4)

# Assortativity
assortativity <- assortativity(net.g4, types1 = "numeric")

# Diameter
diameter <- diameter(net.g4)

# Density
density <- graph.density(net.g4)

# Centralization
centralization <- centralization.degree(net.g4)

# Average clustering coefficient
average_clustering_coef <- average.path.length(net.g4)
                                           
