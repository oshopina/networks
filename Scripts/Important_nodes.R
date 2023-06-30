##Important nodes detection 
load('Data/Shared_data/Basic_metrics.RData')


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