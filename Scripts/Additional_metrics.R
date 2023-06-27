## Script for calculating new metrics for network analysis 

## Load required libraries
library(SpiecEasi)
library(igraph)
library(centiserve)
library(dplyr)

## List RDS files in the directory
files <- list.files(path = "Data/SpiecEasi_results_16S_ITS/", pattern = "^net*", full.names = TRUE)
data_list <- list()

## Read RDS files and store in a list
for (file in files) {
  data <- readRDS(file)
  file_name <- basename(file)
  file_name <- tools::file_path_sans_ext(file_name)
  data_list[[file_name]] <- data
}

## Calculate centrality metrics for each network
otu_metrics <- list()

## Loop through the files
for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  # Extract the number after "g" in the filename
  g_number <- i
  
  # Construct the variable names based on the g_number
  net_var <- paste0("net_g", g_number)
  net_dist_var <- paste0("net_g", g_number, "_dist")
  net_abs_var <- paste0("net_g", g_number, "_abs")
  
  # Get the network objects using the constructed variable names
  net_g <- data_list[[net_var]]
  net_g_dist <- data_list[[net_dist_var]]
  net_g_abs <- data_list[[net_abs_var]]
  
  # Calculate the centrality metrics
  result <- data.frame(
    degree = degree(net_g),
    alpha_centrality = alpha.centrality(net_g),
    strength = strength(net_g_abs),
    betweenness = betweenness(net_g_dist, v = V(net_g_dist)),
    closeness = closeness(net_g_dist),
    transitivity = transitivity(net_g, type = 'localundirected'),
    eigen_centrality = eigen_centrality(net_g_dist)$vector,
    page_rank = page.rank(net_g_dist)$vector,
    bottleneck = bottleneck(net_g_dist, v = V(net_g_dist)),
    authority_score = authority_score(net_g)$vector,
    hub_score = hub_score(net_g)$vector
  )
  
  # Save the resulting list for the file
  otu_metrics[[as.character(i)]] <- result
}

## Merge tables by metric
tables_by_metric <- list()

# Iterate over the column names of the '4' table in otu_metrics
for (i in colnames(otu_metrics$'4')) {
  table_names <- c('4', '45', '5', '55', '6', '65', '7')
  
  # Extract the column of interest from each table
  table_list <- lapply(table_names, function(table_num) {
    otu_metrics[[table_num]] %>% select(matches(i))
  })
  
  # Add row names as a column to each table
  table_list <- lapply(table_list, function(table) {
    table$row_names <- rownames(table)
    table
  })
  
  # Set the names of the new list based on the original table names
  names(table_list) <- table_names
  
  # Rename the metric column in each table with the corresponding suffix
  table_list <- lapply(names(table_list), function(table_name) {
    suffix <- gsub("[^0-9]", "", table_name)  # Extract the numeric suffix from the table name
    table <- table_list[[table_name]]
    colnames(table)[colnames(table) == i] <- paste0(i, "_", table_name)
    table
  })
  
  # Merge the columns of all tables by row names with custom suffixes for duplicate columns
  merged_table <- Reduce(function(x, y) merge(x, y, by = "row_names", all = TRUE),
                         table_list)
  
  # Store the merged table in tables_by_metric
  tables_by_metric[[i]] <- merged_table
}




##Important nodes



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

a = cluster_spinglass(net.g4)


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
                                           
