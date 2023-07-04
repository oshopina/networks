## Classify all the notes in a network according to roles

load('Data/Shared_data/Basic_metrics.RData')
clustering_results <- readRDS('Data/Shared_data/clustering_results.rds')

library(igraph)

pH_indices <- c(4, 45, 5, 55, 6, 65, 7)

##Z
# Initialize an empty list to store the within-module degree z-scores
within_module_z_scores_all <- list()

# Iterate over the desired module indices

for (i in pH_indices) {
  # Get the graph and clustering results
  net_g_var <- paste0("net_g", i, "_dist")
  net_g <- data_list[[net_g_var]]
  cluster <- clustering_results[['louvain']][[as.character(i)]]
  
  # Get module membership vector
  module_membership <- membership(cluster)
  
  # Calculate within-module degree z-scores for all nodes
  within_module_z_scores <- rep(0, vcount(net_g))
  names(within_module_z_scores) <- V(net_g)$name
  
  for (module_index in 1:max(module_membership)) {
    # Subset the graph to include only nodes in the selected module
    subgraph <- induced_subgraph(net_g, which(module_membership == module_index))
    
    # Calculate the within-module degree for each node in the module
    module_degrees <- degree(subgraph, mode = "in")
    
    module_avg_degree <- mean(module_degrees)
    module_std_degree <- sd(module_degrees)
    
    module_nodes <- which(module_membership == module_index)
    
    within_module_z_scores[module_nodes] <- (module_degrees - module_avg_degree) / module_std_degree
  }
  
  # Store the within-module degree z-scores for the current module
  within_module_z_scores_all[[as.character(i)]] <- within_module_z_scores
}

within_module_z_scores_all <- lapply(within_module_z_scores_all, function(scores) {
  scores[is.na(scores)] <- 0
  return(scores)
})


rm(cluster, net_g, subgraph, i, module_avg_degree, module_degrees, module_index, module_membership,
   module_nodes, module_std_degree, net_g_var, within_module_z_scores)

##P

participation_coefficients_all = list()

for (i in pH_indices) {
  net_g_var <- paste0("net_g", i, "_dist")
  net_g <- data_list[[net_g_var]]
  cluster <- clustering_results[['louvain']][[as.character(i)]]
  
  # Get module membership vector
  module_membership <- membership(cluster)
  
  # Calculate the total degree of each node
  node_degrees <- degree(net_g)
  
  # Create a list to store the participation coefficient for each node
  participation_coefficients <- numeric(length(V(net_g)))
  names(participation_coefficients) <- V(net_g)$name
  
  # Calculate the participation coefficient for each node
  for (node_index in 1:length(V(net_g))) {
    node_module <- module_membership[node_index]
    
    # Get the neighbors of the node
    node_neighbors <- neighbors(net_g, node_index, mode = "total")
    
    # Calculate the number of links to nodes in different modules
    links_to_different_modules <-
      table(module_membership[node_neighbors])
    
    # Calculate the participation coefficient
    participation_coefficient <-
      1 - sum((links_to_different_modules / node_degrees[node_index]) ^ 2)
    
    participation_coefficients[node_index] <-
      participation_coefficient
  }
  participation_coefficients_all[[as.character(i)]] = participation_coefficients
}

rm(cluster, net_g, i, links_to_different_modules, module_membership, net_g_var,
   node_degrees, node_index, node_module, node_neighbors, participation_coefficient,
   participation_coefficients)


##Classification

library(ggplot2)

# Initialize lists to store the classification and percentage results
classification_all <- list()
percentage_all <- list()
scatter_plots <- list()

# Iterate over the desired pH groups
for (pH_group in pH_indices) {
  # Perform classification for each pH group
  classification <- vector("character", length(V(data_list[[paste0("net_g", pH_group, "_dist")]])))
  names(classification) <- V(data_list[[paste0("net_g", pH_group, "_dist")]])$name
  
  for (node_index in 1:length(V(data_list[[paste0("net_g", pH_group, "_dist")]]))) {
    z_score <- within_module_z_scores_all[[as.character(pH_group)]][node_index]
    participation_coefficient <- participation_coefficients_all[[as.character(pH_group)]][node_index]
    
    if (is.na(z_score) || is.na(participation_coefficient)) {
      classification[node_index] <- "Unknown"
    } else if (z_score <= 2.5 && participation_coefficient <= 0.62) {
      classification[node_index] <- "Peripheral"
    } else if (z_score <= 2.5 && participation_coefficient > 0.62) {
      classification[node_index] <- "Connector"
    } else if (z_score > 2.5 && participation_coefficient <= 0.62) {
      classification[node_index] <- "Module Hub"
    } else if (z_score > 2.5 && participation_coefficient > 0.62) {
      classification[node_index] <- "Network Hub"
    }
  }
  
  # Store the classification results in the list
  classification_all[[as.character(pH_group)]] <- classification
  
  # Calculate the percentage of each category
  percentage <- prop.table(table(classification)) * 100
  
  # Store the percentage results in the list
  percentage_all[[as.character(pH_group)]] <- percentage
  
  # Create a data frame for plotting
  df <- data.frame(Classification = classification,
                   Z_Score = within_module_z_scores_all[[as.character(pH_group)]],
                   Participation_Coefficient = participation_coefficients_all[[as.character(pH_group)]])
  
  # Create the scatter plot with bigger point size
  plot <- ggplot(df, aes(x = Participation_Coefficient, y = Z_Score, color = Classification)) +
    geom_point(size = 3) +
    labs(x = "Participation Coefficient (P)", y = "Z-Score (z)", color = "Classification") +
    scale_color_manual(values = c("Peripheral" = "blue",
                                  "Connector" = "green",
                                  "Module Hub" = "red",
                                  "Network Hub" = "purple",
                                  "Unknown" = "gray")) +
    theme_classic() +
    ggtitle(paste("pH Group", pH_group)) +
    scale_x_continuous(limits = c(-0.0000001, 0.9)) +
    scale_y_continuous(limits = c(-3, 5))
  
  # Store the scatter plot in the list
  scatter_plots[[as.character(pH_group)]] <- plot
}

rm(node_index, participation_coefficient, percentage, pH_group, z_score, df, plot)

##Find hubs
module_hubs_all <- list()
network_hubs_all <- list()

# Iterate over the pH groups
for (pH_group in pH_indices) {
  # Get the classification result for the current pH group
  classification <- classification_all[[as.character(pH_group)]]
  
  # Find the indices of module hubs and network hubs
  module_hub_indices <- which(classification == "Module Hub")
  network_hub_indices <- which(classification == "Network Hub")
  
  # Get the names of module hubs and network hubs
  module_hubs <- names(classification)[module_hub_indices]
  network_hubs <- names(classification)[network_hub_indices]
  
  # Store the module hubs and network hubs in the lists
  module_hubs_all[[as.character(pH_group)]] <- module_hubs
  network_hubs_all[[as.character(pH_group)]] <- network_hubs
}

rm(module_hub_indices, module_hubs, network_hub_indices, network_hubs, pH_group)

# # Create a PDF device to save the plots
# pdf("Figures/classification_plots.pdf", width = 8, height = 6)
# for (i in c(4, 45, 5, 55, 6, 65, 7)) {
#   plot_obj <- scatter_plots[[as.character(i)]]
#   print(plot_obj)
# }
# 
# # Close the PDF device
# dev.off()

