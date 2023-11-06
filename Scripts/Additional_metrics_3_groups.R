##Metrics with one number per sample, degree distribution
## Load necessary libraries and data
load('Data/Shared_data/Basic_metrics_3_groups.RData')
library(igraph)
library(ggplot2)
library(gridExtra)

################## Degree distribution #######################

# Calculate degree distribution for different pH values
degree_histograms <- list()

for (i in 1:3) {
  net_g <- all_networks[["net"]][[i]]
  
  # Degree distribution
  degree_distribution <- degree.distribution(net_g)
  
  # Create a histogram plot
  plot <- hist(degree_distribution, breaks = "FD", col = "skyblue", border = "black",
               xlim = c(0, 0.15), ylim = c(0, 50),
               xlab = "Degree", ylab = "Frequency", 
               main = paste("pH", i, "Degree Distribution Histogram"))
  
  degree_histograms[[as.character(i)]] <- plot
}
rm(net_g, plot, degree_distribution, i)

## Plot degree distribution
# Set the layout to display plots vertically
par(mfrow = c(length(degree_histograms), 1))
par(mar = c(2, 2, 1, 1))

# Loop through the list of plots and display them
for (i in 1:3) {
  plot(degree_histograms[[i]], breaks = "FD", col = "skyblue", border = "black",
       xlim = c(0, 0.20), ylim = c(0, 50),
       xlab = "Degree", ylab = "Frequency", 
       main = paste("pH", i, "Degree Distribution Histogram"))
}
rm(i)

########################## One per sample metrics ##############################

## Calculate metrics
path_all <- c()
diameter_all <- c()
density_all <- c()
clustering_all = c()

for (i in 1:3) {
  net_g <- all_networks[["dist"]][[i]]
  
  # Calculate average path length
  path <- average.path.length(net_g)
  path_all[as.character(i)] <- path 
  
  # Calculate diameter
  diameter <- diameter(net_g)
  diameter_all[as.character(i)] <- diameter
  
  # Calculate density
  density <- graph.density(net_g)
  density_all[as.character(i)] <- density
  
  # Calculate clustering coefficient
  clustering = transitivity(net_g, "global")
  clustering_all[as.character(i)] = clustering 
}

one_per_sample <- data.frame(path_all, diameter_all, density_all, clustering_all,
                             rownames = names(path_all))

rm(density, density_all, diameter, diameter_all, i, path, path_all, net_g, 
   clustering_all, clustering)

## Plot metrics

# Create plots for each metric
path_plot <- ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = path_all), stat = "identity") +
  labs(title = paste("Average path length")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

diameter_plot <- ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = diameter_all), stat = "identity") +
  labs(title = paste("Diameter")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

density_plot <- ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = density_all), stat = "identity") +
  labs(title = paste("Density")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

clustering_plot = ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = clustering_all), stat = "identity") +
  labs(title = paste("Clustering coefficient")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

# Arrange plots in a grid
one_per_sample_plot = grid.arrange(path_plot, diameter_plot, density_plot, clustering_plot,
                                   nrow = 4)

rm(density_plot, diameter_plot, path_plot, clustering_plot)
