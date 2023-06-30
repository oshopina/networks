##Metrics with one number per sample, degree distribution
load('Data/Shared_data/Basic_metrics.RData')

library(igraph)
library(ggplot2)
library(gridExtra)

################## Degree distribution #######################

degree_histograms = list()

for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  net_g_var <- paste0("net_g", i, "_dist")
  net_g <- data_list[[net_g_var]]
  
  # Degree distribution
  degree_distribution <- degree.distribution(net_g)
  
  plot = hist(degree_distribution, breaks = "FD", col = "skyblue", border = "black",
              xlim = c(0, 0.15), ylim = c(0, 50),
              xlab = "Degree", ylab = "Frequency", 
              main = paste("pH", i, "Degree Distribution Histogram"))
  
  degree_histograms[[as.character(i)]] = plot
}

# Set the layout to display plots vertically
par(mfrow = c(length(degree_histograms), 1))
par(mar = c(2, 2, 1, 1))

# Loop through the list of graphs and plot them
for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  plot(degree_histograms[[as.character(i)]], breaks = "FD", col = "skyblue", border = "black",
       xlim = c(0, 0.20), ylim = c(0, 50),
       xlab = "Degree", ylab = "Frequency", 
       main = paste("pH", i, "Degree Distribution Histogram"))
}

########################## One per sample metrics ##############################

path_all = c()
diameter_all = c()
density_all = c()

for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  net_g_var <- paste0("net_g", i, "_dist")
  net_g <- data_list[[net_g_var]]
  
  path = average.path.length(net_g)
  path_all[as.character(i)] = path 

  diameter = diameter(net_g)
  diameter_all[as.character(i)] = diameter
  
  density = graph.density(net_g)
  density_all[as.character(i)] = density
}

one_per_sample = data.frame(path_all, diameter_all, density_all, rownames = names(path_all))

par(mfrow = c(3, 1))
par(mar = c(2, 2, 1, 1))

path_plot = ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = path_all), stat = "identity") +
  labs(title = paste("Average path length")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

diameter_plot = ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = diameter_all), stat = "identity") +
  labs(title = paste("Diameter")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

density_plot = ggplot(one_per_sample, aes(x = rownames)) +
  geom_bar(aes(y = density_all), stat = "identity") +
  labs(title = paste("Density")) +
  scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4'))

grid.arrange(path_plot, diameter_plot, density_plot, nrow = 3)
