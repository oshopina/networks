library(vegan)
library(dplyr)
library(ggplot2)


otu16S = read.csv('Data/otu_16S_4500.RA.csv', row.names = 1)
env16S = read.csv('Data/env_16S.csv')
env16S = env16S[order(env16S$pH),]
env16S = env16S[env16S$Sample %in% colnames(otu16S),]

### remove outlier
env16S <- env16S[-which(env16S$Sample %in% c("H093")), ]
otu16S <- otu16S[,-which(colnames(otu16S) %in% c("H093")) ]

otu16S = otu16S[, env16S$Sample]
otu16S = as.data.frame(t(otu16S))
otu16S = otu16S[, apply(otu16S, 2, max) > 22]

## Beta-diversity measures

beta_methods <- c("w", "-1", "c", "wb", "r", "I", "e", "t", "me", "j", "sor", "m", "-2", "co", "cc", "g", "-3", "l", "19", "hk", "rlb", "sim", "gl", "z")

heatmap_list <- list()

for (method in beta_methods) {
  beta_diversity <- betadiver(otu16S, method)
  
  # Convert the beta diversity data to a distance matrix
  dist_matrix <- as.matrix(beta_diversity)
  
  # Convert the distance matrix to a long format suitable for ggplot
  dist_long <- reshape2::melt(dist_matrix)
  
  # Create a custom color palette that goes through white
  my_palette <- colorRampPalette(c("black", "green"))(100)
  
  # Create the heatmap using ggplot2
  heatmap_plot <- ggplot(dist_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = my_palette) +
    labs(title = paste("Beta Diversity Heatmap - Method", method), x = "pH Value", y = "pH Value") +
    theme_minimal() +
    scale_x_discrete(labels = env16S$pH) +  # Replace sample names with pH values on x-axis
    scale_y_discrete(labels = env16S$pH) + # Replace sample names with pH values on y-axis
    theme(axis.text.x = element_text(angle=90))
  
  heatmap_list[[method]] <- heatmap_plot
}

heatmap_list$j

# Save the heatmaps or further manipulate them as needed