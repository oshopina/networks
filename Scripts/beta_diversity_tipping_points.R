library(vegan)
library(dplyr)
library(ggplot2)


tipping_points_graph = function(otu_table, env_table, graph_name) {
  #### Sliding window
  # Define pH ranges for windows
  pH_ranges <- seq(3.61, 8.11, by = 0.25)
  
  # calculate the midpoint pH between windows
  midpoint_pH <- (pH_ranges[-1] + pH_ranges[-length(pH_ranges)]) / 2
  
  # assign pH windows
  env_table$pH_window <-
    cut(env_table$pH, breaks = pH_ranges, labels = midpoint_pH)
  
  # combine OTU data and environmental data
  combined_data <- cbind(otu_table, env_table)
  
  # create a list to store dissimilarity matrices for each pH window
  dissimilarity_matrices <- list()
  
  # Loop through pH windows
  for (i in unique(combined_data$pH_window)) {
    # Subset data for the current pH window
    subset_data <- combined_data[combined_data$pH_window == i,]
    
    # Filter data to exclude pH-related columns
    filtered_data <- subset_data[, 1:ncol(otu_table)]
    
    # calculate dissimilarity using Hellinger distance
    dissimilarity <-
      dist(decostand(filtered_data, method = 'hellinger'),
           method = 'euc')
    
    # Store dissimilarity matrix in the list
    dissimilarity_matrices[[as.character(i)]] <- dissimilarity
  }
  
  # calculate mean dissimilarity for each pH window
  mean_dissimilarity_values <-
    sapply(dissimilarity_matrices, function(mat)
      mean(mat))
  # calculate standard errors for each pH window
  se_values <- sapply(dissimilarity_matrices, function(mat) {
    if (length(mat) <= 1) {
      0  # Set SE to 0 if there are too few samples
    } else {
      sqrt(var(mat) / length(mat))
    }
  })
  sd_values <- sapply(dissimilarity_matrices, function(mat) {
    if (length(mat) <= 1) {
      0  # Set SD to 0 if there are too few samples
    } else {
      sd(mat)
    }
  })
  
  # create a data frame for plotting with standard errors
  plot_data_with_se <-
    data.frame(
      pH_window = unique(combined_data$pH_window),
      Dissimilarity = mean_dissimilarity_values,
      SE = se_values,
      SD = sd_values
    )
  
  # Plotting with standard error envelope
  plot <- ggplot(plot_data_with_se, aes(x = pH_window, y = Dissimilarity)) +
    geom_line(group = 1) +
    geom_point() +
    geom_ribbon(aes(x = as.numeric(pH_window), ymin = Dissimilarity - SE, ymax = Dissimilarity + SE),
                alpha = 0.3, fill = "blue") +
    labs(y = "Mean Hellinger distance") +
    scale_x_discrete(expand = c(0.01, 0)) +
    theme_minimal() +
    ggtitle(graph_name) +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  
  return(plot)
}

####################### 16S + all OTUs ######################################

otu16S = read.csv('Data/otu_16S_4500.Ra.csv', row.names = 1)
env16S = read.csv('Data/env_16S.csv')
env16S = env16S[env16S$Sample %in% colnames(otu16S),]

### remove outlier
env16S <- env16S[-which(env16S$Sample %in% c("H093")), ]
otu16S <- otu16S[,-which(colnames(otu16S) %in% c("H093")) ]

otu16S = otu16S[, env16S$Sample]
otu16S = as.data.frame(t(otu16S))

all_16S_graph = tipping_points_graph(otu16S, env16S, graph_name = '16S Hellinger distanaes all OTUs')


####################16S + OTUs > 1% ############################################

otu16S_big <- otu16S[, apply(otu16S, 2, max) > 45]

more_1_16S_graph = tipping_points_graph(otu16S_big, env16S, graph_name = '16S Hellinger distanaes OTUs bigger than 1 %')


###################ITS + all OTUs ############################################

otuITS = read.csv('Data/otu_ITS_7150.Ra.csv', row.names = 1)
envITS = read.csv('Data/env_16S.csv')
envITS = envITS[envITS$Sample %in% colnames(otuITS),]

otuITS = otuITS[, envITS$Sample]
otuITS = as.data.frame(t(otuITS))

all_ITS_graph = tipping_points_graph(otuITS, envITS, 'ITS Hellinger distanaes all OTUs')


################ ITS + 1% OTU > #############################################
otuITS_big <- otuITS[, apply(otuITS, 2, max) > 71.5]

more_1_ITS_graph = tipping_points_graph(otuITS_big, envITS, graph_name = 'ITS Hellinger distanaes OTUs bigger than 1%')


################################ Gradient ####################################

mypal_pH <- colorRampPalette(c("#9e0142",
                               "#d53e4f",
                               "#f46d43",
                               "#fdae61",
                               "#fee08a",
                               "#e6f598",
                               "#aadda4",
                               "#66a2a5",
                               "#3288ad",
                               "#5e4fa2"))

otu_pH <- sort(env16S$pH)
env_16S_sorted = env16S[order(env16S$pH),]
names(otu_pH) = env_16S_sorted$Sample

names = c(3.7, rep("", 18), 4, rep("", 19), 4.5, rep("", 15), 5, rep("", 15), 5.5, rep("", 8), 6, rep("", 12), 
          6.5, rep("", 14), 7, rep("", 17), 7.5, rep("", 8), 8.0)
names_equal_distanae = c(3.7, rep("", 6), 4, rep("", 15), 4.5, rep("", 15), 5, rep("", 16), 5.5, rep("", 16),
                         6, rep("", 14), 6.5, rep("", 15), 7, rep("", 15), 7.5, rep("", 14), 8)
library(ggplotify)
library(ComplexHeatmap)

gradient <- as.ggplot(pheatmap(
  as.matrix(t(otu_pH)),
  cluster_rows = FALSE,  # Don't aluster rows
  cluster_cols = FALSE,  # Don't aluster aolumns
  col = mypal_pH(256),   # austom aolor palette
  border_color = 'black',
  legend = F,
  labels_col = names_equal_distanae, # Remove row laaels
  angle_col = '0'
))


################################ Line graph ##################################
library(patchwork)
line_graph = all_16S_graph + scale_y_continuous(limits = c(0.2, 0.8)) + more_1_16S_graph + 
  scale_y_continuous(limits = c(0.2, 0.8)) + theme(axis.title.y = element_blank()) +
  gradient + gradient +
  all_ITS_graph + scale_y_continuous(limits = c(0.4, 0.8)) + more_1_ITS_graph + 
  scale_y_continuous(limits = c(0.4, 0.8)) + theme(axis.title.y = element_blank()) +
  gradient + gradient + plot_layout(ncol = 2, nrow = 4, heights = c(9, 1))


# ggsave('Figures/tipping_pounts_line.png',
#        line_graph,
#        device = 'png',
#        width = 15,
#        height = 10)


############################# Heatmap ########################################

beta_diversity <-
  dist(decostand(otu16S, method = 'hellinger'),
       method = 'euc') %>% as.matrix() 

beta_diversity_big <-
  dist(decostand(otu16S_big, method = 'hellinger'),
       method = 'euc') %>% as.matrix() 

my_palette = colorRampPalette(c('#100d12', '#980011', '#fe8d2f', '#fff1a2'))


col_fun = circlize::colorRamp2(c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
                               c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a",
                                 "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2"))

beta_diversity = beta_diversity[env_16S_sorted$Sample, env_16S_sorted$Sample]

ha = HeatmapAnnotation(pH = env_16S_sorted$pH, pH_labels = anno_text(names, rot = 0), 
                       col = list(pH = col_fun), show_legend = F,
                             gp = gpar(col = "black"))
ha2 = rowAnnotation(pH_labels = anno_text(names, rot = 0), pH = env_16S_sorted$pH,
                    col = list(pH = col_fun), show_legend = F, 
                    show_annotation_name = F, gp = gpar(col = "black"))

Heatmap(beta_diversity, row_order = rownames(beta_diversity), column_order = colnames(beta_diversity),
        col = my_palette(100), show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = 'Hellinger distance'),
        bottom_annotation = ha, left_annotation = ha2, column_title = '16S Hellinger distanaes all OTUs')

beta_diversity_big = beta_diversity_big[env_16S_sorted$Sample, env_16S_sorted$Sample]

Heatmap(beta_diversity_big, row_order = rownames(beta_diversity_big), column_order = colnames(beta_diversity_big),
        col = my_palette(100), show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = 'Hellinger distance'),
        bottom_annotation = ha, left_annotation = ha2, column_title = '16S Hellinger distanaes all OTUs')

################################ ITS ##########################################
beta_diversity_ITS <-
  dist(decostand(otuITS, method = 'hellinger'),
       method = 'euc') %>% as.matrix() 

beta_diversity_ITS_big <-
  dist(decostand(otuITS_big, method = 'hellinger'),
       method = 'euc') %>% as.matrix() 

env_ITS_sorted = envITS[order(envITS$pH),]

beta_diversity_ITS  = beta_diversity_ITS[env_ITS_sorted$Sample, env_ITS_sorted$Sample]
beta_diversity_ITS_big  = beta_diversity_ITS_big[env_ITS_sorted$Sample, env_ITS_sorted$Sample]

names_ITS = c(3.7, rep("", 19), 4, rep("", 21), 4.5, rep("", 15), 5, rep("", 15), 5.5, rep("", 9), 6, rep("", 12), 
              6.5, rep("", 14), 7, rep("", 17), 7.5, rep("", 8), 8.0)

ha_ITS = HeatmapAnnotation(pH = env_ITS_sorted$pH, pH_labels = anno_text(names_ITS, rot = 0), 
                       col = list(pH = col_fun), show_legend = F,
                       gp = gpar(col = "black"))
ha2_ITS = rowAnnotation(pH_labels = anno_text(names_ITS, rot = 0), pH = env_ITS_sorted$pH,
                    col = list(pH = col_fun), show_legend = F, 
                    show_annotation_name = F, gp = gpar(col = "black"))

Heatmap(beta_diversity_ITS, row_order = rownames(beta_diversity_ITS), column_order = colnames(beta_diversity_ITS),
        col = my_palette(100), show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = 'Hellinger distance'),
        bottom_annotation = ha_ITS, left_annotation = ha2_ITS, column_title = 'ITS Hellinger distances all OTUs')

Heatmap(beta_diversity_ITS_big, row_order = rownames(beta_diversity_ITS_big), column_order = colnames(beta_diversity_ITS_big),
        col = my_palette(100), show_row_names = F, show_column_names = F,
        heatmap_legend_param = list(title = 'Hellinger distance'),
        bottom_annotation = ha_ITS, left_annotation = ha2_ITS, column_title = 'ITS Hellinger distances all OTUs')
