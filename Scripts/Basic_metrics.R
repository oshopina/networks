## Script for calculating basic metrics for network analysis 

## Load required libraries
library(SpiecEasi)
library(igraph)
library(centiserve)
library(dplyr)
library(ggplot2)

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
rm(data, file_name, file)

## Calculate centrality metrics for each network
otu_metrics <- list()

## Loop through the files
for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  
  # Construct the variable names based on the g_number
  net_var <- paste0("net_g", i)
  net_dist_var <- paste0("net_g", i, "_dist")
  net_abs_var <- paste0("net_g", i, "_abs")
  
  # Get the network objects using the constructed variable names
  net_g <- data_list[[net_var]]
  net_g_dist <- data_list[[net_dist_var]]
  net_g_abs <- data_list[[net_abs_var]]
  
  # Calculate the centrality metrics
  result <- data.frame(
    degree = degree(net_g),
    alpha_centrality = alpha.centrality(net_g),
    strength = strength(net_g),
    betweenness = betweenness(net_g_dist, v = V(net_g_dist)),
    closeness = closeness(net_g),
    transitivity = transitivity(net_g, type = 'localundirected'),
    eigen_centrality = eigen_centrality(net_g_dist)$vector,
    page_rank = page.rank(net_g_dist)$vector,
    bottleneck = bottleneck(net_g_dist, v = V(net_g_dist)),
    authority_score = authority_score(net_g)$vector,
    hub_score = hub_score(net_g)$vector,
    centralization = centralization.degree(net_g)$res
  )
  
  # Save the resulting list for the file
  otu_metrics[[as.character(i)]] <- result
}
rm(net_var, net_dist_var, net_abs_var, net_g, net_g_dist, net_g_abs, result, i)

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
rm(table_names, table_list, merged_table, i)

violin_plots = list()
for (i in names(tables_by_metric)) {
  df_long <- tidyr::pivot_longer(tables_by_metric[[i]], 
                                 cols = -row_names, names_to = "sample", values_to = "value")
  anova_result = anova(lm(value ~ sample, data = df_long))
  p_value = anova_result$`Pr(>F)`[1]
  
  # Round the p-value and check if it's equal to 0
  rounded_p_value <- round(p_value, 4)
  if (rounded_p_value == 0) {
    rounded_p_value <- "< 0.001"
  } else {
    rounded_p_value <- paste("= ", rounded_p_value)
  }
  # Plot the violin plot
  plot <- ggplot(df_long, aes(x = sample, y = value, fill = sample)) +
    geom_violin(trim = 0) +
    xlab("pH") +
    ylab(i) +
    theme_classic() +
    scale_fill_manual(values = c("#B51945", '#D13A4C', "#F89151",
                                 "#FDCF7D", "#AEDEA1", "#68C2A3",
                                 "#388FB8")) +
    scale_x_discrete(labels = c('3.9', '4.2', '4.8', '5.4', '6.5', '6.9', '7.4')) +
    guides(fill = 'none') +
    annotate("text", x = 1, y = max(df_long$value, na.rm = T),
             label = paste("p-value", rounded_p_value, "(ANOVA)"),
             hjust = -0.1, vjust = 1, color = "black", size = 4)
  
  # Add the plot to the list
  violin_plots[[i]] <- plot
  }
rm(df_long, plot, anova_result, i, p_value, rounded_p_value)
    
# pdf("Figures/violin_plots.pdf")
# for (plot in violin_plots) {
#   print(plot)
# }
# dev.off()
