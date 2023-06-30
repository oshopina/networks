##Important node selection with different methods
## Load necessary libraries and data
load('Data/Shared_data/Basic_metrics.RData')
library(dplyr)

################### Hub and Page-rank scores ########################## 

# Calculate authority scores
authority <- tables_by_metric$authority_score
authority_nodes <- list()

for (i in colnames(authority[, 2:8])) {
  # Calculate threshold for current column
  authority_threshold <- mean(authority[[i]], na.rm = TRUE) +
    3 * sd(authority[[i]], na.rm = TRUE)
  
  # Filter rows based on the threshold for current column
  imp <- authority %>%
    filter(get(i) > authority_threshold) %>%
    select(row_names)
  
  authority_nodes[[i]] <- imp
}

rm(authority, imp, authority_threshold)
#writexl::write_xlsx(authority_nodes, 'Data/authority_nodes.xlsx')

# Calculate PageRank scores
page_rank <- tables_by_metric$page_rank
page_rank_nodes <- list()

for (i in colnames(page_rank[, 2:8])) {
  # Calculate threshold for current column
  column_values <- page_rank[[i]]
  page_rank_threshold <- mean(column_values, na.rm = TRUE) +
    2 * sd(column_values, na.rm = TRUE)
  
  # Filter rows based on the threshold for current column
  imp <- page_rank %>%
    filter(.data[[i]] > page_rank_threshold) %>%
    select(row_names)
  
  page_rank_nodes[[i]] <- imp
}

rm(page_rank, page_rank_threshold, imp)
#writexl::write_xlsx(page_rank_nodes, 'Data/page_rank_nodes.xlsx')

###################### Combination of centrality measures ####################################

all_keystone_otu <- list()

for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  degree_var <- paste0("degree_", i)
  betweenness_var <- paste0("betweenness_", i)
  closeness_var <- paste0("closeness_", i)
  
  # Create a data frame with network measures
  network_measures <- data.frame(
    OTU = tables_by_metric$degree$row_names,
    Degree = tables_by_metric[['degree']][[degree_var]],
    Betweenness = tables_by_metric[['betweenness']][[betweenness_var]],
    Closeness = tables_by_metric[['closeness']][[closeness_var]]
  )
  
  # Define the thresholds for each network measure
  degree_threshold <- mean(network_measures$Degree, na.rm = TRUE) +
    sd(network_measures$Degree, na.rm = TRUE)
  betweenness_threshold <-
    mean(network_measures$Betweenness, na.rm = TRUE) +
    sd(network_measures$Betweenness, na.rm = TRUE)
  closeness_threshold <-
    mean(network_measures$Closeness, na.rm = TRUE) +
    sd(network_measures$Closeness, na.rm = TRUE)
  
  # Identify keystone OTUs based on high network measures
  keystone_otus <- network_measures[network_measures$Degree >= degree_threshold &
                                      network_measures$Betweenness >= betweenness_threshold &
                                      network_measures$Closeness >= closeness_threshold, ] %>%
    na.omit()
  
  all_keystone_otu[[as.character(i)]] <- keystone_otus
}
rm(betweenness_threshold, betweenness_var, closeness_threshold, closeness_var, degree_var, degree_threshold,
   i, keystone_otus)
writexl::write_xlsx(all_keystone_otu, 'Data/keystone_nodes.xlsx')
