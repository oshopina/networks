##Important nodes detection 
load('Data/Shared_data/Basic_metrics.RData')

library(dplyr)

################### Hub and authority scores ########################## 
authority = tables_by_metric$authority_score
authority_threshold = mean(as.matrix(authority[,c(2:8)]), na.rm = T) +
  3 * sd(as.matrix(authority[,c(2:8)]), na.rm = T)

authority_nodes <- list()
for (i in colnames(authority[, 2:8])) {
  imp <- authority %>%
    filter(!!sym(i) > authority_threshold) %>%
    select(row_names)
  authority_nodes[[i]] <- imp
}

page_rank = tables_by_metric$page_rank
page_rank_threshold = mean(as.matrix(page_rank[,c(2:8)]), na.rm = T) +
  sd(as.matrix(page_rank[,c(2:8)]), na.rm = T)

page_rank_nodes <- list()
for (i in colnames(page_rank[, 2:8])) {
  imp <- page_rank %>%
    filter(!!sym(i) > page_rank_threshold) %>%
    select(row_names)
  page_rank_nodes[[i]] <- imp
}

###################### Combination of centrality measures ####################################

all_keystone_otu = list()

for (i in c(4, 45, 5, 55, 6, 65, 7)) {
  degree_var = paste0("degree_", i)
  betweenness_var = paste0("betweenness_", i)
  closeness_var = paste0("closeness_", i)
  
  network_measures <- data.frame(
    OTU = tables_by_metric$degree$row_names,
    Degree = tables_by_metric[['degree']][[degree_var]],
    Betweenness = tables_by_metric[['betweenness']][[betweenness_var]],
    Closeness = tables_by_metric[['closeness']][[closeness_var]]
  )
  
  # Define the thresholds for each network measure
  degree_threshold <- mean(network_measures$Degree, na.rm = T) +
    sd(network_measures$Degree, na.rm = T)
  betweenness_threshold <-
    mean(network_measures$Betweenness, na.rm = T) +
    sd(network_measures$Betweenness, na.rm = T)
  closeness_threshold <-
    mean(network_measures$Closeness, na.rm = T) +
    sd(network_measures$Closeness, na.rm = T)
  
  # Identify keystone OTUs based on high network measures
  keystone_otus <- network_measures[network_measures$Degree >= degree_threshold &
                                      network_measures$Betweenness >= betweenness_threshold &
                                      network_measures$Closeness >= closeness_threshold, ] %>%
    na.omit()
  
  all_keystone_otu[[as.character(i)]] = keystone_otus
}

