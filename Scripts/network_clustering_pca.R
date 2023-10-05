library(dplyr)

##Instead of running the loop above load ready to go results of it
clustering_results = readRDS('Data/Shared_data/clustering_results.rds')
clustering_results = clustering_results$louvain

library(reshape2)

presence_tables = list()

for (i in c('4', '45', '5', '55', '6', '65', '7')) {
    cluster_result = clustering_results[[i]]
    df = data.frame(cluster = cluster_result$membership, OTU = cluster_result$names)
    df$value = 1
    df = dcast(df, OTU ~ cluster, value.var = "value", fill = 0)
    rownames(df) = df$OTU
    df = df[,-1]
    colnames(df) <- paste0("cluster_", i, "_", colnames(df))
    df$OTU = rownames(df)
    presence_tables[[i]] = df
}

combined_table = Reduce(function(df1, df2) merge(df1, df2, by = "OTU", all = TRUE), presence_tables)
combined_table[is.na(combined_table)] <- 0

library("FactoMineR")
library("factoextra")
pca_table = combined_table[,-1]
rownames(pca_table) = combined_table$OTU
pca_table = as.data.frame(t(pca_table))

network = stringr::str_extract_all(rownames(pca_table), "(?<=_)[0-9]+(?=_)") %>% unlist()


pca = PCA(pca_table, scale.unit = F)
pca_table$pc1 <- pca$ind$coord[, 1] # indexing the first column
pca_table$pc2 <- pca$ind$coord[, 2]  # indexing the second column

full_plot = ggplot(data = pca_table, aes(
  x = pc1,
  y = pc2,
  color = network,
  label = rownames(pca_table)
)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  geom_text()

small_plot = ggplot(data = pca_table, aes(
  x = pc1,
  y = pc2,
  color = network,
  label = rownames(pca_table)
)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  geom_text() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5))

full_plot_no_labels = ggplot(data = pca_table, aes(
  x = pc1,
  y = pc2,
  color = network
)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8) +
  theme_minimal() 

small_plot_no_labels = ggplot(data = pca_table, aes(
  x = pc1,
  y = pc2,
  color = network
)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  xlim(c(-2.5, 2.5)) +
  ylim(c(-2.5, 2.5))

# pdf("Figures/pca_network_clusters.pdf", width = 10, height = 6)
# print(full_plot)
# print(full_plot_no_labels)
# print(small_plot)
# print(small_plot_no_labels)
# dev.off()

