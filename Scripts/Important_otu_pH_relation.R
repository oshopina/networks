
######################### Combine tables ######################################
otu_table_16s_RA = read.csv2(
  'Data/otu_16S_4500.RA.csv',
  header = T,
  sep = ',',
  row.names = 1
)
otu_table_16s_ABS = read.csv2(
  'Data/otu_16S_4500.ABS.csv',
  header = T,
  sep = ',',
  row.names = 1
)
otu_table_ITS_RA = read.csv2(
  'Data/otu_ITS_7150.RA.csv',
  header = T,
  sep = ',',
  row.names = 1
)
otu_table_ITS_RA = otu_table_ITS_RA[, colnames(otu_table_ITS_RA) %in% colnames(otu_table_16s_RA)]
otu_table_ITS_ABS = read.csv2(
  'Data/otu_ITS_7150.ABS.csv',
  header = T,
  sep = ',',
  row.names = 1
)
otu_table_ITS_ABS = otu_table_ITS_ABS[, colnames(otu_table_ITS_ABS) %in% colnames(otu_table_16s_ABS)]

env_data = read.csv2(
  'Data/env_16S.csv',
  header = T,
  sep = ',',
  row.names = 1
)

env_data = env_data[rownames(env_data) %in% colnames(otu_table_16s_RA), ]
env_data = env_data[order(env_data$pH), ]
ordered_row_names <- rownames(env_data)
env_data$pH = as.numeric(env_data$pH)

otu_table_16s_RA <- otu_table_16s_RA[, ordered_row_names]
otu_table_16s_ABS <- otu_table_16s_ABS[, ordered_row_names]
otu_table_ITS_RA <- otu_table_ITS_RA[, ordered_row_names]
otu_table_ITS_ABS <- otu_table_ITS_ABS[, ordered_row_names]

rownames(otu_table_16s_RA) = paste0("B.", rownames(otu_table_16s_RA))
rownames(otu_table_16s_ABS) = paste0("B.", rownames(otu_table_16s_ABS))
rownames(otu_table_ITS_RA) = paste0("F.", rownames(otu_table_ITS_RA))
rownames(otu_table_ITS_ABS) = paste0("F.", rownames(otu_table_ITS_ABS))


RA = rbind(otu_table_16s_RA, otu_table_ITS_RA)
ABS = rbind(otu_table_16s_ABS, otu_table_ITS_ABS)

############################# Choose only important OTUs ##################### 

keystones = readxl::read_xlsx('Data/Results/Important_taxa_with_tax.xlsx')
RA = RA[rownames(RA) %in% keystones$OTU,]
ABS = ABS[rownames(ABS) %in% keystones$OTU,]
otu_tables = list(RA = RA, ABS = ABS)

########################## Correlation coefficients ############################

cor_results_all = list()

for (j in names(otu_tables)) {
  df <- as.data.frame(t(otu_tables[[j]]))
  
  cor_results <- list()
  for (i in rownames(RA)) {
    cor <- cor.test(as.numeric(df[[i]]), env_data$pH, method = 'pearson')
    cor_results[[i]] <- cor
  }
  cor_results_all[[j]] <- cor_results
}

############################ Make graphs #######################################
library(ggplot2)

otu_plots_all = list()

for (j in names(otu_tables)) {
df <- as.data.frame(t(otu_tables[[j]]))
df$pH <- env_data$pH
colnames(df)[colnames(df) == "env_data$pH"] <- "pH"
cor_results = cor_results_all[[j]]
otu_plots <- list()

for (i in rownames(RA)) {
  p = cor_results[[i]][['p.value']]
  r = cor_results[[i]][['estimate']]
  rounded_p_value <- round(p, 4)
  if (rounded_p_value < 0.001) {
    rounded_p_value <- "< 0.001"
  } else {
    rounded_p_value <- rounded_p_value
  }

  if (p > 0.05) {
    color = "gray"
  } else if (r > 0) {
    color = "red"
  } else {
    color = "blue"
  }
  
  plot <- ggplot(df, aes(x = pH, y = .data[[i]])) +
    geom_point() +
    geom_smooth(method = "loess", color = color) +
    annotate("text", x = 4, y = max(df[[i]], na.rm = T),
             label = paste0("p-value = ", rounded_p_value, "; r = ", round(r, 2)),
             hjust = -0.1, vjust = 1, color = "black", size = 4)
  
  otu_plots[[i]] <- plot
}
otu_plots_all[[j]] = otu_plots
}

