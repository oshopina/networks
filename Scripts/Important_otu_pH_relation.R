
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

############################ Make graphs #######################################
library(ggplot2)

df <- as.data.frame(t(RA))
df$pH <- env_data$pH
colnames(df)[colnames(df) == "env_data$pH"] <- "pH"
df$pH <- as.numeric(df$pH)

otu_plots <- list()

for (i in rownames(RA)) {
  plot <- ggplot(df, aes(x = pH, y = .data[[i]])) +
    geom_point() +
    geom_smooth(method = "loess")
  otu_plots[[i]] <- plot
}


