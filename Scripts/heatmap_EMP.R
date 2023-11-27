library(dplyr)
library(ggplot2)
library(vegan)

#library(qiime2R)
####################### 16S + all OTUs #########################################
# otu16S <- read_q2biom('../../EMP/emp_cr_silva_16S_123.qc_filtered.rare_10000.biom')
# env16S <- read.csv('../../EMP/emp_qiime_mapping_qc_filtered.tsv', sep = '\t')
# 
# # Filter samples to soil samples
# 
# env16S = env16S[env16S$empo_3 == "Soil (non-saline)",]
# env16S = env16S[!is.na(env16S$ph),]
# env16S <- env16S[order(env16S$ph), ]
# otu16S = as.data.frame(otu16S)
# otu16S = otu16S[,env16S$X.SampleID]
# otu16S <- as.data.frame(t(otu16S))

load('Data/Shared_data/EMP.RData')

env16S = env16S[env16S$X.SampleID %in% rownames(otu16S),]
#env16S = env16S[env16S$depth_m < 0.1,]
otu16S = otu16S[rownames(otu16S) %in% env16S$X.SampleID,]

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
env16S$colors = mypal_pH(256)[as.numeric(cut(env16S$ph, breaks=256))]
env16S = mutate(env16S, pH_group = case_when(
  ph < 4.5 ~ "Group 1",
  between(ph, 4.5, 6) ~ "Group 2",
  ph >= 6 ~ "Group 3"
))
env16S$pH_group = as.factor(env16S$pH_group)
env16S$study_id = as.factor(env16S$study_id)

otu16S = otu16S[,!colSums(otu16S) == 0]
#otu16S <- otu16S[, apply(otu16S, 2, max) >= 100]

# hellinger_diversity = function(otu_table) {
#   dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
#                      method = 'euc') %>% as.matrix()
#   return(dist_matrix)
# }
# 
# beta_diversity_16S = hellinger_diversity(otu16S)
#  
# # save.image('EMP2.RData')
# 
# library(ComplexHeatmap)
# 
# load('Data/Shared_data/EMP2.RData') ##Data for all OTUs
# 
# beta_diversity_16S = beta_diversity_16S[env16S$X.SampleID,]
# beta_diversity_16S = beta_diversity_16S[,env16S$X.SampleID]
# 
# my_palette = colorRampPalette(c('#100d12', '#980011', '#fe8d2f', '#fff1a2'))
# Heatmap(beta_diversity_16S, row_order = rownames(beta_diversity_16S), column_order = colnames(beta_diversity_16S),
#         col = my_palette(100), show_row_names = F, show_column_names = F, show_heatmap_legend = F)
# 
# otu_heat = vegan::decostand(otu16S[, apply(otu16S, 2, max) >= 1000], 'standardize', 1)
# otu_heat = otu_heat[env16S$X.SampleID,]
# heatmap(t(otu_heat), revC = TRUE, 
#         Colv = NA, Rowv = NA)

DCA = decorana(otu16S, iweigh = 1)
ordiplot(DCA, display = 'sites', type = 'p')
points(DCA, col = env16S$study_id, pch = 20, cex = 2)
legend("topright", legend = levels(env16S$study_id), col = unique(env16S$study_id), pch = 20, title = "Study id")
text(DCA, labels = env16S$ph)
