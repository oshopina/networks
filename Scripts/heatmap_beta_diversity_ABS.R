library(dplyr)
library(ggplot2)

####################### 16S + all OTUs ######################################

otu16S <- read.csv('Data/otu_16S_4500.ABS.csv', row.names = 1)
env16S <- read.csv('Data/env_16S.csv')

# Filter and remove outliers in one step
env16S <- env16S[env16S$Sample %in% colnames(otu16S) & env16S$Sample != "H093", ]
env16S <- env16S[order(env16S$pH), ]

# Subset otu16S based on env16S$Sample
otu16S <- otu16S[, colnames(otu16S) %in% env16S$Sample]
otu16S <- as.data.frame(t(otu16S))

###################ITS + all OTUs ############################################

otuITS <- read.csv('Data/otu_ITS_7150.ABS.csv', row.names = 1)
otuITS[otuITS < 0] = 0
otuITS = otuITS[, !colnames(otuITS) %in% c("H057","H046", "H047", "H041")] ## remove samples with only negative values
envITS <- read.csv('Data/env_16S.csv')

# Filter envITS and order by pH
envITS <- envITS[envITS$Sample %in% colnames(otuITS), ]
envITS <- envITS[order(envITS$pH), ]

# Subset otuITS based on envITS$Sample
otuITS <- otuITS[, colnames(otuITS) %in% envITS$Sample]
otuITS <- as.data.frame(t(otuITS))


############################# Beta_diversity ###################################

hellinger_diversity = function(otu_table) {
  dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
       method = 'euc') %>% as.matrix() 
  return(dist_matrix)
}

beta_diversity_16S = hellinger_diversity(otu16S)
beta_diversity_ITS = hellinger_diversity(otuITS)

beta_diversity_16S = beta_diversity_16S[env16S$Sample, env16S$Sample]
beta_diversity_ITS = beta_diversity_ITS[envITS$Sample, envITS$Sample]

############################### Heatmap #######################################
library(ComplexHeatmap)


distance_heatmap = function(dist_matrix, env_table, names_for_pH, graph_name) {
  
  my_palette = colorRampPalette(c('#100d12', '#980011', '#fe8d2f', '#fff1a2'))
  col_fun = circlize::colorRamp2(c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
                                 c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a",
                                   "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2"))
  
  ha = HeatmapAnnotation(empty = anno_empty(border = FALSE, height = unit(3, "mm")),
                         pH = env_table$pH, empty2 = anno_empty(border = FALSE, height = unit(3, "mm")),
                         pH_labels = anno_text(names_for_pH, rot = 0),
                         col = list(pH = col_fun), show_legend = F,
                         gp = gpar(col = "black"))
  ha2 = rowAnnotation(pH_labels = anno_text(names_for_pH, rot = 0), pH = env_table$pH,
                      col = list(pH = col_fun), show_legend = F, 
                      show_annotation_name = F, gp = gpar(col = "black"))
  
  heatmap = Heatmap(dist_matrix, row_order = rownames(dist_matrix), column_order = colnames(dist_matrix),
                    col = my_palette(100), show_row_names = F, show_column_names = F, show_heatmap_legend = F,
                    bottom_annotation = ha, left_annotation = ha2, column_title = graph_name)
  
  return(heatmap)
}

names_16S = c(3.7, rep("", 18), 4, rep("", 19), 4.5, rep("", 15), 5, rep("", 15), 5.5, rep("", 8), 6, rep("", 12), 
          6.5, rep("", 14), 7, rep("", 17), 7.5, rep("", 8), 8.0)
names_ITS = c(3.7, rep("", 19), 4, rep("", 21), 4.5, rep("", 15), 5, rep("", 15), 5.5, rep("", 9), 6, rep("", 12), 
              6.5, rep("", 14), 7, rep("", 17), 7.5, rep("", 4), 8.0)

library(ggplotify)

heatmap_16S = as.ggplot(distance_heatmap(beta_diversity_16S, env16S, names_16S, "Bacteria (all OTUs)"))
heatmap_ITS = as.ggplot(distance_heatmap(beta_diversity_ITS, envITS, names_ITS, "Fungi (all OTUs)"))

library(patchwork)

legend = Legend(col_fun = circlize::colorRamp2(c(0, 0.5, 1, 1.5), c(
  '#100d12', '#980011', '#fe8d2f', '#fff1a2'
)),
title = "Hellinger distance")
gglegend = grid.grabExpr(draw(legend))

layout <- layout <- '
AAA#
AAA#
AAAC
BBBC
BBB#
BBB#
'


plot_heat = heatmap_16S + heatmap_ITS  + gglegend +
  plot_layout(design = layout)

ggsave('Figures/tipping_pounts_heatmap_ABS.svg', plot_heat,
       device = 'svg', height = 10, width = 7)
