library(SpiecEasi)
library(igraph)

env = read.csv('Data/env_16S.csv')
rownames(env) = env$Sample
bacteria = read.csv('Data/otu_16S_4500.RA.csv')
fungi = read.csv('Data/otu_ITS_7150.RA.csv')
rownames(bacteria) = bacteria$X
rownames(fungi) = fungi$X

bacteria = bacteria[,-1]
fungi = fungi[,-1]

fungi = fungi[,colnames(bacteria)]

rownames(bacteria) = paste0("B.", rownames(bacteria))
rownames(fungi) = paste0("F.", rownames(fungi))

################################ SpiecEasi #####################################

# Define a list of groups and their criteria

groups <- list(
  list(x16s = 2, its = 2, name = "45"),
  list(x16s = 3, its = 1, name = "65"),
  list(x16s = 1, its = 3, name = "7")
)

bacteria_group_tables = list()
fungi_group_tables = list()
env_group_tables = list()

#Loop through the groups
for (i in 1:length(groups)) {
  group = groups[[i]]
  
  env_group <- env[env$X16S_clustering == group$x16s & env$ITS_clustering == group$its,]
  
  bacteria_group <- bacteria[, colnames(bacteria) %in% env_group$Sample]
  bacteria_group <- bacteria_group[which(apply(bacteria_group, 1, max) > 22.5),] |> t()
  
  fungi_group <- fungi[, colnames(fungi) %in% env_group$Sample]
  fungi_group <- fungi_group[which(apply(fungi_group, 1, max) > 35.75),] |> t()
  
  bacteria_group_tables[[i]] = bacteria_group
  fungi_group_tables[[i]] = fungi_group
  env_group_tables[[i]] = env_group
  
  # Perform the network analysis and save the result (uncomment this section when needed)
  # spieceasi <- spiec.easi(list(bacteria_group, fungi_group), 
  #                         method = 'mb', lambda.min.ratio = 1e-2, nlambda = 19, 
  #                         icov.select.params = list(rep.num = 50, ncores = 10))
  # 
  # saveRDS(spieceasi, paste0("network_", group$name, ".rds"))
}


############################ Network preparation ###############################
matrixes = list()
nets =  list()
nets_dist = list()
nets_abs = list()


for (i in 1:length(groups)) {
  group = groups[[i]]
  bacteria_group = bacteria_group_tables[[i]]
  fungi_group = fungi_group_tables[[i]]
  
  spiec = readRDS(paste0('Data/SpiecEasi_3_groups/network_', group$name, '_05.rds'))
  matrix = symBeta(getOptBeta(spiec), mode='maxabs') |> as.matrix()
  colnames(matrix) = c(colnames(bacteria_group), colnames(fungi_group))
  rownames(matrix) = c(colnames(bacteria_group), colnames(fungi_group))
  
  net = graph.adjacency(matrix,mode = "undirected", weighted = TRUE, diag = FALSE)
  V(net)$name = colnames(matrix)
  
  net_dist <- net
  E(net_dist)$weight <- 1 - abs(E(net_dist)$weight)
  
  net_abs <- net
  E(net_abs)$weight <- abs(E(net_abs)$weight)
  
  E(net)[weight>0]$color <- "#68C2A3"
  E(net)[weight<0]$color <- "#B51945"
  
  V(net)$color <- ifelse(substr(V(net)$name, 1, 1) == "B", "black", "white")
  
  matrixes[[i]] = matrix
  nets[[i]] =  net
  nets_dist[[i]] = net_dist
  nets_abs[[i]] = net_abs
  
}

############################ Network graph #####################################
library(grid)
library(gridGraphics)
library(ggplotify)
library(patchwork)
library(ggplot2)
 
l = layout_with_mds(nets_dist[[1]])
  
 plot(nets[[1]], edge.color = E(nets[[1]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[1]])$colors, vertex.label = "", rescale=F, layout=l*0.02)
 
 grid.echo()
 p1 <- grid.grab()
 p1 = ggplotify::as.ggplot(p1)

l = layout_with_mds(nets_dist[[2]])

 plot(nets[[2]], edge.color = E(nets[[2]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[2]])$colors, vertex.label = "", rescale=F, layout=l*0.02)
 
 grid.echo()
 p2 <- grid.grab()
 p2 = ggplotify::as.ggplot(p2)

l = layout_with_mds(nets_dist[[3]])

 plot(nets[[3]], edge.color = E(nets[[3]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[3]])$colors, vertex.label = "", rescale=F, layout=l*0.02)

 grid.echo()
 p3 <- grid.grab()
 p3 = ggplotify::as.ggplot(p3)

########################## Fake legend ########################################

 data <- data.frame(
   Year = 2000:2010,
   Value_A = cumsum(runif(11)),
   Value_B = cumsum(runif(11))
 )

dummy_plot = ggplot(data, aes(x = Year)) +
  geom_line(aes(y = Value_A, color = 'Positive')) +
  geom_line(aes(y = Value_B, color = 'Negative')) + 
  scale_color_manual(values = c("#B51945", "#68C2A3"), name = 'Association') +
  theme_bw()
 
 l1 = cowplot::get_legend(dummy_plot)
 l1 = as.ggplot(l1)
 
 dummy_plot2 <- ggplot(data, aes(x = Year)) +
   geom_point(aes(y = Value_A, fill = 'Bacteria'), shape = 21, size = 3) +
   geom_point(aes(y = Value_B, fill = 'Fungi'), shape = 21, size = 3) +
   scale_fill_manual(values = c('Bacteria' = 'black', 'Fungi' = 'white'), name = 'Guild') +
   theme_bw()
 
 l2 = cowplot::get_legend(dummy_plot2)
 l2 = as.ggplot(l2)
 
############################### Final plot ####################################

layout = c(
  area(t = 1, b = 10, l = 1, r = 10),
  area(t = 1, b = 10, l = 10, r = 20),
  area(t = 1, b = 10, l = 20, r = 30),
  area(t = 7, b = 9, l = 29, r = 30)
)

l = c(area(t = 1, b = 15, l = 1, r = 5),
       area(t = 10, b = 15, l = 1, r = 5))
 
legend = l1 + theme(panel.spacing = unit(0, "lines")) + l2 + plot_layout(design = l) 
 
final_plot = p1 + ggtitle('pH 3.7~4.5') + theme(plot.title = element_text(face = "bold", size = 20)) + 
   p2 + ggtitle('pH 4.5~6.1') + theme(plot.title = element_text(face = "bold", size = 20)) + 
   p3 + ggtitle('pH 6.1~8.0') + theme(plot.title = element_text(face = "bold", size = 20)) +
   legend + plot_layout(design = layout)

# saveRDS(final_plot, 'Data/Shared_data/3_plots_05.rds')

# ggsave('Figures/NETWORK_05.png', 
#        final_plot, 
#        device = 'png', 
#        width = 21,
#        height = 10)
