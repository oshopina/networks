library(dplyr)
library(igraph)
library(ggplot2)

all_networks = readRDS('Data/Shared_data/full_networks_igraph.rds')

edge_weights = data.frame()

for (i in 1:3) {
  network = all_networks[['net']][[i]]
  edge_weight = data_frame(weight = E(network)$weight, network = i)
  edge_weights = rbind(edge_weights, edge_weight)
}

rm(edge_weight, network, i)

edge_weights$Association = ifelse(edge_weights$weight >= 0, "Positive", "Negative")

summary_data <- edge_weights %>%
  group_by(network, Association) %>%
  summarise(count = n())

positive = summary_data %>% tidyr::pivot_wider(names_from = Association, values_from = count, values_fill = 0) %>%
mutate(Positive_Percentage = Positive / (Positive + Negative) * 100)

a <-  (10000 - 0) / (70 - 50) 
b <- 0 - a * 50

final_plot = ggplot() +
  geom_bar(data = summary_data, aes(x = as.factor(network), y = count, fill = Association), stat = "identity", position = "stack", color = 'black') +
  geom_line(data = positive, aes(x = network, y = a * Positive_Percentage + b), linewidth = 1.2, color = '#951c2a') +
  scale_fill_manual(values = c("Positive" = "#eb66c1", "Negative" = "#5caeda")) +
  scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(~ (. - b) / a, name = "Positive Interactions (%)")) +
  labs(x = "pH", y = 'No. interactions') +
  theme_classic() +
  scale_x_discrete(labels = c('3.7~4.5', '4.5~6.1', '6.1~8.0'))

#ggsave('Figures/interactions.png', final_plot, device = 'png')

