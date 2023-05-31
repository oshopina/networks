##Filter fungi from taxonomy 

##libraries
library(dplyr)

##load data 
otu_with_tax = read.csv2('Data/otu_18S_4350_with_tax_Kenneth_update.csv', header = T, sep = ',', 
                             row.names = 1)
otu_with_tax = otu_with_tax[,c(1:141)]

otu_without_fungi = otu_with_tax %>% filter(!grepl('Fungi', taxonomy))

write.csv2(otu_without_fungi, 'Data/otu_18S_4350_without_fungi.csv')
