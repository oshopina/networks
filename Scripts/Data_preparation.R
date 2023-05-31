## SpiecEasi workflow for 16S, ITS, 18S relative abundance 

##libraries
library(dplyr)

##load data 

otu_table_16s_RA = read.csv2('Data/otu_16S_4500.RA.csv', header = T, sep = ',', 
                             row.names = 1)
otu_table_ITS_RA = read.csv2('Data/otu_ITS_7150.RA.csv', header = T, sep = ',', 
                             row.names = 1)
otu_table_18s_RA = read.csv2('Data/otu_18S_4350.RA.csv', header = T, sep = ',', 
                             row.names = 1)
env_data = read.csv2('Data/env_16S.csv', header = T, sep = ',', 
                     row.names = 1)


##filter tables to the same number of samples

otu_table_16s_RA = select(otu_table_16s_RA, !"H093") ##H093 is an outlier
otu_table_ITS_RA = select(otu_table_ITS_RA, !c("H154","H058","H082","H093"))
otu_table_18s_RA = select(otu_table_18s_RA, !c("H154","H058","H082","H093"))
env_data = env_data[rownames(env_data) %in% colnames(otu_table_16s_RA),]

##prepare data for the analysis
##normilise data
otu_table_16s_RA = otu_table_16s_RA/4500
otu_table_16s_RA = t(otu_table_16s_RA)
rowSums(otu_table_16s_RA)

otu_table_ITS_RA = otu_table_ITS_RA/7150
otu_table_ITS_RA = t(otu_table_ITS_RA)
rowSums(otu_table_ITS_RA)

otu_table_18s_RA = otu_table_18s_RA/4350
otu_table_18s_RA = t(otu_table_18s_RA)
rowSums(otu_table_18s_RA)

##sort data by pH

env_data = env_data[order(env_data$pH),]
ordered_row_names <- rownames(env_data)

otu_table_16s_RA <- otu_table_16s_RA[ordered_row_names, ]
otu_table_ITS_RA <- otu_table_ITS_RA[ordered_row_names, ]
otu_table_18s_RA <- otu_table_18s_RA[ordered_row_names, ]

##choose OTUs higher than 5%

otu_table_16s_RA_slim <- otu_table_16s_RA[,which(apply(otu_table_16s_RA, 2, max)>0.005)]
otu_table_ITS_RA_slim <- otu_table_ITS_RA[,which(apply(otu_table_ITS_RA, 2, max)>0.005)]
otu_table_18s_RA_slim <- otu_table_18s_RA[,which(apply(otu_table_18s_RA, 2, max)>0.005)]

##change otu names so they are unique
colnames(otu_table_16s_RA_slim) <- paste("B.", colnames(otu_table_16s_RA_slim), sep = "")
colnames(otu_table_ITS_RA_slim) <- paste("F.", colnames(otu_table_ITS_RA_slim), sep = "")
colnames(otu_table_18s_RA_slim) <- paste("E.", colnames(otu_table_18s_RA_slim), sep = "")


####################test#################################

##subset to pH 4
env_data_g4 = filter(env_data, pH.group == 'g4')
g4_samples = rownames(env_data_g4)
otu_table_16s_RA_g4 <- otu_table_16s_RA_slim[g4_samples, ]
otu_table_ITS_RA_g4 <- otu_table_ITS_RA_slim[g4_samples, ]
otu_table_18s_RA_g4 <- otu_table_18s_RA_slim[g4_samples, ]

g4_test = list(otu_table_16s_RA_g4, otu_table_ITS_RA_g4, otu_table_18s_RA_g4)

saveRDS(g4_test, 'Data/g4_test.RDS')

