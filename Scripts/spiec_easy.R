bacteria = read.csv('Data/otu_16S_4500.RA.csv')
fungi = read.csv('Data/otu_ITS_7150.RA.csv')
rownames(bacteria) = bacteria$X
rownames(fungi) = fungi$X

bacteria = bacteria[,-1]
fungi = fungi[,-1]

fungi = fungi[,colnames(bacteria),]
