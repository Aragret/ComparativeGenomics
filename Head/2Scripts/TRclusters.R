library(stats)
library(gdata)

tr = read.table('../../Body/3Results/TandRepInfo.txt', header=TRUE, sep = '\t')
GenLength = read.xls('../../Body/1Raw/GenerationLengthForMammals.xlsx')

GenLength$Species = gsub(' ', '_', GenLength$Scientific_name)
GenLength = GenLength[, c(14,16)]

mamm = tr[tr$TAXON == 'Mammalia',]
mamm$Consensus = as.character(mamm$Consensus)

pdf('~/Desktop/MammClustersTr.pdf', width = 150, height = 150)

dist_matrix = dist(mamm[, c('ConsensusLength', 'FullLength', 'CopyNumber', 'PercentMatches',
                            'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', "fr_C_repeat", 'Consensus')])

clusters = hclust(dist_matrix)

plot(clusters, labels = mamm$Species)

clusterCut <- cutree(clusters, 5)

clusters_table = cbind(mamm, clusterCut)

dev.off()

############################################

table(clusters_table$clusterCut)

data = merge(clusters_table, GenLength, by="Species")

cor.test(data[data$clusterCut == 1,]$FullLength, data[data$clusterCut == 1,]$GenerationLength_d,
         method = 'spearman')

cor.test(data[data$clusterCut == 2,]$FullLength, data[data$clusterCut == 2,]$GenerationLength_d,
         method = 'spearman')

cor.test(data[data$clusterCut == 3,]$FullLength, data[data$clusterCut == 3,]$GenerationLength_d,
         method = 'spearman')

cor.test(data[data$clusterCut == 4,]$FullLength, data[data$clusterCut == 4,]$GenerationLength_d,
         method = 'spearman')

cor.test(data[data$clusterCut == 5,]$FullLength, data[data$clusterCut == 5,]$GenerationLength_d,
         method = 'spearman')


third = data[data$clusterCut == 3,]

summary(third)
