rm(list = ls())
# pairwise correlations

tr = read.table('../../Body/3Results/TandRepInfo.txt', header=TRUE, sep = '\t')

data = tr[, c('Species', 'ConsensusLength', 'CopyNumber', 'PercentMatches', 'fr_A_repeat',
              'fr_T_repeat', 'fr_G_repeat', 'fr_C_repeat')]

summary(data)

pdf('../../Body/4Figures/TRclusters.R.pdf')
plot(data[, -1])

dev.off()

######################################################################################
########## PCA, clusterization

library(stats)
library(gdata)

tr = read.table('../../Body/3Results/TandRepInfo.txt', header=TRUE, sep = '\t')
GenLength = read.xls('../../Body/1Raw/GenerationLengthForMammals.xlsx')

GenLength$Species = gsub(' ', '_', GenLength$Scientific_name)
GenLength = GenLength[, c(14,16)]

mamm = tr[tr$TAXON == 'Mammalia',]
mamm$Consensus = as.character(mamm$Consensus)

pdf('../../Body/4Figures/MammClustersTr.pdf', width = 500, height = 500)

dist_matrix = dist(mamm[, c('ConsensusLength', 'FullLength', 'CopyNumber', 'PercentMatches',
                            'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', "fr_C_repeat", 'Consensus')])

clusters = hclust(dist_matrix)

plot(clusters, labels = mamm$Species)

clusterCut <- cutree(clusters, 5)

clusters_table = cbind(mamm, clusterCut)

dev.off()

###### play with distance matrix

dist_matrix = as.matrix(dist(mamm[, c('ConsensusLength', 'CopyNumber', 'PercentMatches',
                            'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', "fr_C_repeat")]))



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

#######################################################################################
### PCA 

for_pca = clusters_table[, c('Species', 'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', 'fr_C_repeat',
                             'FullLength', 'CopyNumber', 'InDloop', 'PercentMatches')]

# consensus length instead of full length

a = prcomp(for_pca[, -1])

plot(a, type = "l")

summary(a)

for_pca$Pca1 = a$x[, 1]
for_pca$Pca2 = a$x[, 2]
for_pca$Pca3 = a$x[, 3]

for_pca[for_pca$Pca1 < quantile(for_pca$Pca1, 0.05),]$Species 
for_pca[for_pca$Pca1 > quantile(for_pca$Pca1, 0.95),]$Species 

summary(for_pca$Pca1)
summary(for_pca$Pca2)
summary(for_pca$Pca3)

data = merge(for_pca, GenLength, by='Species')

cor.test(data$GenerationLength_d, data$Pca1, method='spearman')
cor.test(data$GenerationLength_d, data$Pca2, method='spearman')
cor.test(data$GenerationLength_d, data$Pca3, method='spearman')

pdf('../../Body/4Figures/PCA.pdf')

par(mfrow=c(2,3))
summary(a)
plot(a)
plot(data$Pca1, data$Pca2)
plot(data$Pca2, data$Pca3)
# plot(PCA$x[,1],MATRIX$GenerationLength_d); cor.test(PCA$x[,1],MATRIX$GenerationLength_d, method = 'spearman') # nothing  - First mutagen signature! Body mass normalized BMR!
biplot(a, col = c('grey','black'), cex = 0.5)
biplot(a, choices=c(2,3), col = c('grey','black'), cex = 0.5) #  biplot(princomp(USArrests),choices=c(1,3))

dev.off()
