
reps = read.table('../../Body/3Results/TandRepInfo.txt', header = TRUE, sep='\t')

MATRIX = as.matrix(reps[, c('CopyNumber', 'PercentMatches', 'FullLength', 'ConsensusLength',
                            'fr_A_cons', 'fr_T_cons', 'fr_G_cons', 'fr_C_cons',
                            'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', 'fr_C_repeat')])

row.names(MATRIX) = reps$Species

MATRIX = scale(MATRIX)

PCA = prcomp(MATRIX) 
print(PCA)  
summary(PCA)

PCA$rotation


MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]

pdf('../../Body/4Figures/TandRepeats.PCA.pdf')

biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))

plot(MATRIX$Pca1, MATRIX$Pca2)
plot(MATRIX$Pca2, MATRIX$Pca3)
plot(MATRIX$Pca3, MATRIX$Pca4)


dev.off()
