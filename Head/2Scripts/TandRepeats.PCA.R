library(kernlab)

reps = read.table('../../Body/3Results/TandRepInfo.txt', header = TRUE, sep='\t')

MATRIX = as.matrix(reps[, c('CopyNumber', 'PercentMatches', 'FullLength', 'ConsensusLength',
                            'fr_A_cons', 'fr_T_cons', 'fr_G_cons', 'fr_C_cons',
                            'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', 'fr_C_repeat',
                            'Number', 'InDloop')])

PCA = prcomp(MATRIX, center = TRUE) 
print(PCA)  
summary(PCA)

PCA$rotation

# PC1
# FullLength      0.8738773488
# fr_A_repeat     0.3121346054
# fr_T_repeat     0.2949352413

# PC2
# ConsensusLength 0.8271678067
# fr_A_cons 0.2765043819
# fr_T_cons 0.2386549047

MATRIX$Pca1 = PCA$x[,1]
MATRIX$Pca2 = PCA$x[,2]
MATRIX$Pca3 = PCA$x[,3]
MATRIX$Pca4 = PCA$x[,4]

pdf('../../Body/4Figures/TandRepeats.PCA.pdf')

biplot(PCA, choices=c(1,2), col = c('white','black'), cex = 0.8) #  biplot(princomp(USArrests),choices=c(1,3))

plot(MATRIX$Pca1, MATRIX$Pca2)

dev.off()
