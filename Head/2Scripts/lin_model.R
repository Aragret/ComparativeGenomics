rm(list = ls())
library(gdata)

# pairwise correlations

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
GenLength = read.xls('../../Body/1Raw/GenerationLengthForMammals.xlsx')

GenLength$Species = gsub(' ', '_', GenLength$Scientific_name)
GenLength = GenLength[, c(14,16)]

data = merge(CHOR, GenLength, by='Species')

summary(data$REP.DirRepLength)

data = data[, c('Species', 'A', 'T', 'G', 'C', 'REP.DirRepLength', 'GenomeLength', 'GenerationLength_d')]

data$FrA = data$A / data$GenomeLength
data$FrT = data$T / data$GenomeLength
data$FrG = data$G / data$GenomeLength
data$FrC = data$C / data$GenomeLength

data$DRCoverage = data$REP.DirRepLength / data$GenomeLength

summary(lm(log(data$GenerationLength_d) ~ data$DRCoverage + data$FrA + data$FrT + data$FrG))

################################################################################

cor.test(data$GenerationLength_d, data$DRCoverage, method = 'spearman')
# rho -0.1519978, pvalue 5.011e-05

summary(lm(data$DRCoverage ~ data$FrA + data$FrT + data$FrG))

# (Intercept)  -0.5108     0.1480  -3.452 0.000590 ***
#   data$FrA      0.9762     0.2995   3.260 0.001170 ** 
#   data$FrT      0.7906     0.1495   5.287 1.66e-07 ***
#   data$FrG      1.5478     0.4581   3.379 0.000768 ***

par(mfrow=c(2, 2))
hist(data$FrA, breaks = 50)
hist(data$FrT, breaks = 50)
hist(data$FrG, breaks = 50)
hist(data$FrC, breaks = 50)

summary(data$FrA)
summary(data$FrT)
summary(data$FrG)
summary(data$FrC)
