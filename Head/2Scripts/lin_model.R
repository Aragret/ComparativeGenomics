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
