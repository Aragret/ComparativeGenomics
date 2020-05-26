rm(list = ls(all=TRUE))

if (!require(seqinr)) install.packages("seqinr")

library(seqinr)

rep = read.table('../../Body/3Results/TandRepInfo.txt', header = TRUE, sep='\t')

toSave = rep[,c('Species', 'Consensus')]

write.fasta(as.list(toSave$Consensus), toSave$Species, file.out = '../../Body/2Derived/TrConsensuses.fasta')
