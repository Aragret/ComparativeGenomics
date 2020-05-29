rm(list = ls(all=TRUE))

if (!require(seqinr)) install.packages("seqinr")

library(seqinr)

rep = read.table('../../Body/3Results/TandRepInfo.txt', header = TRUE, sep='\t')

toSave = rep[,c('Species', 'Consensus')]

write.fasta(as.list(toSave$Consensus), toSave$Species, file.out = '../../Body/2Derived/TrConsensuses.fasta')

################################################################

names(rep)

mamm_rep = rep[rep$TAXON == 'Mammalia',]
mamm_rep$GC_cons = (mamm_rep$fr_G_cons + mamm_rep$fr_C_cons) / (mamm_rep$fr_G_cons + mamm_rep$fr_C_cons +
                                                             mamm_rep$fr_A_cons + mamm_rep$fr_T_cons)

for(i in 1:nrow(mamm_rep)){
  if(mamm_rep[i,]$CopyNumber > median(mamm_rep$CopyNumber)){
    mamm_rep[i, 'CopyNumberDummy'] = 1
  }
  if(mamm_rep[i,]$CopyNumber <= median(mamm_rep$CopyNumber)){
    mamm_rep[i, 'CopyNumberDummy'] = 0
  }
  if(mamm_rep[i,]$PercentMatches > median(mamm_rep$PercentMatches)){
    mamm_rep[i, 'PercentMatchesDummy'] = 1
  }
  if(mamm_rep[i,]$PercentMatches <= median(mamm_rep$PercentMatches)){
    mamm_rep[i, 'PercentMatchesDummy'] = 0
  }
  if(mamm_rep[i,]$GC_cons > median(mamm_rep$GC_cons)){
    mamm_rep[i, 'GC_consDummy'] = 1
  }
  if(mamm_rep[i,]$GC_cons <= median(mamm_rep$GC_cons)){
    mamm_rep[i, 'GC_consDummy'] = 0
  }
  if(mamm_rep[i,]$ConsensusLength > median(mamm_rep$ConsensusLength)){
    mamm_rep[i, 'ConsensusLengthDummy'] = 1
  }
  if(mamm_rep[i,]$ConsensusLength <= median(mamm_rep$ConsensusLength)){
    mamm_rep[i, 'ConsensusLengthDummy'] = 0
  }
}

write.table(mamm_rep, '../../Body/2Derived/TandRepInfoMammals.txt', sep='\t', quote = FALSE, row.names = FALSE)