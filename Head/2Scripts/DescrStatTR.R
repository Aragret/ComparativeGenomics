rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

tr = read.table('../../Body/2Derived/TRFinder.txt', header=TRUE, sep='\t')
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header = TRUE, sep='\t')

tr$FullLength = tr$End - tr$Start
tr$ConsensusLength = as.numeric(lapply(as.character(tr$Consensus), nchar))

### nucleotide freq

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

for(i in 1:nrow(tr)){
  tr$fr_A_cons[i] = countCharOccurrences('A', as.character(tr$Consensus[i]))
  tr$fr_T_cons[i] = countCharOccurrences('T', as.character(tr$Consensus[i]))
  tr$fr_G_cons[i] = countCharOccurrences('G', as.character(tr$Consensus[i]))
  tr$fr_C_cons[i] = countCharOccurrences('C', as.character(tr$Consensus[i]))
  tr$fr_A_repeat[i] = countCharOccurrences('A', as.character(tr$RepeatsRegion[i]))
  tr$fr_T_repeat[i] = countCharOccurrences('T', as.character(tr$RepeatsRegion[i]))
  tr$fr_G_repeat[i] = countCharOccurrences('G', as.character(tr$RepeatsRegion[i]))
  tr$fr_C_repeat[i] = countCharOccurrences('C', as.character(tr$RepeatsRegion[i]))
}


### get species specific data

tr$Number = 1

sums = aggregate(tr[, c('FullLength', 'ConsensusLength', 'fr_A_cons', 'fr_T_cons', 'fr_G_cons',
                    'fr_C_cons', 'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', 'fr_C_repeat', 'Number')],
                 by=list('Species'), sum)

agg = aggregate(. ~ Species, sum, data=tr[, c('Species', 'FullLength', 'ConsensusLength', 'fr_A_cons', 'fr_T_cons', 'fr_G_cons',
                                              'fr_C_cons', 'fr_A_repeat', 'fr_T_repeat', 'fr_G_repeat', 'fr_C_repeat', 'Number')])

agg2 = aggregate(. ~ Species, mean, data = tr[, c('Species', 'CopyNumber', 'PercentMatches')])

data = merge(agg, agg2, by='Species')

# create conflict
