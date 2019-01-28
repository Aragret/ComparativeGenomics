rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

library(dplyr) # install.packages("dplyr")

tr = read.table('../../Body/2Derived/TRFinder.txt', header=TRUE, sep='\t')
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header = TRUE, sep='\t')
neutr_nuc = read.table('../../Body/2Derived/AllGenesCodonUsageNoOverlap.txt', header=TRUE, sep='\t')
dloops_tr = read.table('../../Body/2Derived/TRinDloops.txt', header=TRUE, sep='\t')

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

agg_neutr = aggregate(. ~ Species, sum, data=neutr_nuc[, c('Species', 'NeutralA', 'NeutralT',
                                                           'NeutralG', 'NeutralC')])
data = merge(tr, agg_neutr, by='Species', all.x=TRUE)

data = merge(data, CHOR[, c('Species', 'GenomeLength', 'A', 'T', 'G', 'C', 'TAXON', 'taxonomy')],
             by = 'Species')

####### in CR or not
dloops_tr$Species = sub('_dloop', '', dloops_tr$Species)
dloops_tr$Species = sub('_control_region', '', dloops_tr$Species)


temp_data = anti_join(tr, dloops_tr[, -c(2, 3)])

nrow(anti_join(temp_data, dloops_tr[, -c(2, 3)]))
nrow(anti_join(dloops_tr[, -c(2, 3)], temp_data))

nrow(anti_join(temp_data, tr))
nrow(anti_join(tr, temp_data))

temp_data$InDloop = 0

final = merge(tr, temp_data, all.x = TRUE)

final[is.na(final$InDloop),]$InDloop = 1

final = merge(final, CHOR[, c('Species', 'GenomeLength', 'A', 'T', 'G', 'C', 'TAXON', 'taxonomy')],
              by = 'Species')

write.table(final, '../../Body/3Results/TandRepInfo.txt', sep='\t', quote=FALSE, row.names=FALSE)
