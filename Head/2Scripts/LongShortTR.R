### correlation of TR length and GT in mammals

rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

pdf("../../Body/4Figures/LongShortTR.R.o1.pdf")

library(gdata) # install.packages("gdata")

GenLength <- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header = TRUE, sep='\t')
tr = read.table('../../Body/2Derived/TRFinder.txt', sep='\t', header = TRUE)

GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

### get average for species with several lines in GL!!
tr = tr[tr$Species %in% GenLength$Species,]
tr$FullLength = tr$End - tr$Start
tr$ConsensusLength = as.numeric(lapply(as.character(tr$Consensus), nchar))

par(mfrow=c(2,2))
hist(tr$ConsensusLength, breaks = 50, main = 'majority of TRs has short consensus')
hist(tr$CopyNumber, breaks = 50, main = 'majority of TRs are low-copied')
plot(tr$ConsensusLength,tr$CopyNumber, main = paste('TRs with long consensus are less-copied','than TRs with short consensus', sep = '\n'))
plot(tr$ConsensusLength,tr$FullLength) # there is a separate spot!!! who they are?

par(mfrow=c(1,1))

hist(tr$FullLength, breaks = 50)
tr$number = 1; 

All = aggregate(list(tr$number,tr$FullLength), by = list(tr$Species), FUN = sum)
names(All) = c('Species','NumberOfAllTR','TotalLengthOfAllTR')
All = merge(All,GenLength)
All = All[All$TotalLengthOfAllTR > 0,]
plot(All$GenerationLength_d,All$TotalLengthOfAllTR,xlim = c(0,10500), ylim = c(0,2000)) # triangle => maybe fisher test?

for (threshold in 1:9)
{ # threshold = 1
LongTr = tr[tr$ConsensusLength >= quantile(tr$ConsensusLength,threshold/10),]; 
LongAGG = aggregate(list(LongTr$number,LongTr$FullLength), by = list(LongTr$Species), FUN = sum)
names(LongAGG) = c('Species','NumberOfLongTR','TotalLengthOfLongTR')
Long = merge(LongAGG,GenLength)
Long = Long[Long$TotalLengthOfLongTR > 0,]
title = paste('length of TR >= ',threshold/10, ' percentile')
plot(Long$GenerationLength_d,Long$TotalLengthOfLongTR, main = title, ylim = c(0,2000), xlim = c(0,10500)) # triangle => maybe fisher test?
}

dev.off()
