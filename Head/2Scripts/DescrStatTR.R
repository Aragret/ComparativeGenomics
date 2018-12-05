rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 
library(gdata) # install.packages("gdata")

tr = read.table('../../Body/2Derived/TRFinder.txt', header=TRUE, sep='\t')
GenLength <- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header = TRUE, sep='\t')

GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

tr$FullLength = tr$End - tr$Start
tr$ConsensusLength = as.numeric(lapply(as.character(tr$Consensus), nchar))


#### Perfect, non perfect TRs

hist(tr$PercentMatches, breaks = 50, main = 'Percent of matches', xlab = 'Percent')
plot(tr$PercentMatches, tr$FullLength)
plot(tr$PercentMatches, tr$ConsensusLength)
plot(tr$PercentMatches, tr$CopyNumber)

### With generation length
# percent matches average ?

tr$number = 1; 

All = aggregate(list(tr$number,tr$FullLength), by = list(tr$Species), FUN = sum)
names(All) = c('Species','NumberOfAllTR','TotalLengthOfAllTR')
All = merge(All,GenLength)
All = All[All$TotalLengthOfAllTR > 0,]


for (threshold in 1:9)
{ # threshold = 1
  PerfectTr = tr[tr$PercentMatches >= quantile(tr$PercentMatches, threshold/10),]; 
  PerfectAGG = aggregate(list(PerfectTr$number, PerfectTr$FullLength), by = list(PerfectTr$Species), FUN = sum)
  names(PerfectAGG) = c('Species','NumberOfPerfectTR','TotalLengthOfPerfectTR')
  Perfect = merge(PerfectAGG,GenLength)
  Perfect = Perfect[Perfect$TotalLengthOfLongTR > 0,]
  title = paste('Perfectness of TR >= ', threshold/10, ' percentile')
  plot(Perfect$GenerationLength_d, Perfect$TotalLengthOfLongTR, main = title, ylim = c(0,2000), xlim = c(0,10500)) # triangle => maybe fisher test?
}

