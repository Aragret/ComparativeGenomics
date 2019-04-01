rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

library(ggplot2)
library(dplyr)

Rep = read.table('../../Body/3Results/RepsCountPseudoStep.csv', sep = ';', header = FALSE)
Rep[is.na(Rep)] <- 0
Rep$TotalLengthOfDirRepeats = Rep$V3*Rep$V4 + Rep$V5*Rep$V6  + Rep$V7*Rep$V8  + Rep$V9*Rep$V10 + Rep$V11*Rep$V12 + Rep$V13*Rep$V14  + Rep$V15*Rep$V16  + Rep$V17*Rep$V18   + Rep$V19*Rep$V20   + Rep$V21*Rep$V22   + Rep$V23*Rep$V24   + Rep$V25*Rep$V26   + Rep$V27*Rep$V28   + Rep$V29*Rep$V30   + Rep$V31*Rep$V32   + Rep$V33*Rep$V34   + Rep$V35*Rep$V36   + Rep$V37*Rep$V38   + Rep$V39*Rep$V40   + Rep$V41*Rep$V42   + Rep$V43*Rep$V44   + Rep$V45*Rep$V46   + Rep$V47*Rep$V48   + Rep$V49*Rep$V50   + Rep$V51*Rep$V52   + Rep$V53*Rep$V54
Rep$TotalDirRepCoverage = Rep$TotalLengthOfDirRepeats / 16000

pdf('../../Body/4Figures/NaiveSimulation.RepeatsAsAFunctionOfNuclFractions.R.01.pdf')

getFrA = function(x){
  return(as.numeric(strsplit(as.character(x), ' ')[[1]][6]))
}
Rep$FrA = sapply(Rep$V2, FUN=getFrA)

grRepeats = group_by(Rep, FrA)

# summ = summarise(grRepeats, meanTotalDrCoverage = mean(TotalDirRepCoverage),
#                  y_max = mean(TotalDirRepCoverage) + 
#                    1.96 * sd(TotalDirRepCoverage) / sqrt(length(TotalDirRepCoverage)), 
#                  y_min = mean(TotalDirRepCoverage) - 
#                    1.96 * sd(TotalDirRepCoverage) / sqrt(length(TotalDirRepCoverage)))

summ = summarise(grRepeats, meanTotalDrCoverage = mean(TotalDirRepCoverage),
                 y_max = max(TotalDirRepCoverage), 
                 y_min = min(TotalDirRepCoverage))


ggplot(summ, aes(FrA, meanTotalDrCoverage)) + 
  geom_point(size = 2) + 
  geom_line(color = 'red') +
  geom_errorbar(aes(ymin = y_min, ymax = y_max), width = 0.01) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('Fraction of A') + ylab('Fraction of direct repeats')


boxplot(Rep$TotalLengthOfDirRepeats[1:10],Rep$TotalLengthOfDirRepeats[11:20],Rep$TotalLengthOfDirRepeats[21:30],Rep$TotalLengthOfDirRepeats[31:40],Rep$TotalLengthOfDirRepeats[41:50],Rep$TotalLengthOfDirRepeats[51:60],Rep$TotalLengthOfDirRepeats[61:70],Rep$TotalLengthOfDirRepeats[71:80],Rep$TotalLengthOfDirRepeats[81:90])

# dev.off()

###########################################################################
# new simulations 

Rep = read.table('../../Body/2Derived/RepsCount.csv', sep = ',', header = FALSE)[-1,]

# write.csv(Rep, '../../Body/2Derived/RepsCount.csv', row.names = FALSE, quote = FALSE)
# Rep = read.table('../../Body/2Derived/RepsCount.csv', sep = ';')[, -1]

for(i in 3:ncol(Rep)){
  Rep[, i] = as.integer(as.character(Rep[, i]))
}

Rep[is.na(Rep)] <- 0
Rep$TotalLengthOfDirRepeats = Rep$V3*Rep$V4 + Rep$V5*Rep$V6  + Rep$V7*Rep$V8  + Rep$V9*Rep$V10 + Rep$V11*Rep$V12 + Rep$V13*Rep$V14  + Rep$V15*Rep$V16  + Rep$V17*Rep$V18   + Rep$V19*Rep$V20   + Rep$V21*Rep$V22   + Rep$V23*Rep$V24   + Rep$V25*Rep$V26   + Rep$V27*Rep$V28   + Rep$V29*Rep$V30   + Rep$V31*Rep$V32   + Rep$V33*Rep$V34   + Rep$V35*Rep$V36   + Rep$V37*Rep$V38   + Rep$V39*Rep$V40   + Rep$V41*Rep$V42   + Rep$V43*Rep$V44   + Rep$V45*Rep$V46   + Rep$V47*Rep$V48   + Rep$V49*Rep$V50   + Rep$V51*Rep$V52   + Rep$V53*Rep$V54
Rep$TotalDirRepCoverage = Rep$TotalLengthOfDirRepeats / 16000

getFrA = function(x){
  return(as.numeric(strsplit(as.character(x), ' ')[[1]][6]))
}
Rep$FrA = sapply(Rep$V2, FUN=getFrA)

grRepeats = group_by(Rep, FrA)

summ = summarise(grRepeats, meanTotalDrCoverage = mean(TotalDirRepCoverage),
                 y_max = max(TotalDirRepCoverage), 
                 y_min = min(TotalDirRepCoverage))

ggplot(summ, aes(FrA, meanTotalDrCoverage)) + 
  geom_point(size = 2) + 
  geom_line(color = 'red') +
  geom_errorbar(aes(ymin = y_min, ymax = y_max), width = 0.01) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlab('Fraction of A') + ylab('Fraction of direct repeats')


boxplot(Rep$TotalLengthOfDirRepeats[1:100],Rep$TotalLengthOfDirRepeats[101:200],Rep$TotalLengthOfDirRepeats[201:300],Rep$TotalLengthOfDirRepeats[301:400],Rep$TotalLengthOfDirRepeats[401:500],Rep$TotalLengthOfDirRepeats[501:600],Rep$TotalLengthOfDirRepeats[601:700],Rep$TotalLengthOfDirRepeats[701:800],Rep$TotalLengthOfDirRepeats[801:900])

dev.off()

