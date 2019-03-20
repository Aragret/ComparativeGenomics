rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

Rep = read.table('../../Body/3Results/RepsCountPseudoStep.csv', sep = ';', header = FALSE)
Rep[is.na(Rep)] <- 0
Rep$TotalLengthOfDirRepeats = Rep$V3*Rep$V4 + Rep$V5*Rep$V6  + Rep$V7*Rep$V8  + Rep$V9*Rep$V10 + Rep$V11*Rep$V12 + Rep$V13*Rep$V14  + Rep$V15*Rep$V16  + Rep$V17*Rep$V18   + Rep$V19*Rep$V20   + Rep$V21*Rep$V22   + Rep$V23*Rep$V24   + Rep$V25*Rep$V26   + Rep$V27*Rep$V28   + Rep$V29*Rep$V30   + Rep$V31*Rep$V32   + Rep$V33*Rep$V34   + Rep$V35*Rep$V36   + Rep$V37*Rep$V38   + Rep$V39*Rep$V40   + Rep$V41*Rep$V42   + Rep$V43*Rep$V44   + Rep$V45*Rep$V46   + Rep$V47*Rep$V48   + Rep$V49*Rep$V50   + Rep$V51*Rep$V52   + Rep$V53*Rep$V54

pdf('../../Body/4Figures/NaiveSimulation.RepeatsAsAFunctionOfNuclFractions.R.01.pdf')
boxplot(Rep$TotalLengthOfDirRepeats[1:10],Rep$TotalLengthOfDirRepeats[11:20],Rep$TotalLengthOfDirRepeats[21:30],Rep$TotalLengthOfDirRepeats[31:40],Rep$TotalLengthOfDirRepeats[41:50],Rep$TotalLengthOfDirRepeats[51:60],Rep$TotalLengthOfDirRepeats[61:70],Rep$TotalLengthOfDirRepeats[71:80],Rep$TotalLengthOfDirRepeats[81:90])
dev.off()