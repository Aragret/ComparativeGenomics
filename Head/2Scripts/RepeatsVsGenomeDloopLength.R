##### Q: which one of repeat types correlates better with genome length ?

rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
CHOR = CHOR[CHOR$TAXON != 'AncientFish',]


cor.test(CHOR$GenomeLength, CHOR$REP.LengthOfTandemRepeats, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.DirRepLength, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.SymmRepLength, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.ComplRepLength, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.InvRepLength, method = 'spearman')

taxon = unique(CHOR$TAXON)

for (i in taxon) {
  tempData = CHOR[CHOR$TAXON == i,]
  print(i)
  print(c(cor.test(tempData$GenomeLength, tempData$REP.LengthOfTandemRepeats, method = 'spearman'),
        cor.test(tempData$GenomeLength, tempData$REP.DirRepLength, method = 'spearman'),
        cor.test(tempData$GenomeLength, tempData$REP.SymmRepLength, method = 'spearman'),
        cor.test(tempData$GenomeLength, tempData$REP.ComplRepLength, method = 'spearman'),
        cor.test(tempData$GenomeLength, tempData$REP.InvRepLength, method = 'spearman')))
}

a = summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength + CHOR$REP.SymmRepLength +
       CHOR$REP.ComplRepLength + CHOR$REP.InvRepLength)); a
a = summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength + CHOR$REP.SymmRepLength +
                 CHOR$REP.InvRepLength)); a
a = summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength +
                 CHOR$REP.InvRepLength)); a

