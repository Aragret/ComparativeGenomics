##### Q: which one of repeat types correlates better with genome length ?

rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')


cor.test(CHOR$GenomeLength, CHOR$REP.LengthOfTandemRepeats, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.DirRepLength, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.SymmRepLength, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.ComplRepLength, method = 'spearman')
cor.test(CHOR$GenomeLength, CHOR$REP.InvRepLength, method = 'spearman')



a = summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength + CHOR$REP.SymmRepLength +
       CHOR$REP.ComplRepLength + CHOR$REP.InvRepLength)); a
a = summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength + CHOR$REP.SymmRepLength +
                 CHOR$REP.InvRepLength)); a
a = summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength +
                 CHOR$REP.InvRepLength)); a

