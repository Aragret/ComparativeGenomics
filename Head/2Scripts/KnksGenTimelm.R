CHOR = read.table('../../Body/2Derived/MitGenomics.txt', sep = '\t', header = TRUE)
kaks = read.table('../../Body/1Raw/KaKsRL_data.txt', header = TRUE, sep='\t')

data = merge(kaks, CHOR[, c('Species', 'REP.NumberOfTandemRepeats')], by.x = 'SPECIES',
             by.y = 'Species')

summary(lm(data$lnDN_DS ~ data$lnBM_g + data$REP.NumberOfTandemRepeats))
