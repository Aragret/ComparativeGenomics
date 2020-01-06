CHOR = read.table('../../Body/2Derived/MitGenomics.txt', sep = '\t', header = TRUE)
kaks = read.table('../../Body/1Raw/KaKsRL_data.txt', header = TRUE, sep='\t')

data = merge(kaks, CHOR[, c('Species', 'REP.NumberOfTandemRepeats')], by.x = 'SPECIES',
             by.y = 'Species')

summary(lm(data$lnDN_DS ~ data$lnBM_g + data$REP.NumberOfTandemRepeats))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    -3.623693   0.100810 -35.946  < 2e-16 ***
#   data$lnBM_g                     0.045442   0.009604   4.732 4.23e-06 ***
#   data$REP.NumberOfTandemRepeats -0.061748   0.022294  -2.770  0.00614 ** 
