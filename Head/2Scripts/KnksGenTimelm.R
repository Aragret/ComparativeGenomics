library(mlr)

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

summary(data$StatusRL)

for(i in 1:nrow(data)){
  # i = 1
  if(as.character(data[i, 'StatusRL']) == 'LC'){
    data[i, 'LeastConcern'] = 1
  }
  else{data[i, 'LeastConcern'] = 0}
}

data$LeastConcern = as.factor(data$LeastConcern)
summary(data$LeastConcern)

summary(lm(data$lnDN_DS ~ data$LeastConcern + data$lnBM_g + data$REP.NumberOfTandemRepeats))

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    -3.526970   0.120534 -29.261  < 2e-16 ***
#   data$LeastConcern1             -0.111004   0.076322  -1.454  0.14742    
# data$lnBM_g                     0.041366   0.009979   4.145 5.04e-05 ***
#   data$REP.NumberOfTandemRepeats -0.059023   0.022310  -2.646  0.00881 ** 


a = createDummyFeatures(data, cols='StatusRL')

summary(lm(a$lnDN_DS ~ a$lnBM_g + a$REP.NumberOfTandemRepeats + a$StatusRL.CR + a$StatusRL.EN +
             a$StatusRL.LC + a$StatusRL.NT + a$StatusRL.VU))

# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                 -3.57081    0.14616 -24.431  < 2e-16 ***
#   a$lnBM_g                     0.03943    0.01027   3.839 0.000167 ***
#   a$REP.NumberOfTandemRepeats -0.05865    0.02255  -2.601 0.010007 *  
#   a$StatusRL.CR                0.15062    0.20593   0.731 0.465401    
# a$StatusRL.EN                0.13619    0.13787   0.988 0.324493    
# a$StatusRL.LC               -0.05215    0.10841  -0.481 0.631036    
# a$StatusRL.NT                0.02332    0.14861   0.157 0.875491    
# a$StatusRL.VU                     NA         NA      NA       NA    

