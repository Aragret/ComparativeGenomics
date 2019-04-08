rm(list = ls())
library(gdata)

# pairwise correlations

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
GenLength = read.xls('../../Body/1Raw/GenerationLengthForMammals.xlsx')

GenLength$Species = gsub(' ', '_', GenLength$Scientific_name)
GenLength = GenLength[, c(14,16)]

data = merge(CHOR, GenLength, by='Species')

summary(data$REP.DirRepLength)

data = data[, c('Species', 'A', 'T', 'G', 'C', 'REP.DirRepLength', 'GenomeLength', 'GenerationLength_d')]

data$FrA = data$A / data$GenomeLength
data$FrT = data$T / data$GenomeLength
data$FrG = data$G / data$GenomeLength
data$FrC = data$C / data$GenomeLength

data$DRCoverage = data$REP.DirRepLength / data$GenomeLength

summary(lm(log(data$GenerationLength_d) ~ data$DRCoverage + data$FrA + data$FrT + data$FrG))

# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      17.9367     1.4528  12.346  < 2e-16 ***
#   data$DRCoverage  -0.0174     0.3675  -0.047    0.962    
#   data$FrA        -13.3814     2.9379  -4.555 6.18e-06 ***
#   data$FrT        -17.9156     1.4846 -12.067  < 2e-16 ***
#   data$FrG         -6.7906     4.4961  -1.510    0.131    

summary(lm(log(data$DRCoverage) ~ log(data$GenerationLength_d) + data$FrA + data$FrT + data$FrG))
# (Intercept)                  -3.845506   0.545880  -7.045 4.45e-12 ***
#   log(data$GenerationLength_d) -0.006054   0.012950  -0.468 0.640277    
# data$FrA                      3.357545   1.014683   3.309 0.000984 ***
#   data$FrT                      2.713117   0.550534   4.928 1.04e-06 ***
#   data$FrG                      4.175280   1.531670   2.726 0.006571 ** 

summary(lm(log(data$DRCoverage) ~ log(data$GenerationLength_d)))

################################################################################

cor.test(data$GenerationLength_d, data$DRCoverage, method = 'spearman')
# rho -0.1519978, pvalue 5.011e-05

summary(lm(data$DRCoverage ~ data$FrA + data$FrT + data$FrG))

# (Intercept)  -0.5108     0.1480  -3.452 0.000590 ***
#   data$FrA      0.9762     0.2995   3.260 0.001170 ** 
#   data$FrT      0.7906     0.1495   5.287 1.66e-07 ***
#   data$FrG      1.5478     0.4581   3.379 0.000768 ***

par(mfrow=c(2, 2))
hist(data$FrA, breaks = 50)
hist(data$FrT, breaks = 50)
hist(data$FrG, breaks = 50)
hist(data$FrC, breaks = 50)

summary(data$FrA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2930  0.3227  0.3302  0.3295  0.3359  0.3603 
summary(data$FrT)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2274  0.2640  0.2739  0.2764  0.2890  0.3388 
summary(data$FrG)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1136  0.1282  0.1318  0.1327  0.1363  0.1582 
summary(data$FrC)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1917  0.2500  0.2613  0.2614  0.2704  0.3275 
#####################################################################################
### pairwise correlations

nucl_freqs = data[, c('FrA', 'FrT', 'FrG', 'FrC')]

round(cor(nucl_freqs), 2)

pdf('../../Body/4Figures/lin_model.R.NucleotidesPairwiseCorrs.pdf')
plot(nucl_freqs)
dev.off()

#####################################################################################
############################## PICs

library(ape) # install.packages('ape') 

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

data = data[which(as.character(data$Species) %in% tree$tip.label),]
data = data[-439,]
row.names(data) = data$Species

df_vec <- as.character(data$Species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

tree2 <- drop.tip(tree, b)

#################################################################################
######## without log

library(pacman)
p_load(tibble, dplyr, magrittr, purrr)
contrasts <- data %>% 
  select(GenerationLength_d, FrA, FrT, FrG, FrC, DRCoverage) %>% 
  map(pic, tree2)

# summary(contrasts$GenerationLength_d)

contrasts = as.data.frame(contrasts)

cor.test(contrasts$GenerationLength_d, contrasts$DRCoverage, method = 'spearman')
plot(contrasts$GenerationLength_d, contrasts$DRCoverage)

summary(lm(DRCoverage ~ GenerationLength_d + FrA, contrasts))
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)        -3.487e-02  1.982e-02  -1.759  0.07895 . 
# GenerationLength_d  3.275e-06  1.501e-06   2.183  0.02938 * 
#   FrA                 7.427e-01  2.678e-01   2.774  0.00569 **

summary(lm(DRCoverage ~ GenerationLength_d + FrT, contrasts))
# (Intercept)        -2.794e-02  1.984e-02  -1.408 0.159519    
# GenerationLength_d  4.663e-06  1.610e-06   2.896 0.003892 ** 
#   FrT                 4.356e-01  1.316e-01   3.311 0.000977 ***

summary(lm(DRCoverage ~ GenerationLength_d + FrG, contrasts))
# (Intercept)        -3.369e-02  1.995e-02  -1.689   0.0917 .
# GenerationLength_d  2.467e-06  1.520e-06   1.623   0.1051  
# FrG                 2.940e-02  4.112e-01   0.071   0.9430  

summary(lm(DRCoverage ~ GenerationLength_d + FrC, contrasts))
# (Intercept)        -2.789e-02  1.964e-02  -1.420 0.156091    
# GenerationLength_d  5.787e-06  1.612e-06   3.590 0.000354 ***
#   FrC                -6.324e-01  1.319e-01  -4.795 1.99e-06 ***

summary(lm(DRCoverage ~ GenerationLength_d, contrasts))
# (Intercept)        -3.362e-02  1.991e-02  -1.689   0.0917 .
# GenerationLength_d  2.491e-06  1.481e-06   1.683   0.0929 .


summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrA, contrasts))
# Estimate Std. Error t value Pr(>|t|)   
# GenerationLength_d 3.350e-06  1.502e-06    2.23  0.02604 * 
#   FrA                7.319e-01  2.681e-01    2.73  0.00649 **

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrT, contrasts))
# GenerationLength_d 4.811e-06  1.608e-06   2.993 0.002860 ** 
#  FrT                4.516e-01  1.312e-01   3.443 0.000609 ***

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrG, contrasts))
# GenerationLength_d  2.577e-06  1.521e-06   1.694   0.0906 .
# FrG                -3.231e-03  4.113e-01  -0.008   0.9937 

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrC, contrasts))
# GenerationLength_d  5.915e-06  1.611e-06   3.673 0.000258 ***
# FrC                -6.438e-01  1.317e-01  -4.887 1.27e-06 ***

summary(lm(DRCoverage ~ 0 + GenerationLength_d, contrasts))
# GenerationLength_d 2.575e-06  1.482e-06   1.737   0.0827 .

#################################################################################
############### with log

contrasts <- data %>% 
  select(GenerationLength_d, FrA, FrT, FrG, FrC, DRCoverage) %>% 
  mutate_if(is.numeric, log2) %>% 
  map(pic, tree2)

# summary(contrasts$GenerationLength_d)

contrasts = as.data.frame(contrasts)

cor.test(contrasts$GenerationLength_d, contrasts$DRCoverage, method = 'spearman')

summary(lm(DRCoverage ~ GenerationLength_d + FrA, contrasts))
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)        -0.214368   0.095317  -2.249   0.0248 * 
#   GenerationLength_d -0.005679   0.013070  -0.435   0.6641   
# FrA                 0.947124   0.293243   3.230   0.0013 **

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrA, contrasts))
# GenerationLength_d -0.005294   0.013107  -0.404  0.68642   
# FrA                 0.930128   0.293992   3.164  0.00162 **

summary(lm(DRCoverage ~ GenerationLength_d + FrT, contrasts))
# (Intercept)        -0.189388   0.096025  -1.972   0.0490 *
# GenerationLength_d  0.001839   0.014794   0.124   0.9011  
# FrT                 0.272742   0.128884   2.116   0.0347 *

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrT, contrasts))
# GenerationLength_d 0.003499   0.014801   0.236   0.8132  
# FrT                0.294066   0.128694   2.285   0.0226 *

summary(lm(DRCoverage ~ GenerationLength_d + FrG, contrasts))
# Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        -0.20438    0.09608  -2.127   0.0337 *
#   GenerationLength_d -0.01170    0.01351  -0.866   0.3868  
# FrG                -0.08649    0.18640  -0.464   0.6428  

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrG, contrasts))
# GenerationLength_d -0.01082    0.01353  -0.799    0.424
# FrG                -0.10474    0.18667  -0.561    0.575

summary(lm(DRCoverage ~ GenerationLength_d + FrC, contrasts))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -0.18158    0.09471  -1.917   0.0556 .  
# GenerationLength_d  0.02153    0.01482   1.452   0.1469    
# FrC                -0.54527    0.11781  -4.628 4.39e-06 ***

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrC, contrasts))
# GenerationLength_d  0.02280    0.01484   1.536    0.125    
# FrC                -0.55808    0.11784  -4.736 2.64e-06 ***

summary(lm(DRCoverage ~ GenerationLength_d, contrasts))
# (Intercept)        -0.20643    0.09592  -2.152   0.0317 *
#  GenerationLength_d -0.01349    0.01293  -1.044   0.2970  

summary(lm(DRCoverage ~ 0 + GenerationLength_d, contrasts))
# GenerationLength_d -0.01299    0.01296  -1.002    0.317
