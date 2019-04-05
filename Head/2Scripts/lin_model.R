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

library(pacman)
p_load(tibble, dplyr, magrittr, purrr)
contrasts <- data %>% 
  select(GenerationLength_d, FrA, FrT, FrG, FrC, DRCoverage) %>% 
  map(pic, tree2)

# summary(contrasts$GenerationLength_d)

contrasts = as.data.frame(contrasts)

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

contrasts <- data %>% 
  select(GenerationLength_d, FrA, FrT, FrG, FrC, DRCoverage) %>% 
  mutate_if(is.numeric, log2) %>% 
  map(pic, tree2)

# summary(contrasts$GenerationLength_d)

contrasts = as.data.frame(contrasts)

summary(lm(DRCoverage ~ GenerationLength_d + FrA, contrasts))
# Estimate Std. Error t value Pr(>|t|)   
# (Intercept)        -0.214368   0.095317  -2.249   0.0248 * 
#   GenerationLength_d -0.005679   0.013070  -0.435   0.6641   
# FrA                 0.947124   0.293243   3.230   0.0013 **

summary(lm(DRCoverage ~ 0 + GenerationLength_d + FrT, contrasts))
# GenerationLength_d 0.003499   0.014801   0.236   0.8132  
# FrT                0.294066   0.128694   2.285   0.0226 *

summary(lm(DRCoverage ~ GenerationLength_d + FrG, contrasts))
# Estimate Std. Error t value Pr(>|t|)  
# (Intercept)        -0.20438    0.09608  -2.127   0.0337 *
#   GenerationLength_d -0.01170    0.01351  -0.866   0.3868  
# FrG                -0.08649    0.18640  -0.464   0.6428  

summary(lm(DRCoverage ~ GenerationLength_d + FrC, contrasts))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -0.18158    0.09471  -1.917   0.0556 .  
# GenerationLength_d  0.02153    0.01482   1.452   0.1469    
# FrC                -0.54527    0.11781  -4.628 4.39e-06 ***

summary(lm(DRCoverage ~ GenerationLength_d, contrasts))
# (Intercept)        -0.20643    0.09592  -2.152   0.0317 *
#  GenerationLength_d -0.01349    0.01293  -1.044   0.2970  

