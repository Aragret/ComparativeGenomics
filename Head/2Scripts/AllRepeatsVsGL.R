rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')

summary(lm(CHOR$GenomeLength ~ CHOR$REP.LengthOfTandemRepeats + CHOR$REP.DirRepLength +
             CHOR$REP.SymmRepLength + CHOR$REP.ComplRepLength + CHOR$REP.InvRepLength))

#Coefficients:
#  Estimate Std. Error  t value Pr(>|t|)    
#(Intercept)                     1.660e+04  1.344e+01 1235.445  < 2e-16 ***
#  CHOR$REP.LengthOfTandemRepeats  6.419e-01  1.888e-02   33.994  < 2e-16 ***
#  CHOR$REP.DirRepLength           5.673e-02  3.027e-03   18.740  < 2e-16 ***
#  CHOR$REP.SymmRepLength         -2.496e-02  4.354e-03   -5.732 1.06e-08 ***
#  CHOR$REP.ComplRepLength        -1.358e-02  1.498e-02   -0.906    0.365    
#  CHOR$REP.InvRepLength          -6.523e-02  1.122e-02   -5.815 6.54e-09 ***

# Residual standard error: 495.7 on 3948 degrees of freedom
# Multiple R-squared:  0.4848,	Adjusted R-squared:  0.4842 
# F-statistic: 743.1 on 5 and 3948 DF,  p-value: < 2.2e-16


### Repeats coverage 
CHOR$TRCoverage = CHOR$REP.LengthOfTandemRepeats / CHOR$GenomeLength
CHOR$DirRepCoverage = CHOR$REP.DirRepLength / CHOR$GenomeLength
CHOR$SymmRepCoverage = CHOR$REP.SymmRepLength / CHOR$GenomeLength
CHOR$ComplRepCoverage = CHOR$REP.ComplRepLength / CHOR$GenomeLength
CHOR$InvRepCoverage = CHOR$REP.InvRepLength / CHOR$GenomeLength

summary(lm(scale(CHOR$GenomeLength) ~ scale(CHOR$TRCoverage) + scale(CHOR$DirRepCoverage) + scale(CHOR$SymmRepCoverage) +
             scale(CHOR$ComplRepCoverage) + scale(CHOR$InvRepCoverage)))

# scale(CHOR$TRCoverage)        4.994e-01  1.454e-02  34.347  < 2e-16 ***
#  scale(CHOR$DirRepCoverage)    3.874e-01  1.932e-02  20.051  < 2e-16 ***
#  scale(CHOR$SymmRepCoverage)  -1.668e-01  1.750e-02  -9.533  < 2e-16 ***
#  scale(CHOR$ComplRepCoverage) -5.680e-02  2.537e-02  -2.239   0.0252 *  
#  scale(CHOR$InvRepCoverage)   -1.386e-01  2.472e-02  -5.606 2.22e-08 ***


####################################################################################
### Perfect TR

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
tr = read.table('../../Body/2Derived/TRFinder.txt', header=TRUE, sep='\t')

### get perfect repeats and merge them to previous data

perfect_tr = tr[tr$PercentMatches == 100,]
length(unique(perfect_tr$Species)) #785

VEC = unique(perfect_tr$Species); length(VEC)
for (i in 1:length(VEC))
{  # i = 1
  Species = VEC[i];
  TEMP = perfect_tr[perfect_tr$Species == Species,]
  vec_all = c(1); vec_all = vec_all[-1];
  NumberOfTandemRepeats = nrow(TEMP);
  for (j in 1:nrow(TEMP))
  {
    vec = seq(TEMP$Start[j],TEMP$End[j],1);
    vec_all = c(vec_all,vec)
  }
  LengthOfPerfectTandemRepeats = length(vec_all); vec_all = unique(vec_all); LengthOfTandemRepeatsWithoutOverlaps = length(vec_all);
  result_line = data.frame(Species, LengthOfPerfectTandemRepeats)
  if (i == 1) {FINAL = result_line}
  if (i >  1) {FINAL = rbind(FINAL,result_line)}
}

data = merge(CHOR, FINAL, by='Species')

### repeats coverage

data$PerfectTRCoverage = data$LengthOfPerfectTandemRepeats / data$GenomeLength
data$DirRepCoverage = data$REP.DirRepLength / data$GenomeLength
data$SymmRepCoverage = data$REP.SymmRepLength / data$GenomeLength
data$ComplRepCoverage = data$REP.ComplRepLength / data$GenomeLength
data$InvRepCoverage = data$REP.InvRepLength / data$GenomeLength

summary(lm(scale(data$GenomeLength) ~ scale(data$PerfectTRCoverage) + scale(data$DirRepCoverage) + 
             scale(data$SymmRepCoverage) + scale(data$ComplRepCoverage) + scale(data$InvRepCoverage)))


#scale(data$PerfectTRCoverage)  3.098e-01  2.946e-02  10.513  < 2e-16 ***
#  scale(data$DirRepCoverage)     8.740e-01  4.738e-02  18.445  < 2e-16 ***
#  scale(data$SymmRepCoverage)   -2.163e-01  4.467e-02  -4.843 1.54e-06 ***
#  scale(data$ComplRepCoverage)  -6.337e-01  1.126e-01  -5.630 2.52e-08 ***
#  scale(data$InvRepCoverage)     3.785e-02  1.005e-01   0.377    0.707    






####################################################################################
###  PIC - don't look, I haven't finished it yet



library(ape)

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
tree = read.tree('../../Body/1Raw/mtalign.aln.treefile.rooted')

### merge tree and table to use PIC

data = CHOR[which(as.character(CHOR$Species) %in% tree$tip.label),]
df_vec <- as.character(CHOR$Species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

row.names(data) = data$Species

tree2 <- drop.tip(tree, b)

features <- c('GenomeLength', 'REP.LengthOfTandemRepeatsWithoutOverlaps', 
              'REP.DirRepLength', 'REP.SymmRepLength', 'REP.ComplRepLength',
              'REP.InvRepLength')

row.names(data) = data$Species
taxon_vec = unique(data$TAXON)

one_line = c()
for(k in taxon_vec){
  TempData <- data[data$TAXON == k, features]
  for (i in 1:length(colnames(TempData))){
    TempData = TempData[!is.na(TempData[,i]),]
  }
  temp_diff = setdiff(tree2$tip.label, row.names(TempData))
  temp_tree <- drop.tip(tree2, temp_diff)
  
  contrasts <- apply(TempData, 2, pic, temp_tree)
  
  a = summary(lm(contrasts[,1] ~ 0 + contrasts[,2] +  contrasts[,3] +
                   contrasts[,4] + contrasts[,5] + contrasts[,6]))
  one_line = rbind(one_line, c(k, a$coefficients[,1], a$coefficients[,4]))
}

contrast_table <- as.data.frame(one_line)
names(contrast_table) = c('TAXON', 'SlopeTandemLength', 'SlopeDirLength', 'SlopeSymmLength',
                          'SlopeComplLength', 'SlopeInvLength', 'PValueTandLength', 'PValueDirLength', 
                          'PValueSymmLength', 'PValueComplLength', 'PValueInvLength')
