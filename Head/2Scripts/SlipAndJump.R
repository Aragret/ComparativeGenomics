rm(list=ls(all=TRUE))

library(raster) # install.packages("raster")
plots_dir <- normalizePath(file.path("../../Body/4Figures/"))

pdf("../../Body/4Figures/SlipAndJump.R.01.pdf")
  
##### 1: READ microhomology from pair-wise alignments
homol = read.table("../../Body/2Derived/HeatMaps/100x100.csv", sep = ';', header = TRUE)
row.names(homol)=homol$X; homol = homol[,-1]
# make long vertical table from the matrix
for (i in 1:nrow(homol))
{
  for (j in 1:ncol(homol))
  { # i  = 2; j = 1
    FirstWindow = as.character(row.names(homol)[i])
    SecondWindow = as.character(names(homol)[j])
    Score = as.numeric(homol[i,j])
    OneLine = data.frame(FirstWindow,SecondWindow,Score)
    if (i == 1 & j == 1) {Final = OneLine}
    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
MicroHomology = Final

##### 2: READ density of direct repeats per window

DirectRepDensity = read.table("../../Body/2Derived/HeatMaps/Link_matrix_direct_major_activ_left.csv", sep = ';', header = TRUE)
DirectRepDensity = DirectRepDensity[,-1]

# make long vertical table from the matrix
for (i in 1:nrow(DirectRepDensity))
{
  for (j in 1:ncol(DirectRepDensity))
  { # i  = 2; j = 1
    FirstWindow = as.character(row.names(DirectRepDensity)[i])
    SecondWindow = as.character(names(DirectRepDensity)[j])
    Score = as.numeric(DirectRepDensity[i,j])
    OneLine = data.frame(FirstWindow,SecondWindow,Score)
    if (i == 1 & j == 1) {Final = OneLine}
    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
DirectRepDensity = Final

##### 3: correlate MicroHomology and  DirectRepDensity, derive HomologyAndRepeats dataset
DirectRepDensity = DirectRepDensity[order(DirectRepDensity$FirstWindow,DirectRepDensity$SecondWindow),]
MicroHomology = MicroHomology[order(MicroHomology$FirstWindow,MicroHomology$SecondWindow),]
cor.test(DirectRepDensity$Score,MicroHomology$Score, method = 'spearman') # p-value = 1.698e-06, spearman rho = 0.06796994 
nrow(DirectRepDensity) # 4950
plot(DirectRepDensity$Score,MicroHomology$Score)

cor.test(DirectRepDensity[DirectRepDensity$Score > 0,]$Score,MicroHomology[DirectRepDensity$Score > 0,]$Score, method = 'spearman')
plot(DirectRepDensity[DirectRepDensity$Score > 0,]$Score,MicroHomology[DirectRepDensity$Score > 0,]$Score)

HomologyAndRepeats = data.frame(DirectRepDensity$FirstWindow,DirectRepDensity$SecondWindow,DirectRepDensity$Score,MicroHomology$Score)
names(HomologyAndRepeats) = c('FirstWindow','SecondWindow','DirectRepeatsDensity','MicroHomologyScore')

##### 4: READ MITOBREAK AND FILTER (KEEP ONLY MAJOR ARC DELETIONS):

breaks = read.table("../../Body/1Raw/MitoBreakDB_12122019.csv", sep = ',', header = TRUE)
breaks$X5..breakpoint = as.numeric(as.character(breaks$X5..breakpoint)); summary(breaks$X5..breakpoint)
breaks$X3..breakpoint = as.numeric(as.character(breaks$X3..breakpoint)); summary(breaks$X3..breakpoint)
breaks = breaks[!is.na(breaks$X3..breakpoint) & !is.na(breaks$X5..breakpoint),]

breaks$FirstWindowBreakpoint = breaks$X3..breakpoint
breaks$SecondWindowBreakpoint = breaks$X5..breakpoint
breaks = breaks[breaks$FirstWindowBreakpoint > 5781 & breaks$FirstWindowBreakpoint < 16569 & breaks$SecondWindowBreakpoint > 5781 & breaks$SecondWindowBreakpoint < 16569,] # can make it better!! take in ot account 0-100?
# поскольку координаты не такие простые (см ниже) - чтобы не париться можно взять все точки разрыва что больше чем 5781 и меньше чем 16569
# OH: 110-441
# OL: 5721-5781

HomologyAndRepeats$Deletion = 0
for (i in 1:nrow(HomologyAndRepeats))
{ # i = 1
  FirstWindow  = HomologyAndRepeats$FirstWindow[i]
  SecondWindow = HomologyAndRepeats$SecondWindow[i]
  TempBreaks = breaks[breaks$FirstWindowBreakpoint >= FirstWindow & breaks$FirstWindowBreakpoint < (FirstWindow + 100) & breaks$SecondWindowBreakpoint >= SecondWindow & breaks$SecondWindowBreakpoint < (SecondWindow + 100),]
  if (nrow(TempBreaks) > 0) {HomologyAndRepeats$Deletion[i] = 1}
}
  
table(HomologyAndRepeats$Deletion)
# 0    1 
# 4466  484 

a<- glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$DirectRepeatsDensity + HomologyAndRepeats$MicroHomologyScore, family = binomial) 
summary(a)

a<- glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore), family = binomial) 
summary(a) 
#  Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  -2.25264    0.04914 -45.838  < 2e-16 ***
#  scale(HomologyAndRepeats$MicroHomologyScore)  0.27442    0.04697   5.843 5.13e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# (Dispersion parameter for binomial family taken to be 1)
#
# Null deviance: 3169.7  on 4949  degrees of freedom
# Residual deviance: 3135.9  on 4948  degrees of freedom
# AIC: 3139.9
# 
# Number of Fisher Scoring iterations: 5
nrow(HomologyAndRepeats)

### may be add perfect repeats of Orlov - yes or no for each deletion? and see, if microhomology still important?


##### 5: READ GLOBAL FOLDING:

GlobalFolding = read.table("../../Body/2Derived/HeatMaps/Link_matrix1000_major.csv", sep = ';', header = TRUE)
row.names(GlobalFolding) = GlobalFolding$X 
GlobalFolding = GlobalFolding[,-1]

# make long vertical table from the matrix
for (i in 1:nrow(GlobalFolding))
{
  for (j in 1:ncol(GlobalFolding))
  { # i  = 2; j = 1
    FirstWindow = as.character(row.names(GlobalFolding)[i])
    SecondWindow = as.character(names(GlobalFolding)[j])
    Score = as.numeric(GlobalFolding[i,j])
    OneLine = data.frame(FirstWindow,SecondWindow,Score)
    if (i == 1 & j == 1) {Final = OneLine}
    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is noizy and bold)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
GlobalFolding1000 = Final
GlobalFolding1000 = GlobalFolding1000[order(GlobalFolding1000$FirstWindow,GlobalFolding1000$SecondWindow),]
names(GlobalFolding1000) = c('FirstWindowWholeKbRes','SecondWindowWholeKbRes','GlobalFolding1000Score')

##### 5.1: READ GLOBAL FOLDING WITH WINDOW = 100 bp: (it automaticaly rewrites the GlobalFolding matrix from the previous point 5)

GlobalFolding = read.table("../../Body/2Derived/HeatMaps/Link_matrix100hydra_major.csv", sep = ';', header = TRUE)
row.names(GlobalFolding) = GlobalFolding$X 
GlobalFolding = GlobalFolding[,-1]

# make long vertical table from the matrix
for (i in 1:nrow(GlobalFolding))
{
  for (j in 1:ncol(GlobalFolding))
  { # i  = 2; j = 1
    FirstWindow = as.character(row.names(GlobalFolding)[i])
    SecondWindow = as.character(names(GlobalFolding)[j])
    Score = as.numeric(GlobalFolding[i,j])
    OneLine = data.frame(FirstWindow,SecondWindow,Score)
    if (i == 1 & j == 1) {Final = OneLine}
    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final) 
# Should we delete bold diagonal or erase it to zeroes??? If delete, dimension will be decreased - try this. delete 5 windows next to diagonal (500)
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow + 1000,]; nrow(Final) # 500 or 1000!!!!!! similarly good results but 1000 is a bit better
GlobalFolding = Final
GlobalFolding = GlobalFolding[order(GlobalFolding$FirstWindow,GlobalFolding$SecondWindow),]
names(GlobalFolding)[3] = c('GlobalFoldingScore');

# GlobalFolding - is the whole genome without bold diagonal, not only the major arc!! Keep only major arc in downstream analyses.
# will do it when merge with InvRepDens.

###### 6: READ INVERTED REPEATS WITH STEP 1000

InvRepDens = read.table("../../Body/2Derived/HeatMaps/Link_matrix_1000_invert_major_activ_left.csv", sep = ';', header = TRUE)
row.names(InvRepDens) = InvRepDens$X 
InvRepDens = InvRepDens[,-1]

# make long vertical table from the matrix
for (i in 1:nrow(InvRepDens))
{
  for (j in 1:ncol(InvRepDens))
  { # i  = 2; j = 1
    FirstWindow = as.character(row.names(InvRepDens)[i])
    SecondWindow = as.character(names(InvRepDens)[j])
    Score = as.numeric(InvRepDens[i,j])
    OneLine = data.frame(FirstWindow,SecondWindow,Score)
    if (i == 1 & j == 1) {Final = OneLine}
    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
InvRepDens = Final
InvRepDens = InvRepDens[order(InvRepDens$FirstWindow,InvRepDens$SecondWindow),]

###### 6.1: READ INVERTED REPEATS WITH STEP 100 bp (it automatically rewrites InvRepDens from previous point 6)

InvRepDens = read.table("../../Body/2Derived/HeatMaps/Link_matrix_invert_major_activ_left.modified.csv", sep = '\t', header = TRUE, row.names = 1) # , row.names = NULL)

# make long vertical table from the matrix
for (i in 1:nrow(InvRepDens))
{
  for (j in 1:ncol(InvRepDens))
  { # i  = 2; j = 1
    FirstWindow = as.character(row.names(InvRepDens)[i])
    SecondWindow = as.character(names(InvRepDens)[j])
    Score = as.numeric(InvRepDens[i,j])
    OneLine = data.frame(FirstWindow,SecondWindow,Score)
    if (i == 1 & j == 1) {Final = OneLine}
    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
  }
}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
InvRepDens = Final
InvRepDens = InvRepDens[order(InvRepDens$FirstWindow,InvRepDens$SecondWindow),]
summary(InvRepDens$FirstWindow)  # 6000 15800
summary(InvRepDens$SecondWindow) # 5900 15700
names(InvRepDens)[3] = c('InvRepDensScore');

###### 6.2: READ INVERTED REPEATS WITH STEP 100 bp WITH OVERLAPS (it automatically rewrites InvRepDens from previous point 6.1)
# not nice results - use 6.1 
#InvRepDens = read.table("../../Body/2Derived/HeatMaps/Link_matrix_direct_major_activ_left_cross.KpModifByHand.mtrx", sep = '\t',header = TRUE, row.names = 1) # , row.names = NULL)

# make long vertical table from the matrix
#for (i in 1:nrow(InvRepDens))
#{
#  for (j in 1:ncol(InvRepDens))
#  { # i  = 2; j = 1
#    FirstWindow = as.character(row.names(InvRepDens)[i])
#    SecondWindow = as.character(names(InvRepDens)[j])
#    Score = as.numeric(InvRepDens[i,j])
#    OneLine = data.frame(FirstWindow,SecondWindow,Score)
#    if (i == 1 & j == 1) {Final = OneLine}
#    if (i > 1 | j > 1) {Final = rbind(Final,OneLine)}
#  }
#}

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal)
#Final$SecondWindow = gsub('X','',Final$SecondWindow)
#Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
#nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
#InvRepDens = Final
#InvRepDens = InvRepDens[order(InvRepDens$FirstWindow,InvRepDens$SecondWindow),]
#summary(InvRepDens$FirstWindow)  # 6000 15800
#summary(InvRepDens$SecondWindow) # 5900 15700
#names(InvRepDens)[3] = c('InvRepDensScore');

###### 7: CORRELATE GlobalFolding$Score and InvRepDens$Score - weak positive!
merged = merge(InvRepDens,GlobalFolding, by = c("FirstWindow","SecondWindow"))
summary(merged$FirstWindow)  # diag 500: 6500 15800; diag 1000: 7000 - 15800
summary(merged$SecondWindow) # diag 500: 5900 15200; diag 1000: 5900 14700
cor.test(merged$InvRepDensScore,merged$GlobalFoldingScore, method = 'spearman') # diag 500: rho = 0.04926082, p-value = 0.0009922; diag 1000: rho = 0.04945796, p = 0.001743
nrow(merged) # 4005

###### 8: ADD InfinitySign parameter into HomologyAndRepeats dataset (13 - 16 kb vs 6-9 kb):
HomologyAndRepeats$InfinitySign = 0
for (i in 1:nrow(HomologyAndRepeats))
{
  if (HomologyAndRepeats$FirstWindow[i] >= 13000 & HomologyAndRepeats$FirstWindow[i] <= 16000 & HomologyAndRepeats$SecondWindow[i] >= 6000 & HomologyAndRepeats$SecondWindow[i] <= 9000)
  {HomologyAndRepeats$InfinitySign[i] = 1}    
}
table(HomologyAndRepeats$InfinitySign)
summary(HomologyAndRepeats$FirstWindow)  # 6000 15800
summary(HomologyAndRepeats$SecondWindow) # 5900 15700

## merge HomologyAndRepeats with merged(InvRepDens + GlobalFolding)
dim(HomologyAndRepeats) # 4950
HomologyAndRepeats = merge(HomologyAndRepeats,merged, by = c("FirstWindow","SecondWindow"))
dim(HomologyAndRepeats) # diag 500: 4465; diag 1000:  4005

# is GlobalFoldingScore higher within the cross according to our InfinitySign model? YES!!! 
wilcox.test(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore,HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore)# diag 1000: 3.358e-09
boxplot(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore,HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore, notch = TRUE, names = c('stem','loop'), ylab = 'in silico folding score', outline = FALSE)
t.test(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore,HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore) # diag 1000: 0.002639
summary(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore) # diag 1000: 0.3893
summary(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore) # diag 1000: 0.09755

# we have to link better global folding and InfinitySign model - till now it was done by eye. Clusterisation? One cluster? 
# dev.off()
###### 9: LOGISTIC REGRESSION: HomologyAndRepeats$Deletion as a function of HomologyAndRepeats$MicroHomologyScore and HomologyAndRepeats$InfinitySign:

a<-glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore + HomologyAndRepeats$InfinitySign, family = 'binomial')
summary(a)

a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$InfinitySign), family = 'binomial')
summary(a) # PAPER!!! 0.33 + 0.91

a<-glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore + HomologyAndRepeats$GlobalFoldingScore, family = 'binomial')
summary(a) # non significant - may be I have to take it on bigger scale! (1kb without diagonal, because this is global parameter not precise)
# to reconstruct 100 bp matrix back from 1 kb matrix!!!!! 

# get residuals and correlate them with global matrix

a<-glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore, family = 'binomial')
HomologyAndRepeats$Residuals = residuals(a)
summary(HomologyAndRepeats$Residuals)

cor.test(HomologyAndRepeats$Residuals,HomologyAndRepeats$GlobalFoldingScore, method = 'spearman') # diag 500: rho = 0.04714374, p = 0.001627; diag 1000: rho = 0.05628272, p = 0.0003658

#### reconstruct Global folding 100 bp back from 1kb resolution (GlobalFolding1000) assuming that global folding can work remotely enough.
#### another idea - to use a distance from a given cell to closest contact (from global matrix) - so, infinity sign is not zero or one, but continuos!

# round(6600,-3) = 7000; round(6500,-3) = 6000; 
HomologyAndRepeats$FirstWindowWholeKbRes = round(HomologyAndRepeats$FirstWindow,-3)
HomologyAndRepeats$SecondWindowWholeKbRes = round(HomologyAndRepeats$SecondWindow,-3)
table(HomologyAndRepeats$FirstWindowWholeKbRes)
table(HomologyAndRepeats$SecondWindowWholeKbRes)

nrow(HomologyAndRepeats)  # diag 1000: 4005
HomologyAndRepeats = merge(HomologyAndRepeats,GlobalFolding1000, by = c("FirstWindowWholeKbRes","SecondWindowWholeKbRes"))
nrow(HomologyAndRepeats)  # diag 1000: 4005

a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$GlobalFolding1000Score), family = 'binomial')
summary(a)

# diag 500:
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                      -2.13660    0.05004 -42.698  < 2e-16 ***
# scale(HomologyAndRepeats$MicroHomologyScore)      0.29139    0.04774   6.103 1.04e-09 ***
# scale(HomologyAndRepeats$GlobalFolding1000Score)  0.07605    0.04477   1.699   0.0894 .  

# diag 1000:
# Intercept)                                      -2.03721    0.05031 -40.490  < 2e-16 ***
#  scale(HomologyAndRepeats$MicroHomologyScore)      0.29126    0.04806   6.061 1.35e-09 ***
#  scale(HomologyAndRepeats$GlobalFolding1000Score)  0.09128    0.04444   2.054     0.04 *  

##### derive distance to the strongest contact: 6500 vs 14500 (see heatmap: global folding 1 kb resolution). Check coordinates once more!!!
HomologyAndRepeats$DistanceToContact = 0
for (i in 1:nrow(HomologyAndRepeats))
{ # i = 1
  HomologyAndRepeats$DistanceToContact[i] = pointDistance(c(HomologyAndRepeats$FirstWindow[i],HomologyAndRepeats$SecondWindow[i]), c(14550,6550), lonlat = FALSE)  
}
summary(HomologyAndRepeats$DistanceToContact)
summary(HomologyAndRepeats$DistanceToContact) # the closest: -70; the most distant: -8245
a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$DistanceToContact), family = 'binomial')
summary(a)
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                    -2.56654    0.07383 -34.760  < 2e-16 ***
#  scale(HomologyAndRepeats$MicroHomologyScore)  0.43224    0.05269   8.203 2.34e-16 ***
#  scale(HomologyAndRepeats$DistanceToContact)  -1.24881    0.06520 -19.152  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Dispersion parameter for binomial family taken to be 1)
# Null deviance: 2924.7  on 4004  degrees of freedom
# Residual deviance: 2402.3  on 4002  degrees of freedom
# AIC: 2408.3

### derive distance to common repeat:
HomologyAndRepeats$DistanceToContact = 0
for (i in 1:nrow(HomologyAndRepeats))
{ # i = 1
  HomologyAndRepeats$DistanceToContact[i] = pointDistance(c(HomologyAndRepeats$FirstWindow[i],HomologyAndRepeats$SecondWindow[i]), c(13447,8469), lonlat = FALSE)  #  (8469-8482 - 13447-13459)
}
summary(HomologyAndRepeats$DistanceToContact)
summary(HomologyAndRepeats$DistanceToContact) # the closest: -70; the most distant: -8245
a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$DistanceToContact), family = 'binomial')
summary(a)
#                                             Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                                  -2.35793    0.06536 -36.077  < 2e-16 ***
#  scale(HomologyAndRepeats$MicroHomologyScore)  0.27667    0.05038   5.492 3.98e-08 ***
#  scale(HomologyAndRepeats$DistanceToContact)  -1.00602    0.06878 -14.627  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 2924.7  on 4004  degrees of freedom
#Residual deviance: 2609.1  on 4002  degrees of freedom
#AIC: 2615.1

dev.off()  


###### run many log regr in the loop and find the Contact Point with the best AIC 

HomologyAndRepeats$LogRegr.ContactPoint.Coord1 = 0
HomologyAndRepeats$LogRegr.ContactPoint.Coord2 = 0
HomologyAndRepeats$LogRegr.ContactPoint.PiValue = 0
HomologyAndRepeats$LogRegr.ContactPoint.Coeff = 0
HomologyAndRepeats$LogRegr.ContactPoint.AIC = 0
HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance = 0

for (search in 1:nrow(HomologyAndRepeats))
{ # search = 1
Coord1 = HomologyAndRepeats$FirstWindow[search]+50  # big numbers
Coord2 = HomologyAndRepeats$SecondWindow[search]+50 # small numbers
  
HomologyAndRepeats$TempDistanceToContact = 0
for (i in 1:nrow(HomologyAndRepeats))
{ # i = 1
  HomologyAndRepeats$DistanceToContact[i] = pointDistance(c(HomologyAndRepeats$FirstWindow[i],HomologyAndRepeats$SecondWindow[i]), c(Coord1,Coord2), lonlat = FALSE)  
}
summary(HomologyAndRepeats$DistanceToContact) 
a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$DistanceToContact), family = 'binomial')
summary(a)
Res = as.data.frame(summary(a)$coefficients)
HomologyAndRepeats$LogRegr.ContactPoint.PiValue[search] = Res[3,4]
HomologyAndRepeats$LogRegr.ContactPoint.Coeff[search] = Res[3,1]
HomologyAndRepeats$LogRegr.ContactPoint.AIC[search] = a$aic
HomologyAndRepeats$LogRegr.ContactPoint.Coord1[search] = Coord1
HomologyAndRepeats$LogRegr.ContactPoint.Coord2[search] = Coord2
HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance[search] = a$deviance
}

write.table(HomologyAndRepeats,"../../Body/3Results/SlipAndJump.HomologyAndRepeats.txt", sep = '\t')

HomologyAndRepeats = read.table("../../Body/3Results/SlipAndJump.HomologyAndRepeats.txt", sep = '\t')

HomologyAndRepeats = HomologyAndRepeats[order(HomologyAndRepeats$LogRegr.ContactPoint.AIC),]
names(HomologyAndRepeats)
summary(HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance)
summary(HomologyAndRepeats$LogRegr.ContactPoint.AIC)

temp = HomologyAndRepeats[HomologyAndRepeats$LogRegr.ContactPoint.Coord1 == 11950 & HomologyAndRepeats$LogRegr.ContactPoint.Coord2 == 8950,] 
temp


pdf("../../Body/4Figures/SlipAndJump.R.02.pdf")

par(mfrow=c(2,4))

plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord2,HomologyAndRepeats$LogRegr.ContactPoint.AIC, xlab = '5 prime position', ylab = 'AIC'); abline(v = 9000, col = 'red', lwd = 1); abline(v = 6000, col = 'red', lwd = 1) 
plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord2,HomologyAndRepeats$LogRegr.ContactPoint.Coeff, xlab = '5 prime position', ylab = 'Coefficient');  abline(v = 9000, col = 'red', lwd = 1); abline(v = 6000, col = 'red', lwd = 1) 
plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord2,-log10(HomologyAndRepeats$LogRegr.ContactPoint.PiValue),   xlab = '5 prime position', ylab = '-log10(p-value)');  abline(v = 9000, col = 'red', lwd = 1); abline(v = 6000, col = 'red', lwd = 1) 
plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord2,HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance,   xlab = '5 prime position', ylab = 'ResidualDeviance');  abline(v = 9000, col = 'red', lwd = 1); abline(v = 6000, col = 'red', lwd = 1) 

plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord1,HomologyAndRepeats$LogRegr.ContactPoint.AIC, xlab = '3 prime position', ylab = 'AIC');  abline(v = 13000, col = 'red', lwd = 1); abline(v = 16000, col = 'red', lwd = 1); 
plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord1,HomologyAndRepeats$LogRegr.ContactPoint.Coeff, xlab = '3 prime position', ylab = 'Coefficient');  abline(v = 13000, col = 'red', lwd = 1); abline(v = 16000, col = 'red', lwd = 1); 
plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord1,-log10(HomologyAndRepeats$LogRegr.ContactPoint.PiValue),  xlab = '3 prime position', ylab = '-log10(p-value)');  abline(v = 13000, col = 'red', lwd = 1); abline(v = 16000, col = 'red', lwd = 1); 
plot(HomologyAndRepeats$LogRegr.ContactPoint.Coord1,HomologyAndRepeats$LogRegr.ContactPoint.ResidualDeviance,   xlab = '3 prime position', ylab = 'ResidualDeviance');  abline(v = 13000, col = 'red', lwd = 1); abline(v = 16000, col = 'red', lwd = 1) 

dev.off()  

##### Heatmap merged microhomology and AIC scores #### 
if (!require(tidyverse)) install.packages("tidyverse")
if (!require(ggasym)) install.packages("ggasym")

tib <- HomologyAndRepeats %>% 
  select(
    LogRegr.ContactPoint.Coord1,
    LogRegr.ContactPoint.Coord2,
    MicroHomologyScore,
    LogRegr.ContactPoint.AIC
  ) %>% 
  asymmetrise(., 
              LogRegr.ContactPoint.Coord1, 
              LogRegr.ContactPoint.Coord2)

pltHeatmap_mhAIC <- ggplot(tib,
                           aes(x = LogRegr.ContactPoint.Coord1, 
                               y = LogRegr.ContactPoint.Coord2)) +
  geom_asymmat(aes(fill_tl = LogRegr.ContactPoint.AIC,
                   fill_br = MicroHomologyScore)) +
  scale_fill_tl_distiller(
    type = "seq",
    palette = "Spectral",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 1,
      title.position = "top"
    )
  ) +
  scale_fill_br_distiller(
    type = "seq",
    palette = "RdYlGn",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 2,
      title.position = "top"
    )
  ) +
  labs(fill_tl = "top-left Contact AIC",
       fill_br = "bottom-right Microhomology",
       title = "Model of mtDNA contacts") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  )

if (!require(cowplot)) install.packages("cowplot")
cowplot::save_plot(
  plot = pltHeatmap_mhAIC,
  base_height = 8.316,
  base_width = 11.594,
  file = normalizePath(
    file.path(plots_dir, 'heatmap_microhomology_AIC.pdf')
  )
)

##### Again: Heatmap merged microhomology and AIC scores with actual deledions circles ####
tib <- HomologyAndRepeats %>% 
  select(
    LogRegr.ContactPoint.Coord1,
    LogRegr.ContactPoint.Coord2,
    MicroHomologyScore,
    LogRegr.ContactPoint.AIC,
    Deletion
  ) %>% 
  asymmetrise(., 
              LogRegr.ContactPoint.Coord1, 
              LogRegr.ContactPoint.Coord2)

pltHeatmap_mhAIC_wt_deletions <- ggplot(tib,
                                        aes(x = LogRegr.ContactPoint.Coord1, 
                                            y = LogRegr.ContactPoint.Coord2)) +
  geom_asymmat(aes(fill_tl = LogRegr.ContactPoint.AIC,
                   fill_br = MicroHomologyScore)) +
  geom_point(
    data = subset(tib, Deletion > 0),
    aes(alpha = 0.3),
    shape = 1,
    size = 2.3,
    color = "#041c00",
    show.legend = FALSE,
    na.rm = TRUE
  ) +
  scale_fill_tl_distiller(
    type = "seq",
    palette = "Spectral",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 1,
      title.position = "top"
    )
  ) +
  scale_fill_br_distiller(
    type = "seq",
    palette = "RdYlGn",
    direction = 1,
    na.value = "white",
    guide = guide_colourbar(
      direction = "horizontal",
      order = 2,
      title.position = "top"
    )
  ) +
  labs(fill_tl = "top-left Contact AIC",
       fill_br = "bottom-right Microhomology",
       title = "Model of mtDNA contacts and deletions") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  )

cowplot::save_plot(
  plot = pltHeatmap_mhAIC_wt_deletions,
  base_height = 8.316,
  base_width = 11.594,
  file = normalizePath(
    file.path(plots_dir, 'heatmap_microhomology_AIC_wt_deledions.pdf')
  )
)

