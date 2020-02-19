rm(list=ls(all=TRUE))

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
cor.test(DirectRepDensity$Score,MicroHomology$Score, method = 'spearman')
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

a<- glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$DirectRepeatsDensity*HomologyAndRepeats$MicroHomologyScore, family = binomial) 
summary(a)

a<- glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$DirectRepeatsDensity, family = binomial) 
summary(a) # only density of direct repeats doesn't work. 

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

## the matrix is symmetric - I need to keep only one triangle: X>Y (don't need also diagonal, which is made by '500's)
Final$SecondWindow = gsub('X','',Final$SecondWindow)
Final$FirstWindow = as.numeric(as.character(Final$FirstWindow)); Final$SecondWindow = as.numeric(Final$SecondWindow); 
nrow(Final); Final=Final[Final$FirstWindow > Final$SecondWindow,]; nrow(Final)  
GlobalFolding = Final
GlobalFolding = GlobalFolding[order(GlobalFolding$FirstWindow,GlobalFolding$SecondWindow),]

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
GlobalFolding = Final
GlobalFolding = GlobalFolding[order(GlobalFolding$FirstWindow,GlobalFolding$SecondWindow),]
names(GlobalFolding)[3] = c('GlobalFoldingScore');

# GlobalFolding - is the whole genome, not only the major arc!! Keep only major arc in downstream analyses.
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

###### 7: CORRELATE GlobalFolding$Score and InvRepDens$Score - weak positive!
merged = merge(InvRepDens,GlobalFolding, by = c("FirstWindow","SecondWindow"))
summary(merged$FirstWindow)  # 6000 15800
summary(merged$SecondWindow) # 5900 15700
cor.test(merged$InvRepDensScore,merged$GlobalFoldingScore, method = 'spearman') # rho = 0.04716305, p-value = 0.0009027

###### 8: ADD InfinitySign parameter into HomologyAndRepeats dataset:
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
dim(HomologyAndRepeats) # 4950

# is GlobalFoldingScore higher within the cross according to our InfinitySign model? NOT significant, may be try 1000 windows? but should be the same.
wilcox.test(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore,HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore)
t.test(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore,HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore)
summary(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 1,]$GlobalFoldingScore)
summary(HomologyAndRepeats[HomologyAndRepeats$InfinitySign == 0,]$GlobalFoldingScore)

# probably we have to link better global folding and InfinitySign model - till now it was done by eye. Clusterisation? One cluster? 
# dev.off()
###### 9: LOGISTIC REGRESSION: HomologyAndRepeats$Deletion as a function of HomologyAndRepeats$MicroHomologyScore and HomologyAndRepeats$InfinitySign:

a<-glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore + HomologyAndRepeats$InfinitySign, family = 'binomial')
summary(a)

a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$InfinitySign), family = 'binomial')
summary(a)

a<-glm(HomologyAndRepeats$Deletion ~ scale(HomologyAndRepeats$MicroHomologyScore) + scale(HomologyAndRepeats$GlobalFoldingScore), family = 'binomial')
summary(a) # non significant and negative.


# get residuals and correlate them with global matrix

a<-glm(HomologyAndRepeats$Deletion ~ HomologyAndRepeats$MicroHomologyScore, family = 'binomial')
HomologyAndRepeats$Residuals = residuals(a)
summary(HomologyAndRepeats$Residuals)




dev.off()  
