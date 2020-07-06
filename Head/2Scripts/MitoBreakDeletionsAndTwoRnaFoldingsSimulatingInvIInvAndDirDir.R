rm(list=ls(all=TRUE)) # 19 May 2020
library(dplyr)

# CONTENT:
#### 1: READ MITOBREAK AND KEEP ONLY MAJOR ARC DELETIONS:
#### 2: read RnaFolding results, reflecting Gibbs Energy of folding of two dir repeats (100 bp windows) and add contact zone in the file
#### 3. See if DirDirAbsGibbsEnergy explains the presence of deletions in the corresponding window (Logistic or Poisson regression)
#### 4. read RnaFolding results, reflecting Gibbs Energy of folding of two inv repeats (100 bp windows), write down the file
#### 5 . GQUADRUPLEXES CAN WORK AS REPLICATION-PAUSING AGENTS AND/OR HIGH STABILITY DUPLEXES

pdf("../../Body/4Figures/MitoBreakDeletionsAndTwoRnaFoldingsSimulatingInvIInvAndDirDir.R.pdf", width = 20,height = 20)
par(mfrow=c(2,2))  

#### 1: READ MITOBREAK AND KEEP ONLY MAJOR ARC DELETIONS:
breaks = read.table("../../Body/1Raw/MitoBreakDB_12122019.csv", sep = ',', header = TRUE)
breaks$X5..breakpoint = as.numeric(as.character(breaks$X5..breakpoint)); summary(breaks$X5..breakpoint)
breaks$X3..breakpoint = as.numeric(as.character(breaks$X3..breakpoint)); summary(breaks$X3..breakpoint)
breaks = breaks[!is.na(breaks$X3..breakpoint) & !is.na(breaks$X5..breakpoint),]
  
# чтобы не париться можно взять все точки разрыва что больше чем 5781 и меньше чем 16569
nrow(breaks); breaks = breaks[breaks$X5..breakpoint > 5781 & breaks$X3..breakpoint > 5781,]; nrow(breaks)
summary(breaks$X5..breakpoint) # start
summary(breaks$X3..breakpoint) # end

#### 2: read RnaFolding results, reflecting Gibbs Energy of folding of two dir repeats (100 bp windows) and add contact zone in the file
InputDir = "../../Body/1Raw/VictorFoldings/Reverse1_100bp.x10.Complement2_100bp/result/"
VecOfFiles = list.files(path = InputDir); length(VecOfFiles) 
DirDirFoldings = data.frame()
for (i in 1:length(VecOfFiles))
{ # i = 1
  if (file.info(paste(InputDir,VecOfFiles[i],sep = ''))$size > 0)
  {
  Gibbs = read.table(paste(InputDir,VecOfFiles[i],sep = ''),sep = '\t')
  GibbsEnergy = gsub(".*\\(",'',Gibbs$V1[3]); GibbsEnergy = gsub("\\)",'',GibbsEnergy); 
  Win1Win2 = gsub(".csv",'',VecOfFiles[i])
  Win1 = as.numeric(gsub("_.*",'',Win1Win2))
  Win2 = as.numeric(gsub(".*_",'',Win1Win2))
  OneLine = data.frame(Win1,Win2,as.numeric(GibbsEnergy))
  DirDirFoldings = rbind(DirDirFoldings,OneLine)
  }
}
names(DirDirFoldings)=c('Win1','Win2','DirDirGibbsEnergy')
summary(DirDirFoldings$DirDirGibbsEnergy)
DirDirFoldings$DirDirAbsGibbsEnergy = -1*DirDirFoldings$DirDirGibbsEnergy
summary(DirDirFoldings$DirDirAbsGibbsEnergy)

## add contact zone info
DirDirFoldings$ContactZone = 0
for (i in 1:nrow(DirDirFoldings))
{ # i=1
  if (DirDirFoldings$Win1[i] >= 6000 & DirDirFoldings$Win1[i] <= 9000 & DirDirFoldings$Win2[i] >= 13000 &  DirDirFoldings$Win2[i] <= 16000)
  {DirDirFoldings$ContactZone[i] = 1}
}
table(DirDirFoldings$ContactZone)

# within the contact zone DirDirAbsGibbsEnergy is a bit higher:
summary(DirDirFoldings[DirDirFoldings$ContactZone == 0,]$DirDirAbsGibbsEnergy)
summary(DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DirDirAbsGibbsEnergy)
wilcox.test(DirDirFoldings[DirDirFoldings$ContactZone == 0,]$DirDirAbsGibbsEnergy,DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DirDirAbsGibbsEnergy)
t.test(DirDirFoldings[DirDirFoldings$ContactZone == 0,]$DirDirAbsGibbsEnergy,DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DirDirAbsGibbsEnergy)

##### 3. See if DirDirAbsGibbsEnergy explains the presence of deletions in the corresponding window (Logistic or Poisson regression)
# walk along each line of the DirDirFoldings and check - how many deletions in the corresponding window
DirDirFoldings$NumbOfDeletions = 0
DirDirFoldings$DummyDeletions = 0
for (i in 1:nrow(DirDirFoldings))
{ # i = 1
  Win1 = DirDirFoldings$Win1[i];
  Win2 = DirDirFoldings$Win2[i]; 
  temp = breaks[breaks$X5..breakpoint >= Win1 & breaks$X5..breakpoint < (Win1+100) & breaks$X3..breakpoint >= Win2 & breaks$X3..breakpoint < (Win2+100),]
  DirDirFoldings$NumbOfDeletions[i] = nrow(temp)
  if (nrow(temp)>0) {DirDirFoldings$DummyDeletions[i] = 1}
}
summary(DirDirFoldings$NumbOfDeletions)
summary(DirDirFoldings[DirDirFoldings$NumbOfDeletions > 0,]$NumbOfDeletions)
table(DirDirFoldings$DummyDeletions)

cor.test(DirDirFoldings$DirDirAbsGibbsEnergy,DirDirFoldings$NumbOfDeletions, method = 'spearman') # positive
cor.test(DirDirFoldings$DirDirAbsGibbsEnergy,DirDirFoldings$DummyDeletions, method = 'spearman') # positive

summary(glm(DirDirFoldings$NumbOfDeletions ~ DirDirFoldings$DirDirAbsGibbsEnergy, family = poisson())) # weak
summary(glm(DirDirFoldings$DummyDeletions ~ DirDirFoldings$DirDirAbsGibbsEnergy, family = binomial())) # good

## the final correlation to digest:
summary(glm(DirDirFoldings$DummyDeletions ~ scale(DirDirFoldings$DirDirAbsGibbsEnergy) + scale(DirDirFoldings$ContactZone), family = binomial())) # good

# other correlations:
summary(glm(DirDirFoldings$NumbOfDeletions ~ scale(DirDirFoldings$DirDirAbsGibbsEnergy) + scale(DirDirFoldings$ContactZone), family = poisson())) # weak
summary(glm(DirDirFoldings$DummyDeletions ~ scale(DirDirFoldings$DirDirAbsGibbsEnergy)*scale(DirDirFoldings$ContactZone), family = binomial())) # good
summary(glm(DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DummyDeletions ~ DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DirDirAbsGibbsEnergy, family = binomial())) # with Khrapko we used similar (more broad still...) dataset and results were significant  - do I need to play with different subsets of deletions - Falkenberg's data????

### why poisson regression describes data poorer than binomial?????!!!!!!!
boxplot(DirDirFoldings[DirDirFoldings$DummyDeletions == 0,]$DirDirAbsGibbsEnergy,DirDirFoldings[DirDirFoldings$DummyDeletions == 1,]$DirDirAbsGibbsEnergy, names=c('0','>0'), notch = TRUE)
boxplot(DirDirFoldings[DirDirFoldings$NumbOfDeletions == 0,]$DirDirAbsGibbsEnergy,DirDirFoldings[DirDirFoldings$NumbOfDeletions == 1,]$DirDirAbsGibbsEnergy, DirDirFoldings[DirDirFoldings$NumbOfDeletions == 2,]$DirDirAbsGibbsEnergy, DirDirFoldings[DirDirFoldings$NumbOfDeletions > 2,]$DirDirAbsGibbsEnergy, names=c('0','1','2','>2'), notch = TRUE, outline = FALSE, varwidth =  TRUE)

# probably we need to delete close windows, because there are NO short deletions
# making minimal distance from 1000 to 2000 all results are improving, but from other side,
# increasing distance between windows we artificially enrich our dataset to long deletions (and deletions in contact zone)
summary(DirDirFoldings$Win2 - DirDirFoldings$Win1) # min = 1000. make it 2000 - everything is a bit better
DirDirFoldingsLong = DirDirFoldings[DirDirFoldings$Win2 - DirDirFoldings$Win1 > 2000,]
cor.test(DirDirFoldingsLong$DirDirAbsGibbsEnergy,DirDirFoldingsLong$NumbOfDeletions, method = 'spearman') # positive
cor.test(DirDirFoldingsLong$DirDirAbsGibbsEnergy,DirDirFoldingsLong$DummyDeletions, method = 'spearman') # positive
summary(glm(DirDirFoldingsLong$NumbOfDeletions ~ DirDirFoldingsLong$DirDirAbsGibbsEnergy, family = poisson())) # weak
summary(glm(DirDirFoldingsLong$DummyDeletions ~ DirDirFoldingsLong$DirDirAbsGibbsEnergy, family = binomial())) # good 

###### 4 read RnaFolding results, reflecting Gibbs Energy of folding of two inv repeats (100 bp windows)
InputDir = "../../Body/1Raw/VictorFoldings/Complement1_100bp.x10.Complement2_100bp/result/"
VecOfFiles = list.files(path = InputDir); length(VecOfFiles) 
InvInvFoldings = data.frame()
for (i in 1:length(VecOfFiles))
{ # i = 1
  if (file.info(paste(InputDir,VecOfFiles[i],sep = ''))$size > 0)
  {
    Gibbs = read.table(paste(InputDir,VecOfFiles[i],sep = ''),sep = '\t')
    GibbsEnergy = gsub(".*\\(",'',Gibbs$V1[3]); GibbsEnergy = gsub("\\)",'',GibbsEnergy); 
    Win1Win2 = gsub(".csv",'',VecOfFiles[i])
    Win1 = as.numeric(gsub("_.*",'',Win1Win2))
    Win2 = as.numeric(gsub(".*_",'',Win1Win2))
    OneLine = data.frame(Win1,Win2,as.numeric(GibbsEnergy))
    InvInvFoldings = rbind(InvInvFoldings,OneLine)
  }
}
names(InvInvFoldings)=c('Win1','Win2','InvInvGibbsEnergy')
summary(InvInvFoldings$InvInvGibbsEnergy)
InvInvFoldings$InvInvAbsGibbsEnergy = -1*InvInvFoldings$InvInvGibbsEnergy
summary(InvInvFoldings$InvInvAbsGibbsEnergy)

write.table(InvInvFoldings, "../../Body/2Derived/InvInvFoldings.Complement1_100bp.x10.Complement2_100bp.txt")

## add contact zone info
InvInvFoldings$ContactZone = 0
for (i in 1:nrow(InvInvFoldings))
{ # i=1
  if (InvInvFoldings$Win1[i] >= 6000 & InvInvFoldings$Win1[i] <= 9000 & InvInvFoldings$Win2[i] >= 13000 &  InvInvFoldings$Win2[i] <= 16000)
  {InvInvFoldings$ContactZone[i] = 1}
}
table(InvInvFoldings$ContactZone)

# within the contact zone InvInvAbsGibbsEnergy is a bit higher (not significant):
summary(InvInvFoldings[DirDirFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy)
summary(InvInvFoldings[DirDirFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy)
wilcox.test(InvInvFoldings[DirDirFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy,InvInvFoldings[InvInvFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy)
t.test(InvInvFoldings[InvInvFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy,InvInvFoldings[InvInvFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy)

# test if InvInv directly associated with deletions:
InvInvFoldings$NumbOfDeletions = 0
InvInvFoldings$DummyDeletions = 0
for (i in 1:nrow(DirDirFoldings))
{ # i = 1
  Win1 = InvInvFoldings$Win1[i];
  Win2 = InvInvFoldings$Win2[i]; 
  #temp = breaks[breaks$X5..breakpoint >= Win1 & breaks$X5..breakpoint < (Win1+100) & breaks$X3..breakpoint >= Win2 & breaks$X3..breakpoint < (Win2+100),]
  temp = breaks[breaks$X5..breakpoint >= (Win1-100) & breaks$X5..breakpoint < (Win1) & breaks$X3..breakpoint >= (Win2+100) & breaks$X3..breakpoint < (Win2+200),] # deletions should be before the first and after the second window
  InvInvFoldings$NumbOfDeletions[i] = nrow(temp)
  if (nrow(temp)>0) {InvInvFoldings$DummyDeletions[i] = 1}
}
summary(InvInvFoldings$NumbOfDeletions)
summary(InvInvFoldings[InvInvFoldings$NumbOfDeletions > 0,]$NumbOfDeletions)
table(InvInvFoldings$DummyDeletions)

cor.test(InvInvFoldings$InvInvAbsGibbsEnergy,InvInvFoldings$NumbOfDeletions, method = 'spearman') # NOTHING
cor.test(InvInvFoldings$InvInvAbsGibbsEnergy,InvInvFoldings$DummyDeletions, method = 'spearman') # NOTHING 
summary(glm(InvInvFoldings$NumbOfDeletions ~ InvInvFoldings$InvInvAbsGibbsEnergy, family = poisson())) # WEAK POSITIVE
summary(glm(InvInvFoldings$DummyDeletions ~ InvInvFoldings$InvInvAbsGibbsEnergy, family = binomial())) # NOTHING

# THE SAME WITHIN THE CONTACT ZONE:
cor.test(InvInvFoldings[InvInvFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy,InvInvFoldings[InvInvFoldings$ContactZone == 1,]$NumbOfDeletions, method = 'spearman') # NOTHING
cor.test(InvInvFoldings[InvInvFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy,InvInvFoldings[InvInvFoldings$ContactZone == 1,]$DummyDeletions, method = 'spearman') # NOTHING
summary(glm(InvInvFoldings[InvInvFoldings$ContactZone == 1,]$NumbOfDeletions ~ InvInvFoldings[InvInvFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy, family = poisson())) # POSITIVE
summary(glm(InvInvFoldings[InvInvFoldings$ContactZone == 1,]$DummyDeletions ~ InvInvFoldings[InvInvFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy, family = binomial())) # NOTHING
# THE SAME OUTSIDE THE CONTACT ZONE:
cor.test(InvInvFoldings[InvInvFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy,InvInvFoldings[InvInvFoldings$ContactZone == 0,]$NumbOfDeletions, method = 'spearman') # NOTHING
cor.test(InvInvFoldings[InvInvFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy,InvInvFoldings[InvInvFoldings$ContactZone == 0,]$DummyDeletions, method = 'spearman') # NOTHING
summary(glm(InvInvFoldings[InvInvFoldings$ContactZone == 0,]$NumbOfDeletions ~ InvInvFoldings[InvInvFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy, family = poisson())) # WEAK NEGATIVE
summary(glm(InvInvFoldings[InvInvFoldings$ContactZone == 0,]$DummyDeletions ~ InvInvFoldings[InvInvFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy, family = binomial())) # NOTHING

###### MERGE
## to merge InvInv with DirDir I need to make a shift - to assure that InvInv will be nested within DirDir
InvInvFoldings = select(InvInvFoldings,Win1,Win2,InvInvAbsGibbsEnergy)
InvInvFoldings$Win1 = InvInvFoldings$Win1-100 # allign with DirDir windows:
InvInvFoldings$Win2 = InvInvFoldings$Win2+100

DirDirFoldings = merge(DirDirFoldings,InvInvFoldings, by = c('Win1','Win2'))

##### ANALYSES:

summary(glm(DirDirFoldings$DummyDeletions ~ scale(DirDirFoldings$DirDirAbsGibbsEnergy)+ scale(DirDirFoldings$InvInvAbsGibbsEnergy), family = binomial()))
summary(glm(DirDirFoldings$NumbOfDeletions ~ scale(DirDirFoldings$DirDirAbsGibbsEnergy)+ scale(DirDirFoldings$InvInvAbsGibbsEnergy), family = poisson())) # two are significant!!
# the best result:
summary(glm(DirDirFoldings$NumbOfDeletions ~ scale(DirDirFoldings$DirDirAbsGibbsEnergy)*scale(DirDirFoldings$InvInvAbsGibbsEnergy), family = poisson())) # three are significant - interaction is negative and very significant!!!

### DirDir and InvInv stabilities negatively correlate with each other (result of negative selection???) in the whole major arc, in contact zone and outside the contact zone:
cor.test(DirDirFoldings$DirDirAbsGibbsEnergy,DirDirFoldings$InvInvAbsGibbsEnergy, method = 'spearman') # negative strong!!
plot(log2(DirDirFoldings$DirDirAbsGibbsEnergy+1),log2(DirDirFoldings$InvInvAbsGibbsEnergy+1))
cor.test(DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DirDirAbsGibbsEnergy,DirDirFoldings[DirDirFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy, method = 'spearman') # negative strong!!
plot(log2(DirDirFoldings[DirDirFoldings$ContactZone == 1,]$DirDirAbsGibbsEnergy+1),log2(DirDirFoldings[DirDirFoldings$ContactZone == 1,]$InvInvAbsGibbsEnergy+1))
cor.test(DirDirFoldings[DirDirFoldings$ContactZone == 0,]$DirDirAbsGibbsEnergy,DirDirFoldings[DirDirFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy, method = 'spearman') # negative strong!!
plot(log2(DirDirFoldings[DirDirFoldings$ContactZone == 0,]$DirDirAbsGibbsEnergy+1),log2(DirDirFoldings[DirDirFoldings$ContactZone == 0,]$InvInvAbsGibbsEnergy+1))

## effect of more long deletion and/or contact zone:
summary(glm(NumbOfDeletions ~ scale(DirDirAbsGibbsEnergy)*scale(InvInvAbsGibbsEnergy), family = poisson(), data = DirDirFoldings))

summary(glm(DummyDeletions ~ scale(DirDirAbsGibbsEnergy)*scale(InvInvAbsGibbsEnergy), family = binomial(), data = DirDirFoldings))
summary(glm(DummyDeletions ~ scale(DirDirAbsGibbsEnergy)+scale(InvInvAbsGibbsEnergy), family = binomial(), data = DirDirFoldings))

#### 5 . GQUADRUPLEXES CAN WORK AS REPLICATION-PAUSING AGENTS AND/OR HIGH STABILITY DUPLEXES
# G-quadruplex secondary structures (G4) are formed in nucleic acids by sequences that are rich in guanine.
# They are helical in shape and contain guanine tetrads that can form from one (it should be case in mtDNA), two or four strands.
# Pqsfinder detects DNA and RNA sequence patterns that are likely to fold into an intramolecular G-quadruplex (G4).
# Unlike many other approaches, pqsfinder is able to detect G4s folded from imperfect G-runs containing bulges or mismatches or G4s having long loops. 
# Pqsfinder also assigns an integer score to each hit that was fitted on G4 sequencing data and corresponds to expected stability of the folded G4
# here on pages 6 and 7 there is explanation of outputs

GQ = read.table("../../Body/1Raw/QuadruplexFormingSequences/Homo_sapiens.genome.cut.gff", sep = '\t', header = FALSE)
GQ=GQ[,-c(1,2)]
names(GQ)=c('start','end','score','strand','type')
GQ = GQ[order(GQ$score),]
GQ = GQ[GQ$strand == '-',] # 99% of all quadruplexes are on "-" heavy chain.
GQ$length = GQ$end-GQ$start
summary(GQ$length) # max = 49
plot(GQ$start,GQ$score) # vertical lines mean many GQ with similar start and different scores 
plot(GQ$start,GQ$score,xlim=c(6000,16560)) # vertical lines mean many GQ with similar start and different scores 

# count overlap of each window with Gquadruplex = how many nucleotides are covered by Qquadruplex, what is the average and maximal score in each window?
DirDirFoldings$Win1GQOverlap = 0
DirDirFoldings$Win1GQAverageScore = 0
DirDirFoldings$Win1GQMaxScore = 0

DirDirFoldings$Win2GQOverlap = 0
DirDirFoldings$Win2GQAverageScore = 0
DirDirFoldings$Win2GQMaxScore = 0

for (i in 1:nrow(DirDirFoldings))
{ # i = 1
  # either start or end (max length = 49) should be inside 100 bp window:
  temp = GQ[(GQ$start >= DirDirFoldings$Win1[i] & GQ$start < (DirDirFoldings$Win1[i] + 100)) | (GQ$end >= DirDirFoldings$Win1[i] & GQ$end < (DirDirFoldings$Win1[i] + 100)),]
  if (nrow(temp)>0)
  {
  DirDirFoldings$Win1GQAverageScore[i] = mean(temp$score)
  DirDirFoldings$Win1GQMaxScore[i] = max(temp$score) 
  MaxEnd = max(temp$end); MinStart = min(temp$start)
  window=seq(DirDirFoldings$Win1[i],DirDirFoldings$Win1[i]+100)
  gquadruplex = seq(MinStart,MaxEnd)
  DirDirFoldings$Win1GQOverlap[i] = length(intersect(window,gquadruplex))
  }
  temp = GQ[(GQ$start >= DirDirFoldings$Win2[i] & GQ$start < (DirDirFoldings$Win2[i] + 100)) | (GQ$end >= DirDirFoldings$Win2[i] & GQ$end < (DirDirFoldings$Win2[i] + 100)),]
  if (nrow(temp)>0)
  {
    DirDirFoldings$Win2GQAverageScore[i] = mean(temp$score)
    DirDirFoldings$Win2GQMaxScore[i] = max(temp$score) 
    MaxEnd = max(temp$end); MinStart = min(temp$start)
    window=seq(DirDirFoldings$Win2[i],DirDirFoldings$Win2[i]+100)
    gquadruplex = seq(MinStart,MaxEnd)
    DirDirFoldings$Win2GQOverlap[i] = length(intersect(window,gquadruplex))
  }
}

### GQ may increase in frequency/strength with the mtDNA location because of A>G mutation shift: NOT REALLY - SUPER WEAK NON-SIGNIFICANT TREND
summary(DirDirFoldings$Win1)
summary(DirDirFoldings$Win2) # the second win is shifted a bit on the right.

GqVsLocationWin1 = select(DirDirFoldings,Win1,Win1GQOverlap,Win1GQAverageScore,Win1GQMaxScore)
names(GqVsLocationWin1)=c('Win','WinGQOverlap','WinGQAverageScore','WinGQMaxScore')
GqVsLocationWin1 = unique(GqVsLocationWin1); GqVsLocationWin1 = GqVsLocationWin1[order(GqVsLocationWin1$Win),]

GqVsLocationWin2 = select(DirDirFoldings,Win2,Win2GQOverlap,Win2GQAverageScore,Win2GQMaxScore)
names(GqVsLocationWin2)=c('Win','WinGQOverlap','WinGQAverageScore','WinGQMaxScore')
GqVsLocationWin2 = unique(GqVsLocationWin2); GqVsLocationWin2 = GqVsLocationWin2[order(GqVsLocationWin2$Win),]

GqVsLocation = rbind(GqVsLocationWin1,GqVsLocationWin2)
GqVsLocation = unique(GqVsLocation)

cor.test(GqVsLocation$Win,GqVsLocation$WinGQOverlap, method = 'spearman')
cor.test(GqVsLocation$Win,GqVsLocation$WinGQAverageScore, method = 'spearman')
cor.test(GqVsLocation$Win,GqVsLocation$WinGQMaxScore, method = 'spearman')

cor.test(GqVsLocation[GqVsLocation$WinGQAverageScore > 0,]$Win,GqVsLocation[GqVsLocation$WinGQAverageScore > 0,]$WinGQAverageScore, method = 'spearman')
cor.test(GqVsLocation[GqVsLocation$WinGQAverageScore > 0,]$Win,GqVsLocation[GqVsLocation$WinGQAverageScore > 0,]$WinGQMaxScore, method = 'spearman')
plot(GqVsLocation$Win,GqVsLocation$WinGQMaxScore)
plot(GqVsLocation$Win,GqVsLocation$WinGQOverlap) # dev.off()
plot(GqVsLocation$Win,GqVsLocation$WinGQAverageScore) # dev.off()
plot(GqVsLocation$WinGQOverlap,GqVsLocation$WinGQAverageScore) # dev.off()

### CONTACT ZONE MIGHT BE ENRICHED IN MORE NUMERIC/MORE STRONG GQ: NO
wilcox.test(GqVsLocation[(GqVsLocation$Win >= 6000 & GqVsLocation$Win <= 9000) | (GqVsLocation$Win >= 13000 & GqVsLocation$Win <= 16000),]$WinGQOverlap,GqVsLocation[GqVsLocation$Win >= 9000 & GqVsLocation$Win < 13000,]$WinGQOverlap)
wilcox.test(GqVsLocation[(GqVsLocation$Win >= 6000 & GqVsLocation$Win <= 9000) | (GqVsLocation$Win >= 13000 & GqVsLocation$Win <= 16000),]$WinGQAverageScore,GqVsLocation[GqVsLocation$Win >= 9000 & GqVsLocation$Win < 13000,]$WinGQAverageScore)
wilcox.test(GqVsLocation[(GqVsLocation$Win >= 6000 & GqVsLocation$Win <= 9000) | (GqVsLocation$Win >= 13000 & GqVsLocation$Win <= 16000),]$WinGQMaxScore,GqVsLocation[GqVsLocation$Win >= 9000 & GqVsLocation$Win < 13000,]$WinGQMaxScore)

### IF GQ PER SE ARE ASSOCIATED WITH DELETIONS: ONLY THE START OF DELETIONS (WIN1) STRONGLY DEPENDS ON GQ: ASYMMETRY -> FALKENBERG IS RIGHT
summary(glm(NumbOfDeletions ~ scale(Win1GQOverlap) + scale(Win2GQOverlap), family = poisson(), data = DirDirFoldings))
summary(glm(NumbOfDeletions ~ scale(Win1GQOverlap), family = poisson(), data = DirDirFoldings)) # strong positive !!!! AS IF GQ PAUSES REPLICATION
summary(glm(NumbOfDeletions ~ scale(Win2GQOverlap), family = poisson(), data = DirDirFoldings)) # nothing

plot(DirDirFoldings$Win1GQOverlap,DirDirFoldings$NumbOfDeletions)

ResVec=data.frame()
for (i in 0:99) 
{
  MeanWin1 = mean(DirDirFoldings[DirDirFoldings$Win1GQOverlap > i,]$NumbOfDeletions)
  MeanWin2 = mean(DirDirFoldings[DirDirFoldings$Win2GQOverlap > i,]$NumbOfDeletions)
  MeanWin1Dummy = mean(DirDirFoldings[DirDirFoldings$Win1GQOverlap > i,]$DummyDeletions)
  MeanWin2Dummy = mean(DirDirFoldings[DirDirFoldings$Win2GQOverlap > i,]$DummyDeletions)
  ResVec=rbind(ResVec,data.frame(i,MeanWin1,MeanWin2,MeanWin1Dummy,MeanWin2Dummy)) 
}
names(ResVec)=c('ThresholdOfGqCoverage','MeanNumberOfDeletionsWin1','MeanNumberOfDeletionsWin2','MeanNumberOfDeletionsWin1Dummy','MeanNumberOfDeletionsWin2Dummy')
str(ResVec)
par(mfrow=c(2,2))
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin1)
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin2)
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin1Dummy)
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin2Dummy)

#par(mfrow=c(1,1))
YMin = min(min(ResVec$MeanNumberOfDeletionsWin1),min(ResVec$MeanNumberOfDeletionsWin2))
YMax = max(max(ResVec$MeanNumberOfDeletionsWin1),max(ResVec$MeanNumberOfDeletionsWin2))
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin1, xlim = c(0,100), ylim = c(YMin,YMax), col = 'red', xlab = '', ylab = ''); par(new = TRUE)
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin2, xlim = c(0,100), ylim = c(YMin,YMax), col = 'green', xlab = 'minimal overlap with G-Quadruplexes', ylab = 'average number of deletions')
legend(0,0.40,legend=c('start','end'), col = c('red','green'), pch = 16)

YMin = min(min(ResVec$MeanNumberOfDeletionsWin1Dummy),min(ResVec$MeanNumberOfDeletionsWin2Dummy))
YMax = max(max(ResVec$MeanNumberOfDeletionsWin1Dummy),max(ResVec$MeanNumberOfDeletionsWin2Dummy))
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin1Dummy, xlim = c(0,100), ylim = c(YMin,YMax), col = 'red', xlab = '', ylab = ''); par(new = TRUE)
plot(ResVec$ThresholdOfGqCoverage,ResVec$MeanNumberOfDeletionsWin2Dummy, xlim = c(0,100), ylim = c(YMin,YMax), col = 'green', xlab = 'minimal overlap with G-Quadruplexes', ylab = 'fraction of windows with at least one deletion')
legend(0,0.20,legend=c('start','end'), col = c('red','green'), pch = 16)

### IF FALKENBERG IS RIGHT WE EXPECT GQ BEING LOCATED AFTER THE START OF THE DELETION: SO-SO - think more!!! start of the deletion is more often within the GQ!!!
# find for each deletion (start and end) the closest GQ (start and end)
breaks$TheClosestGQToDelStart = 0
breaks$TheClosestGQToDelEnd = 0
for (i in 1:nrow(breaks))
{ # i = 125
  GQ$DelStart = breaks$X5..breakpoint[i]; GQ$DelEnd = breaks$X3..breakpoint[i];
  # DelStart
  GQ$MinStart1 = abs(GQ$start - GQ$DelStart) 
  GQ$MinStart2 = abs(GQ$end - GQ$DelStart)
  GQ$MinStart = apply(data.frame(GQ$MinStart1,GQ$MinStart2),1,FUN = min)
  # GQ$MinStart = min(GQ$MinStart1,GQ$MinStart2)
  GQ = GQ[order(GQ$MinStart),]
  temp = GQ[1,]
  if (temp$DelStart > temp$start & temp$DelStart > temp$end) {breaks$TheClosestGQToDelStart[i] = temp$end - temp$DelStart}
  if (temp$DelStart > temp$start & temp$DelStart < temp$end) {breaks$TheClosestGQToDelStart[i] = 0}
  if (temp$DelStart < temp$start & temp$DelStart < temp$end) {breaks$TheClosestGQToDelStart[i] = temp$start - temp$DelStart}
  
  # DelEnd: 
  GQ$MinStart1 = abs(GQ$start - GQ$DelEnd) 
  GQ$MinStart2 = abs(GQ$end - GQ$DelEnd)
  GQ$MinStart = apply(data.frame(GQ$MinStart1,GQ$MinStart2),1,FUN = min)
  # GQ$MinStart = min(GQ$MinStart1,GQ$MinStart2)
  GQ = GQ[order(GQ$MinStart),]
  temp = GQ[1,]
  if (temp$DelEnd > temp$start & temp$DelEnd > temp$end) {breaks$TheClosestGQToDelEnd[i] = temp$end - temp$DelEnd}
  if (temp$DelEnd > temp$start & temp$DelEnd < temp$end) {breaks$TheClosestGQToDelEnd[i] = 0}
  if (temp$DelEnd < temp$start & temp$DelEnd < temp$end) {breaks$TheClosestGQToDelEnd[i] = temp$start - temp$DelEnd}
}

hist(breaks$TheClosestGQToDelStart, breaks = 100)
hist(breaks$TheClosestGQToDelEnd, breaks = 100)
wilcox.test(breaks$TheClosestGQToDelStart,breaks$TheClosestGQToDelEnd, paired = TRUE)
t.test(breaks$TheClosestGQToDelStart,breaks$TheClosestGQToDelEnd, paired = TRUE)
mean(breaks$TheClosestGQToDelStart)
mean(breaks$TheClosestGQToDelEnd)

### GQuadruplexes and Foldings (very positively with DirDirAbsGibbsEnergy and very negatively with InvInvAbsGibbsEnergy):
# because many GGGG on one side and reversed many CCCC - strong contact of course 
cor.test(DirDirFoldings$Win1GQOverlap,DirDirFoldings$DirDirAbsGibbsEnergy, method = 'spearman') # very positive
cor.test(DirDirFoldings$Win2GQOverlap,DirDirFoldings$DirDirAbsGibbsEnergy, method = 'spearman') # super positive
cor.test((DirDirFoldings$Win1GQOverlap + DirDirFoldings$Win2GQOverlap)/2,DirDirFoldings$DirDirAbsGibbsEnergy, method = 'spearman') # super positive
cor.test((DirDirFoldings$Win1GQAverageScore + DirDirFoldings$Win1GQAverageScore)/2,DirDirFoldings$DirDirAbsGibbsEnergy, method = 'spearman') # super positive
cor.test((DirDirFoldings$Win1GQMaxScore + DirDirFoldings$Win1GQMaxScore)/2,DirDirFoldings$DirDirAbsGibbsEnergy, method = 'spearman') # super positive

## DirDirFoldings$InvInvAbsGibbsEnergy negatively correlates to GQuadruplexes:
# because GQ are not inverted repeats
cor.test(DirDirFoldings$Win1GQOverlap,DirDirFoldings$InvInvAbsGibbsEnergy, method = 'spearman') # very negative
cor.test(DirDirFoldings$Win2GQOverlap,DirDirFoldings$InvInvAbsGibbsEnergy, method = 'spearman') # very negative

# distribution of GQuadruplexes explains negative correlation between DirDir and InvInv? YES, IT IS!!!!!!!
summary(glm(Win1GQOverlap ~ scale(DirDirAbsGibbsEnergy) + scale(InvInvAbsGibbsEnergy), family = poisson(), data = DirDirFoldings))
summary(glm(Win2GQOverlap ~ scale(DirDirAbsGibbsEnergy) + scale(InvInvAbsGibbsEnergy), family = poisson(), data = DirDirFoldings))

# very high AIC => plot it somehow
DirDirFoldings$Win12GQOverlap = (DirDirFoldings$Win1GQOverlap + DirDirFoldings$Win2GQOverlap)/2
ev3 = DirDirFoldings$Win12GQOverlap; summary(ev3); ev3 = ev3/200;  
# plot(DirDirFoldings$DirDirAbsGibbsEnergy,DirDirFoldings$InvInvAbsGibbsEnergy,pch.cex = DirDirFoldings$Win12GQOverlap)
# with(DirDirFoldings, symbols(x=log2(DirDirFoldings$DirDirAbsGibbsEnergy+1), log2(DirDirFoldings$InvInvAbsGibbsEnergy+1), circles=ev3, inches=1/4,ann=F, fg=rgb(0,0,1,0.1), bg=NULL, xlab = 'DirDir', ylab = 'InvInv'))
par(mfrow=c(2,2))
plot(DirDirFoldings$Win12GQOverlap,DirDirFoldings$DirDirAbsGibbsEnergy)
plot(DirDirFoldings$Win12GQOverlap,DirDirFoldings$InvInvAbsGibbsEnergy)

dev.off()



# put altogether:
# the best result: Win1GQOverlap correlates with DirDirAbsGibbsEnergy 
# => should we interpret GQ as high GibbsEnergy of direct repeats or as pausing of replication only at first window!!!!!?????
summary(glm(log2(NumbOfDeletions+1) ~ log2(Win1GQOverlap+1) + log2(Win2GQOverlap+1), family = poisson(), data = DirDirFoldings))
summary(glm(log2(NumbOfDeletions+1) ~ log2(Win1GQOverlap+1) + log2(InvInvAbsGibbsEnergy+1) + log2(DirDirAbsGibbsEnergy+1), family = poisson(), data = DirDirFoldings))
summary(glm(log2(NumbOfDeletions+1) ~ log2(Win1GQOverlap+1) + log2(InvInvAbsGibbsEnergy+1)*log2(DirDirAbsGibbsEnergy+1), family = poisson(), data = DirDirFoldings))

summary(glm(NumbOfDeletions ~ log2(Win1GQOverlap+1) + log2(InvInvAbsGibbsEnergy+1) + log2(DirDirAbsGibbsEnergy+1), family = poisson(), data = DirDirFoldings))
summary(glm(NumbOfDeletions ~ scale(log2(Win1GQOverlap+1)) + scale(log2(InvInvAbsGibbsEnergy+1))*scale(log2(DirDirAbsGibbsEnergy+1)), family = poisson(), data = DirDirFoldings))

# WORKING HYPOTHESIS
# DirDirAbs ~ +GQ => it makes stability higher!
# InvInvAbs ~ -GQ => because GQ doesn't interacts with another GQ
## Dir-GQ-InvInv-Dir ### THE BEST COMBINATION!!!!!


###### TO DO:
### olivier presentation = add many components just to explain more and more....

dev.off()

