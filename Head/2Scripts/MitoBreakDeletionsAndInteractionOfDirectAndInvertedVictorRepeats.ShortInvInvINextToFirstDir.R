rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.ShortInvInvNextToFirstDir.R01.pdf")
  
##### 1: READ MITOBREAK AND KEEP ONLY MAJOR ARC DELETIONS:
  breaks = read.table("../../Body/1Raw/MitoBreakDB_12122019.csv", sep = ',', header = TRUE)
  breaks$X5..breakpoint = as.numeric(as.character(breaks$X5..breakpoint)); summary(breaks$X5..breakpoint)
  breaks$X3..breakpoint = as.numeric(as.character(breaks$X3..breakpoint)); summary(breaks$X3..breakpoint)
  breaks = breaks[!is.na(breaks$X3..breakpoint) & !is.na(breaks$X5..breakpoint),]
  
  # поскольку координаты не такие простые (см ниже) - чтобы не париться можно взять все точки разрыва что больше чем 5781 и меньше чем 16569
  # OH: 110-441
  # OL: 5721-5781
  for (i in 1:nrow(breaks))
  {  
    if (breaks$X5..breakpoint[i] < 110) {breaks$X5..breakpoint[i] = breaks$X5..breakpoint[i] + 16569}
    if (breaks$X3..breakpoint[i] < 110) {breaks$X3..breakpoint[i] = breaks$X3..breakpoint[i] + 16569}
  }
  summary(breaks$X5..breakpoint)
  summary(breaks$X3..breakpoint)
  
  nrow(breaks); breaks = breaks[breaks$X5..breakpoint > 5781 & breaks$X3..breakpoint > 5781,]; nrow(breaks)
  summary(breaks$X5..breakpoint)
  summary(breaks$X3..breakpoint)

##### 2: read all degraded repeats of Victor (direct and inverted):
  Rep = read.table("../../Body/1Raw/DegRepsHomo_sapiens 2018_07_02.csv", header = TRUE, sep = ';') 
  #Rep = Rep[Rep$id_type == 1 | Rep$id_type == 4,] # 1 = direct, 4 = inverted
  Rep = Rep[Rep$first_start > 5781 & Rep$first_start < 16569 & Rep$second_start > 5781 & Rep$second_start < 16569,] 

##### 3 Shift from the start of each arm to the center of each arm, create DirRep and InvRep WHY????? TO DELETE IT!!!
# common repeat will be in DirRep: 8116	8474.5	13451.5	15 : 8116	8467	13444	15	1	CTACCTCCCTCACCA	CAACCTCCCTCACCA	CAACCTCCCTCACCA	1
# closest to him inverted repeat in InvRep: 27341	8488.0	13358.0	10	4: 27341	8483	13353	10	4	AGCCCATAAA	TTTATGTGCT	AGCACATAAA	1  
  table(Rep$id_type)
  DirRep=Rep[Rep$id_type == 1,]
  ComplRep=Rep[Rep$id_type == 2,]
  MirrRep=Rep[Rep$id_type == 3,]
  InvRep=Rep[Rep$id_type == 4,]
  
  # 1- прямые, 2 - комплиментарные, 3 - зеркальные , 4 - инвертированные
  
  par(mfrow=c(4,1))
  hist(InvRep$first_start)
  hist(InvRep$second_start)
  hist(DirRep$first_start)
  hist(DirRep$second_start)
  
  par(mfrow=c(1,1))
  
  #DirRep$first_start = DirRep$first_start + DirRep$length/2
  #InvRep$second_start = InvRep$second_start + InvRep$length/2

##### 3.5: sort all IndInd repeats by length and see if some of the strongest can shape "CONTACT ZONE"
  
summary(InvRep$length)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.00   10.00   10.00   10.47   11.00   16.00 
InvRepSorted = InvRep[order(-InvRep$length),] 
InvRepSorted$Distance =  InvRepSorted$second_start - InvRepSorted$first_start
summary(InvRepSorted$Distance)
PotentSupplTableForPaper = InvRepSorted[InvRepSorted$length >= 14,]
PotentSupplTableForPaper = PotentSupplTableForPaper[order(-PotentSupplTableForPaper$Distance),] # first two are good => make them bold
  
##### 4: find closest InvInv for each DirDir: form DIID combinations
  #Distance2 = 5 # 5
  #DirRep$second_end = DirRep$second_start + DirRep$length
  #DirRep$PresenceOfDeletion = 0
  DirRep$ClosestInvInvAfterFirstDir.Start = 0
  DirRep$ClosestInvInvAfterFirstDir.End = 0
  DirRep$ClosestInvInvAfterFirstDir.FirstGap = 0
  DirRep$ClosestInvInvAfterFirstDir.SecondGap = 0

  DirRep$ClosestInvInvBeforeSecondDir.Start = 0
  DirRep$ClosestInvInvBeforeSecondDir.End = 0
  DirRep$ClosestInvInvBeforeSecondDir.FirstGap = 0
  DirRep$ClosestInvInvBeforeSecondDir.SecondGap = 0
  
  # DirRep$ClosestInvInvBeforeSecondDir = 0
  for (i in 1:nrow(DirRep))
  { # i=1
    TempInvRep = InvRep[InvRep$first_start > (DirRep$first_start[i] + DirRep$length[i]) & (InvRep$second_start + InvRep$length) < DirRep$second_start[i],]
    TempInvRep$FirstGap  = TempInvRep$first_start - (DirRep$first_start[i] + DirRep$length[i]) 
    TempInvRep$SecondGap  = DirRep$second_start[i] - (TempInvRep$second_start + TempInvRep$length) 
    TempInvRep = TempInvRep[order(TempInvRep$FirstGap),]
    DirRep$ClosestInvInvAfterFirstDir.Start[i] = TempInvRep$first_start[1]
    DirRep$ClosestInvInvAfterFirstDir.End[i] = TempInvRep$second_start[1]
    DirRep$ClosestInvInvAfterFirstDir.FirstGap[i] = TempInvRep$FirstGap[1]
    DirRep$ClosestInvInvAfterFirstDir.SecondGap[i] = TempInvRep$SecondGap[1]
    
    TempInvRep = InvRep[InvRep$first_start > (DirRep$first_start[i] + DirRep$length[i]) & (InvRep$second_start + InvRep$length) < DirRep$second_start[i],]
    TempInvRep$FirstGap  = TempInvRep$first_start - (DirRep$first_start[i] + DirRep$length[i]) 
    TempInvRep$SecondGap  = DirRep$second_start[i] - (TempInvRep$second_start + TempInvRep$length) 
    TempInvRep = TempInvRep[order(TempInvRep$SecondGap),]
    
    DirRep$ClosestInvInvBeforeSecondDir.Start[i] = TempInvRep$first_start[1]
    DirRep$ClosestInvInvBeforeSecondDir.End[i]   = TempInvRep$second_start[1]
    DirRep$ClosestInvInvBeforeSecondDir.FirstGap[i] = TempInvRep$FirstGap[1]
    DirRep$ClosestInvInvBeforeSecondDir.SecondGap[i] = TempInvRep$SecondGap[1]
  }
    
summary(DirRep$ClosestInvInvAfterFirstDir.FirstGap) ## NA = 207!!!! why???
summary(DirRep$ClosestInvInvAfterFirstDir.SecondGap) ## NA = 207!!!! why???

summary(DirRep$ClosestInvInvBeforeSecondDir.FirstGap) ## NA = 207!!!! why???
summary(DirRep$ClosestInvInvBeforeSecondDir.SecondGap) ## NA = 207!!!! why???

## NA - because arms of the direct repeats are too close to each other - we can delete such rows from our dataset
NaDirRep = DirRep[is.na(DirRep$ClosestInvInvBeforeSecondDir.Start),]
nrow(NaDirRep) # 207
nrow(DirRep) # 2957
DirRep = DirRep[!is.na(DirRep$ClosestInvInvBeforeSecondDir.Start),] 
nrow(DirRep) # 2750

##### 4.5 how many combinations of Inv1 and Inv2 are overlapped? and how many among them when first Inv pair equals to the second?
nrow(DirRep) # 2750  = 1235 + 1510 + 5 

## two pairs of inverted repeats without overlap:
nrow(DirRep[DirRep$ClosestInvInvBeforeSecondDir.Start > DirRep$ClosestInvInvAfterFirstDir.End,]) # 1235

## two pairs of inverted repeats with overlap:
nrow(DirRep[DirRep$ClosestInvInvBeforeSecondDir.Start < DirRep$ClosestInvInvAfterFirstDir.End,]) # 1510

## end of the first equals start of the second
nrow(DirRep[DirRep$ClosestInvInvBeforeSecondDir.Start == DirRep$ClosestInvInvAfterFirstDir.End,]) # 5
Wau = DirRep[DirRep$ClosestInvInvBeforeSecondDir.Start == DirRep$ClosestInvInvAfterFirstDir.End,]

##### 4.6 derive Intermediate distance between the FirstGap and the SecondGap and flag "OverlapFlag"

DirRep$IntermediateGap = 0 #
DirRep$OverlapFlag = 0 #

for (i in 1:nrow(DirRep))
{
  if (DirRep$ClosestInvInvBeforeSecondDir.Start[i] >= DirRep$ClosestInvInvAfterFirstDir.End[i]) #  no overlap
  { 
    DirRep$IntermediateGap[i] = DirRep$ClosestInvInvBeforeSecondDir.Start[i] -  DirRep$ClosestInvInvAfterFirstDir.End[i] 
  }
  if (DirRep$ClosestInvInvBeforeSecondDir.Start[i] < DirRep$ClosestInvInvAfterFirstDir.End[i]) #  no overlap
  {
    DirRep$IntermediateGap[i] = min((DirRep$ClosestInvInvAfterFirstDir.End[i] - DirRep$ClosestInvInvBeforeSecondDir.Start[i]),(DirRep$ClosestInvInvBeforeSecondDir.End[i] - DirRep$ClosestInvInvAfterFirstDir.End[i]))
    DirRep$OverlapFlag[i] = 1
  }
}
summary(DirRep$IntermediateGap)
DirRep$EffectiveDistance = DirRep$IntermediateGap + DirRep$ClosestInvInvBeforeSecondDir.SecondGap
summary(DirRep$EffectiveDistance)

## 5: derive control zone - make it as in falkenberg!

DirRep$ContactZone = 0
for (i in 1:nrow(DirRep))
{ # i=1
  if (DirRep$first_start[i] >= 6000 &  DirRep$first_start[i] <= 9000 & DirRep$second_start[i] >= 13000 &  DirRep$second_start[i] <= 16000)
  {DirRep$ContactZone[i] = 1}
}
table(DirRep$ContactZone)

## 6: derive PresenceOfDeletion for each direct repeat:
Distance2 = 20 # 5+loop
DirRep$PresenceOfDeletion = 0
DirRep$PresenceOfDeletionDummy = 0
for (i in 1:nrow(DirRep))
{ # i=1
  first_start = DirRep$first_start[i] 
  second_start = DirRep$second_start[i]
  # breaks with close (+/- Distance) coordinates
  TempBreaks = breaks[(breaks$X5..breakpoint > (first_start-Distance2) & breaks$X5..breakpoint < (first_start+Distance2)) & (breaks$X3..breakpoint > (second_start-Distance2) & breaks$X3..breakpoint < (second_start+Distance2)),] 
  DirRep$PresenceOfDeletion[i] = nrow(TempBreaks)
  if (nrow(TempBreaks)>0) {DirRep$PresenceOfDeletionDummy[i] = 1}
}
table(DirRep$PresenceOfDeletion) 
table(DirRep$PresenceOfDeletionDummy) 

DirRep$DistanceBetweenDD = DirRep$second_start - DirRep$first_start
summary(DirRep$DistanceBetweenDD)
summary(DirRep$ClosestInvInvAfterFirstDir.FirstGap)# 53-8
summary(DirRep$ClosestInvInvBeforeSecondDir.SecondGap)# 54-10

## how to understand - which window to take to associate repeats with deletions? run opposite test and associate each deletion with repeats => majority of deletions should be associated with repeats.

### 7: statistics
# A 
# inverted should be close to the first Dir
summary(glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvInvAfterFirstDir.FirstGap*DirRep$ClosestInvInvAfterFirstDir.SecondGap, family = 'poisson'))
summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap), family = 'poisson'))

# inverted should be close to the second Dir
summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvInvBeforeSecondDir.FirstGap) + scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap), family = 'poisson'))

# both pairs of inverted repeats - together
summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap) + scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap), family = 'poisson'))
#  PAPER
#  (Intercept)                                          -2.54234    0.09333 -27.241  < 2e-16 ***
#  scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap)    -1.04850    0.18254  -5.744 9.25e-09 ***
#  scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap) -0.54681    0.11693  -4.676 2.92e-06 ***

# both pairs of inverted repeats - together: FirstGap, SecondGap + other traits: OverlapFlag, IntermediateGap, EffectiveDistance, ContactZone
summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap) + scale(DirRep$OverlapFlag) +scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap), family = 'poisson'))

summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$IntermediateGap) * scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap) + scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap), family = 'poisson'))

summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap) + scale(DirRep$EffectiveDistance), family = 'poisson'))

### effective distance in the middle is difficult to estimate => density of inverted repeats?

summary(glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap) + scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap) + scale(DirRep$ContactZone), family = 'poisson'))
summary(glm(DirRep$PresenceOfDeletion ~ (scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap) + scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap))*scale(DirRep$ContactZone), family = 'poisson'))
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                                                    -2.57141    0.10135 -25.371  < 2e-16 ***
# scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap)                              -0.66842    0.19769  -3.381 0.000722 ***
# scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap)                           -0.35862    0.13139  -2.729 0.006345 ** 
# scale(DirRep$ContactZone)                                                       0.88128    0.08652  10.185  < 2e-16 ***
# scale(DirRep$ClosestInvInvAfterFirstDir.FirstGap):scale(DirRep$ContactZone)     0.54459    0.16860   3.230 0.001237 ** 
# scale(DirRep$ClosestInvInvBeforeSecondDir.SecondGap):scale(DirRep$ContactZone)  0.32634    0.10320   3.162 0.001566 ** 

# In contact zone Gaps are shorter:
wilcox.test(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvAfterFirstDir.FirstGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvInvAfterFirstDir.FirstGap)
t.test(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvAfterFirstDir.FirstGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvInvAfterFirstDir.FirstGap)
boxplot(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvAfterFirstDir.FirstGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvInvAfterFirstDir.FirstGap, outline = FALSE, notch = TRUE, ylab = 'FirstGap', names = c('OutContZone','InContZone'))
length(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvAfterFirstDir.FirstGap)
length(DirRep[DirRep$ContactZone == 1,]$ClosestInvInvAfterFirstDir.FirstGap)
# violinplot - alinaplot

wilcox.test(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvBeforeSecondDir.SecondGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvInvBeforeSecondDir.SecondGap)
t.test(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvBeforeSecondDir.SecondGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvInvBeforeSecondDir.SecondGap)
boxplot(DirRep[DirRep$ContactZone == 0,]$ClosestInvInvBeforeSecondDir.SecondGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvInvBeforeSecondDir.SecondGap, outline = FALSE, notch = TRUE, ylab = 'SecondGap', names = c('OutContZone','InContZone'))
# violinplot - alinaplot


# why? density of inverted repeats is higher?
#
plot(InvRep$first_start,InvRep$second_start, main = 'InvRep')
library(MASS)

kde = kde2d(InvRep$first_start,InvRep$second_start) # , h, n = 25, lims = c(range(x), range(y)))
image(kde,xlim=c(6000,16500), ylim=c(6000,16500), main = 'InvRep')

kde = kde2d(DirRep$first_start,DirRep$second_start) # , h, n = 25, lims = c(range(x), range(y)))
image(kde,xlim=c(6000,16500), ylim=c(6000,16500), main = 'DirRep')

kde = kde2d(ComplRep$first_start,ComplRep$second_start) # , h, n = 25, lims = c(range(x), range(y)))
image(kde,xlim=c(6000,16500), ylim=c(6000,16500), main = 'ComplRep')

kde = kde2d(MirrRep$first_start,MirrRep$second_start) # , h, n = 25, lims = c(range(x), range(y)))
image(kde,xlim=c(6000,16500), ylim=c(6000,16500), main = 'MirrRep')

### extra 1. DIID can be extended to DIIIID

### extra 2. CONTACT ZONE AS AN OVERLAP OF DISTRIBUTION OF DIRECT AND INVERTED REPEATS

### extra 3. PotentSupplTableForPaper

dev.off()

