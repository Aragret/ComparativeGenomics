rm(list=ls(all=TRUE)) # 19 May 2020

pdf("../../Body/4Figures/MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.DiidForEachDD.R01.pdf")
  
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
  Rep = Rep[Rep$id_type == 1 | Rep$id_type == 4,] # 1 = direct, 4 = inverted
  Rep = Rep[Rep$first_start > 5781 & Rep$first_start < 16569 & Rep$second_start > 5781 & Rep$second_start < 16569,] 

##### 3 Shift from the start of each arm to the center of each arm, create DirRep and InvRep WHY????? TO DELETE IT!!!
# common repeat will be in DirRep: 8116	8474.5	13451.5	15 : 8116	8467	13444	15	1	CTACCTCCCTCACCA	CAACCTCCCTCACCA	CAACCTCCCTCACCA	1
# closest to him inverted repeat in InvRep: 27341	8488.0	13358.0	10	4: 27341	8483	13353	10	4	AGCCCATAAA	TTTATGTGCT	AGCACATAAA	1  
  DirRep=Rep[Rep$id_type == 1,]
  InvRep=Rep[Rep$id_type == 4,]
  
  plot(InvRep$first_start,InvRep$second_start)
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
InvRepSorted = InvRep[order(InvRep$length),] 
  

## 4: find closest InvInv for each DirDir: form DIID combinations
  #Distance2 = 5 # 5
  #DirRep$second_end = DirRep$second_start + DirRep$length
  #DirRep$PresenceOfDeletion = 0
  DirRep$ClosestInvFirstStart = 0
  DirRep$ClosestInvSecondStart = 0
  DirRep$ClosestInvFirstGap = 0
  DirRep$ClosestInvSecondGap = 0
  DirRep$ClosestInvTwoGapsTogether = 0
  DirRep$NumberOfTheClosestIi = 0
  DirRep$ClosestInvFirstGap.FirstShorter = 0
  DirRep$ClosestInvSecondGap.FirstShorter = 0
  DirRep$ClosestInvFirstGap.SecondShorter = 0
  DirRep$ClosestInvSecondGap.SecondShorter = 0
  for (i in 1:nrow(DirRep))
  { # i=1
    TempInvRep = InvRep[InvRep$first_start > (DirRep$first_start[i] + DirRep$length[i]) & (InvRep$second_start + InvRep$length) < DirRep$second_start[i],]
    TempInvRep$FirstGap  = TempInvRep$first_start - (DirRep$first_start[i] + DirRep$length[i]) 
    TempInvRep$SecondGap = DirRep$second_start[i] - (TempInvRep$second_start + TempInvRep$length) 
    TempInvRep$TwoGapsTogether = TempInvRep$FirstGap + TempInvRep$SecondGap
    TempInvRep = TempInvRep[order(TempInvRep$TwoGapsTogether),] # the first match will have the minimal "TwoGapsTogether"
    DirRep$ClosestInvFirstStart[i]  = TempInvRep$first_start[1]
    DirRep$ClosestInvSecondStart[i] = TempInvRep$second_start[1]
    DirRep$ClosestInvFirstGap[i]    = TempInvRep$FirstGap[1]
    DirRep$ClosestInvSecondGap[i]    = TempInvRep$SecondGap[1]
    DirRep$ClosestInvTwoGapsTogether[i]    = TempInvRep$TwoGapsTogether[1]
    DirRep$NumberOfTheClosestIi[i] = nrow(TempInvRep[TempInvRep$TwoGapsTogether < 1000,])
    # now sort so, that among the minimal "TwoGapsTogether" we take only those where FirstGap < SecondGap:
    TempInvRep.FirstShorter = TempInvRep[TempInvRep$FirstGap < TempInvRep$SecondGap,]
    DirRep$ClosestInvFirstGap.FirstShorter[i]    = TempInvRep.FirstShorter$FirstGap[1]
    DirRep$ClosestInvSecondGap.FirstShorter[i]    = TempInvRep.FirstShorter$SecondGap[1]
    # now sort so, that among the minimal "TwoGapsTogether" we take only those where FirstGap > SecondGap:
    TempInvRep.SecondShorter = TempInvRep[TempInvRep$FirstGap > TempInvRep$SecondGap,]
    DirRep$ClosestInvFirstGap.SecondShorter[i]    = TempInvRep.SecondShorter$FirstGap[1]
    DirRep$ClosestInvSecondGap.SecondShorter[i]    = TempInvRep.SecondShorter$SecondGap[1]
  }
  
summary(DirRep$ClosestInvFirstGap)
summary(DirRep$ClosestInvSecondGap)
summary(DirRep$ClosestInvTwoGapsTogether) # this is expected to be the most naive metric

## 5: derive control zone - make it as in falkenberg!

DirRep$ContactZone = 0
for (i in 1:nrow(DirRep))
{ # i=1
  if (DirRep$first_start[i] >= 6000 &  DirRep$first_start[i] <= 9000 & DirRep$second_start[i] >= 13000 &  DirRep$second_start[i] <= 16000)
  {DirRep$ContactZone[i] = 1}
}
table(DirRep$ContactZone)

## 6: derive PresenceOfDeletion for each direct repeat:
Distance2 = 40 # 5+loop
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

DirRep$DistanceBetweenDD = DirRep$second_start - DirRep$first_start
summary(DirRep$DistanceBetweenDD)

DirRep = DirRep[!is.na(DirRep$ClosestInvFirstGap),]
## how to understand - which window to take? run opposite test and associate each deletion with repeats => majority of deletions should be associated with repeats.

DirRep = DirRep[DirRep$DistanceBetweenDD > 2000,]

### 7: statistics
### if the minimal sum of two gaps is important for deletions?

cor.test(DirRep$PresenceOfDeletion,DirRep$ClosestInvTwoGapsTogether, method = 'spearman')
plot(DirRep$PresenceOfDeletion,DirRep$ClosestInvTwoGapsTogether, col = rgb(0.1,0.1,0.1,0.5)) # dev.off()
plot(DirRep[DirRep$ContactZone == 1,]$PresenceOfDeletion,DirRep[DirRep$ContactZone == 1,]$ClosestInvTwoGapsTogether, col = rgb(0.1,0.1,0.1,0.5)) # dev.off()



summary(DirRep[DirRep$ClosestInvTwoGapsTogether < 100,]$PresenceOfDeletion)
summary(DirRep[DirRep$ClosestInvTwoGapsTogether < 200,]$PresenceOfDeletion)
summary(DirRep[DirRep$ClosestInvTwoGapsTogether < 300,]$PresenceOfDeletion)
summary(DirRep[DirRep$ClosestInvTwoGapsTogether < 500,]$PresenceOfDeletion)
summary(DirRep[DirRep$ClosestInvTwoGapsTogether < 1000,]$PresenceOfDeletion)
summary(DirRep[DirRep$ClosestInvTwoGapsTogether < 1500,]$PresenceOfDeletion)

(DirRep$ClosestInvFirstGap<)$PresenceOfDeletion

summary(glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvTwoGapsTogether, family = poisson))
summary(glm(DirRep[DirRep$ContactZone == 1,]$PresenceOfDeletion ~ DirRep[DirRep$ContactZone == 1,]$ClosestInvTwoGapsTogether, family = poisson))
summary(glm(DirRep[DirRep$ContactZone == 1,]$PresenceOfDeletion ~ DirRep[DirRep$ContactZone == 1,]$ClosestInvFirstGap, family = poisson))


# A Logistic regression - a probability to have a deletion as a function of FirstGap, SecondGap and NumberOfTheClosestIi

a = glm(DirRep$PresenceOfDeletionDummy ~ DirRep$ClosestInvFirstGap + DirRep$ClosestInvSecondGap + DirRep$NumberOfTheClosestIi + DirRep$DistanceBetweenDD, family = binomial)
summary(a)
a = glm(DirRep$PresenceOfDeletionDummy ~ DirRep$ClosestInvFirstGap + DirRep$NumberOfTheClosestIi + DirRep$DistanceBetweenDD, family = binomial)
summary(a)
a = glm(DirRep$PresenceOfDeletionDummy ~  DirRep$NumberOfTheClosestIi + DirRep$DistanceBetweenDD, family = binomial)
summary(a)

boxplot(DirRep[DirRep$ContactZone == 0,]$NumberOfTheClosestIi,DirRep[DirRep$ContactZone == 1,]$NumberOfTheClosestIi, names = c('OutCZ','InCZ'), ylab = 'NumberOfTheClosestIi', notch = TRUE)
wilcox.test(DirRep[DirRep$ContactZone == 0,]$NumberOfTheClosestIi,DirRep[DirRep$ContactZone == 1,]$NumberOfTheClosestIi) # nothing

# B Poisson regression - number of deletions as a function of ClosestInvTwoGapsTogether (DirRep$ClosestInvFirstGap + DirRep$ClosestInvSecondGap) and DirRep$NumberOfTheClosestIi 
a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvTwoGapsTogether + DirRep$NumberOfTheClosestIi + DirRep$DistanceBetweenDD, family = "poisson")
summary(a)
a = glm(DirRep$PresenceOfDeletion ~ DirRep$NumberOfTheClosestIi + DirRep$DistanceBetweenDD, family = "poisson")
summary(a)

a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap + DirRep$NumberOfTheClosestIi + DirRep$DistanceBetweenDD, family = "poisson")
summary(a)

a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap + DirRep$ClosestInvSecondGap, family = "poisson")
summary(a)
a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap, family = "poisson")
summary(a)
# Good results if distance2 = 20. When it is 10 or 30 - results are qualitatively similar

boxplot(DirRep[DirRep$PresenceOfDeletion == 0,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 1,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 2,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 3,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 4,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion > 4,]$ClosestInvFirstGap, outline = FALSE, names = c("0","1","2","3","4","5+"))
boxplot(DirRep[DirRep$PresenceOfDeletion == 0,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 1,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 2,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 3,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion > 3,]$ClosestInvFirstGap, outline = FALSE, names = c("0","1","2","3","4+"))
boxplot(DirRep[DirRep$PresenceOfDeletion == 0 & DirRep$ContactZone == 1,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 1 & DirRep$ContactZone == 1,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 2 & DirRep$ContactZone == 1,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 3 & DirRep$ContactZone == 1,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion > 3 & DirRep$ContactZone == 1,]$ClosestInvFirstGap, outline = FALSE, names = c("0","1","2","3","4+"), main = 'InCz')
boxplot(DirRep[DirRep$PresenceOfDeletion == 0 & DirRep$ContactZone == 0,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 1 & DirRep$ContactZone == 0,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 2 & DirRep$ContactZone == 0,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion == 3 & DirRep$ContactZone == 0,]$ClosestInvFirstGap,DirRep[DirRep$PresenceOfDeletion > 3 & DirRep$ContactZone == 0,]$ClosestInvFirstGap, outline = FALSE, names = c("0","1","2","3","4+"), main = 'OutCz')

## C
a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap.FirstShorter + DirRep$ClosestInvSecondGap.FirstShorter, family = "poisson")
summary(a) # the longer the second the more deletions
a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap.SecondShorter + DirRep$ClosestInvSecondGap.SecondShorter, family = "poisson")
summary(a)

a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap + DirRep$DistanceBetweenDD + DirRep$ContactZone, family = "poisson")
summary(a)

## final 
a = glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvFirstGap)*scale(DirRep$ContactZone) + scale(DirRep$DistanceBetweenDD), family = "poisson")
summary(a)

a = glm(DirRep$PresenceOfDeletion ~ DirRep$ClosestInvFirstGap*DirRep$DistanceBetweenDD, family = "poisson")
summary(a)

## DirRep$ClosestInvFirstGap is shorter in contact zone! why it has positive coefficient in poisson regression?
# may be this is because close to each other DirDir have also very close InvInv => we need to delete them

plot(DirRep$DistanceBetweenDD,DirRep$ClosestInvFirstGap)
cor.test(DirRep$DistanceBetweenDD,DirRep$ClosestInvFirstGap, method = 'spearman') # nothing

par(mfrow=c(2,2))

boxplot(DirRep[DirRep$ContactZone == 0,]$ClosestInvFirstGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvFirstGap, names = c('OutCZ','InCZ'), ylab = 'FirstGap', notch = TRUE)
wilcox.test(DirRep[DirRep$ContactZone == 0,]$ClosestInvFirstGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvFirstGap)
boxplot(DirRep[DirRep$ContactZone == 0,]$ClosestInvSecondGap,DirRep[DirRep$ContactZone == 1,]$ClosestInvSecondGap, names = c('OutCZ','InCZ'), ylab = 'SecondGap', notch = TRUE)
boxplot(DirRep[DirRep$ContactZone == 0,]$ClosestInvTwoGapsTogether,DirRep[DirRep$ContactZone == 1,]$ClosestInvTwoGapsTogether, names = c('OutCZ','InCZ'), ylab = 'TwoGaps', notch = TRUE)

### Contact zone only: 
DirRep = DirRep[DirRep$ContactZone == 1,] # just 490, fuck


a = glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvFirstGap) + scale(DirRep$ClosestInvSecondGap), family = "poisson")
summary(a)
a = glm(DirRep$PresenceOfDeletion ~ scale(DirRep$ClosestInvTwoGapsTogether), family = "poisson")
summary(a)



dev.off()

