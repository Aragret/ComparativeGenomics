rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.R01.pdf")
  
## 1: READ MITOBREAK AND KEEP ONLY MAJOR ARC DELETIONS:
  breaks = read.table("../../Body/1Raw/MitoBreakDB_12122019.csv", sep = ',', header = TRUE)
  breaks$X5..breakpoint = as.numeric(as.character(breaks$X5..breakpoint)); summary(breaks$X5..breakpoint)
  breaks$X3..breakpoint = as.numeric(as.character(breaks$X3..breakpoint)); summary(breaks$X3..breakpoint)
  breaks = breaks[!is.na(breaks$X3..breakpoint) & !is.na(breaks$X5..breakpoint),]
  
  par(mfrow=c(2,1))
  hist(breaks$X5..breakpoint, breaks = seq(0, 16600, 100))
  hist(breaks$X3..breakpoint, breaks = seq(0, 16600, 100))
  nrow(breaks); breaks = breaks[breaks$Deletion.of.replication.origins == 'None',]; nrow(breaks)
  breaks = breaks[breaks$Location.of.the.deleted.region == 'Inside the major arc',]; nrow(breaks)
  summary(breaks$X5..breakpoint)
  summary(breaks$X3..breakpoint)
  hist(breaks$X5..breakpoint, breaks = seq(0, 16600, 100))
  hist(breaks$X3..breakpoint, breaks = seq(0, 16600, 100))
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
  hist(breaks$X5..breakpoint, breaks = seq(0, 16700, 100))
  hist(breaks$X3..breakpoint, breaks = seq(0, 16700, 100))
  
## 2: read all degraded repeats of Victor (direct and inverted):
  Rep = read.table("../../Body/1Raw/DegRepsHomo_sapiens 2018_07_02.csv", header = TRUE, sep = ';') 
  Rep = Rep[Rep$id_type == 1 | Rep$id_type == 4,] # 1 = direct, 4 = inverted
  Rep = Rep[Rep$first_start > 5781 & Rep$first_start < 16569 & Rep$second_start > 5781 & Rep$second_start < 16569,] 
  
## 3: derive InvertedInside (diid) and InvertedOutside (iddi) combinations for each direct repeat:
  Distance = 30 # from 30 till 300
  DirRep=Rep[Rep$id_type == 1,]
  InvRep=Rep[Rep$id_type == 4,]
  DirRep$InvertedInside = 0
  DirRep$InvertedOutside = 0
  for (i in 1:nrow(DirRep))
  { # i=1
  first_start = DirRep$first_start[i]
  second_start = DirRep$second_start[i]
  # DirRep$InvertedInside (cases)
  InvTemp = InvRep[(InvRep$first_start > first_start & InvRep$first_start < (first_start+Distance)) & (InvRep$second_start < second_start & InvRep$second_start > (second_start - Distance)),] 
  if (nrow(InvTemp) > 0)   {DirRep$InvertedInside[i] = 1}
  # DirRep$InvertedOutside (controls)
  InvTemp = InvRep[(InvRep$first_start < first_start & InvRep$first_start > (first_start-Distance)) & (InvRep$second_start > second_start & InvRep$second_start < (second_start + Distance)),] 
  if (nrow(InvTemp) > 0)   {DirRep$InvertedOutside[i] = 1}
  }
  # number of DIID combinations is always < than IDDI with different Distance parameter!!!
  # using this difference we can touch the distance of the effect - walk with windows step by step...!
  # if this difference (iddi - diid) is stronger within the contact zone?
  # we can do the same for mammals!!!! and bacteria - longevity!!!
  table(DirRep$InvertedInside) # 1798; 976; 634; 296; 77; 24;
  table(DirRep$InvertedOutside) # 1859; 1099; 669; 326; 80; 36;

## 3: derive PresenceOfDeletion for each direct repeat:
  Distance = 5 # from 30 till 300
  DirRep$PresenceOfDeletion = 0
  for (i in 1:nrow(DirRep))
  { # i=1
    first_start = DirRep$first_start[i] # I have to take into account the length of the repeat and the fact that one arm is always removed!!
    second_start = DirRep$second_start[i]
    # breaks with close (+/- Distance) coordinates
    TempBreaks = breaks[(breaks$X5..breakpoint > (first_start-Distance) & breaks$X5..breakpoint < (first_start+Distance)) & (breaks$X3..breakpoint > (second_start-Distance) & breaks$X3..breakpoint < (second_start+Distance)),] 
    if (nrow(TempBreaks) > 0)   {DirRep$PresenceOfDeletion[i] = 1}
  }
  table(DirRep$PresenceOfDeletion) # 2738/219
  
## 4: compare frequencies:
  
mean(DirRep[DirRep$InvertedInside == 1,]$PresenceOfDeletion)  
mean(DirRep[DirRep$InvertedInside == 0,]$PresenceOfDeletion)  
mean(DirRep[DirRep$InvertedOutside == 1,]$PresenceOfDeletion)  
mean(DirRep[DirRep$InvertedOutside == 0,]$PresenceOfDeletion)  

# if distance 1 == 50 and distance 2 = 5, DiiD have 0.0259 fraction of deletions but not DiiD - 0.0142
# if distance 1 == 30 and distance 2 = 5, DiiD have 0.04   fraction of deletions but not DiiD - 0.0144

  
  


dev.off()  
