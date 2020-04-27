rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.SearchForBestDiid.R01.pdf")
  
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

## 3 Shift from the start of each arm to the center of each arm, create DirRep and InvRep
# common repeat will be in DirRep: 8116	8474.5	13451.5	15 : 8116	8467	13444	15	1	CTACCTCCCTCACCA	CAACCTCCCTCACCA	CAACCTCCCTCACCA	1
# closest to him inverted repeat in InvRep: 27341	8488.0	13358.0	10	4: 27341	8483	13353	10	4	AGCCCATAAA	TTTATGTGCT	AGCACATAAA	1  
  DirRep=Rep[Rep$id_type == 1,]
  InvRep=Rep[Rep$id_type == 4,]
  
  DirRep$first_start = DirRep$first_start + DirRep$length/2
  InvRep$second_start = InvRep$second_start + InvRep$length/2
  
## 4: define control zone
DirRep$ContactZone = 0
  for (i in 1:nrow(DirRep))
  { # i=1
    if (DirRep$first_start[i] >= 6000 &  DirRep$first_start[i] <= 9000 & DirRep$second_start[i] >= 13000 &  DirRep$second_start[i] <= 16000)
    {DirRep$ContactZone[i] = 1}
  }
  table(DirRep$ContactZone)

## 5: derive PresenceOfDeletion for each direct repeat
  Distance2 = 5 # 5
  DirRep$second_end = DirRep$second_start + DirRep$length
  DirRep$PresenceOfDeletion = 0
  for (i in 1:nrow(DirRep))
  { # i=1
    first_start = DirRep$first_start[i] # I have to take into account the length of the repeat and the fact that one arm is always removed!!
    second_start = DirRep$second_start[i]
    # breaks with close (+/- Distance) coordinates
    TempBreaks = breaks[(breaks$X5..breakpoint > (first_start-Distance2) & breaks$X5..breakpoint < (first_start+Distance2)) & (breaks$X3..breakpoint > (second_start-Distance2) & breaks$X3..breakpoint < (second_start+Distance2)),] 
    if (nrow(TempBreaks) > 0)   {DirRep$PresenceOfDeletion[i] = 1}
  }
  table(DirRep$PresenceOfDeletion) # 2930/27 /109 / 218 
  # 109/(2848+109) = 0.036 = random chance that direct repeat is associated with deletion
  # 218/(218 + 2739) = 0.073 = random chance that direct repeat is associated with deletion
  # 27 / (27 + 2930) = 0.009

## 6: derive DIID as a function of DistanceX, Y and Z: Dir-DistanceX-Inv-DistanceY-Inv-DistanceZ-Dir 
DistanceXVec = c(10,50,seq(100,500,100))
DistanceYVec = seq(2000,9000,500)
DistanceZVec = c(10,50,seq(100,500,100))
FirstLoopFlag = 1
  
for (DistanceX in DistanceXVec) 
  { # DistanceX = 50
  for (DistanceY in DistanceYVec) 
    { # DistanceY = 5000
    for (DistanceZ in DistanceZVec) 
      { # DistanceZ = 100
      Temp = DirRep
      Temp$InvertedInside = 0
      Temp$DiidDel = 0
      Temp$DiidNoDel = 0
      for (i in 1:nrow(Temp))
        { # i=1
        first_start = Temp$first_start[i]
        second_start = Temp$second_start[i]
        InvTemp = InvRep[(InvRep$second_start - InvRep$first_start < DistanceY) & (InvRep$first_start > first_start & InvRep$first_start < (first_start+DistanceX)) & (InvRep$second_start < second_start & InvRep$second_start > (second_start - DistanceZ)),] 
        if (nrow(InvTemp) > 0)   {Temp$InvertedInside[i] = 1}
        if (nrow(InvTemp) > 0 & Temp$PresenceOfDeletion[i] == 1)   {Temp$DiidDel[i] = 1}
        if (nrow(InvTemp) > 0 & Temp$PresenceOfDeletion[i] == 0)   {Temp$DiidNoDel[i] = 1}
        }
      if ((sum(Temp$DiidNoDel)+sum(Temp$DiidDel)) == 0) {FrOfDiidWithDel = 0}
      if ((sum(Temp$DiidNoDel)+sum(Temp$DiidDel)) >  0) 
        {
        FrOfDiidWithDel = sum(Temp$DiidDel)/(sum(Temp$DiidNoDel)+sum(Temp$DiidDel))
        }
      NumDiid = nrow(Temp[Temp$InvertedInside == 1,])
      MedianStart = median(Temp[Temp$InvertedInside == 1,]$first_start)
      MedianEnd = median(Temp[Temp$InvertedInside == 1,]$second_start)
      ResLine = data.frame(DistanceX,DistanceY,DistanceZ,NumDiid,FrOfDiidWithDel,MedianStart,MedianEnd)
      if (FirstLoopFlag == 0) {Final = rbind(Final,ResLine)}
      if (FirstLoopFlag == 1) {Final = ResLine; FirstLoopFlag = 0}
      }
    }
  }

Final = Final[order(-Final$FrOfDiidWithDel),]

# I can also walk with bins (not with cumulative curves) -> it will show if there is indeed trend or - only noise

##### DIIIIIIID => walk along each direct repeat and collect all inverted (with different settings) which fit them


