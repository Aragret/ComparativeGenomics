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
  hist(breaks$X5..breakpoint, breaks = seq(0, 16700, 100)) # Note: Shapiro-Wilk Normality Test for X5..breakpoint: p-value = < 0.001
  hist(breaks$X3..breakpoint, breaks = seq(0, 16700, 100)) # Note: Shapiro-Wilk Normality Test for X3..breakpoint: p-value = < 0.001
  
  ggstatsplot::gghistostats(
    data = breaks,
    x = X5..breakpoint,
    effsize.type = "d",
    test.value = 11000,
    test.value.size = TRUE,
    bar.measure = "mix",
    centrality.para = "median",
    test.value.line = TRUE,
    normal.curve = TRUE
  )
  ggstatsplot::gghistostats(
    data = breaks,
    x = X3..breakpoint,
    effsize.type = "d",
    test.value = 11000,
    test.value.size = TRUE,
    bar.measure = "mix",
    centrality.para = "median",
    test.value.line = TRUE,
    normal.curve = TRUE
  )
  
## 2: read all degraded repeats of Victor (direct and inverted):
  Rep = read.table("../../Body/1Raw/DegRepsHomo_sapiens 2018_07_02.csv", header = TRUE, sep = ';') 
  Rep = Rep[Rep$id_type == 1 | Rep$id_type == 4,] # 1 = direct, 4 = inverted
  Rep = Rep[Rep$first_start > 5781 & Rep$first_start < 16569 & Rep$second_start > 5781 & Rep$second_start < 16569,] 

## 3 Shift from the start of each arm to the center of each arm, create DirRep and InvRep
# common repeat will be in DirRep: 8116	8474.5	13451.5	15 : 8116	8467	13444	15	1	CTACCTCCCTCACCA	CAACCTCCCTCACCA	CAACCTCCCTCACCA	1
# closest to him inverted repeat in InvRep: 27341	8488.0	13358.0	10	4: 27341	8483	13353	10	4	AGCCCATAAA	TTTATGTGCT	AGCACATAAA	1  
  DirRep=Rep[Rep$id_type == 1,]
  InvRep=Rep[Rep$id_type == 4,]
  
DirRep_keep <- DirRep
InvRep_keep <- InvRep
Rep_keep <- Rep

  
  DirRep$first_start = DirRep$first_start + DirRep$length/2
  InvRep$second_start = InvRep$second_start + InvRep$length/2
  
## 4: derive control zone - make it as in falkenberg!
DirRep$ContactZone = 0
  for (i in 1:nrow(DirRep))
  { # i=1
    if (DirRep$first_start[i] >= 6000 &  DirRep$first_start[i] <= 9000 & DirRep$second_start[i] >= 13000 &  DirRep$second_start[i] <= 16000)
    {DirRep$ContactZone[i] = 1}
  }
  table(DirRep$ContactZone)
  ggstatsplot::ggpiestats(
    data = as.data.frame(DirRep),
    x = ContactZone, # <<
    slice.label = "both", # <<
    messages = TRUE
  )
  
## 5: derive InvertedInside (diid) and InvertedOutside (iddi) combinations for each direct repeat:
# if inverted: Distance1A = 1000 & Distance1B = 10 => 0	1000	5	0.06149588	2.6172263	2696_22_234_5  
# if acording to Falcenberg: Distance1A = 10 & Distance1B = 1000 => 0	10	5	0.007159735	3.7005539	2677_20_253_7  
  loop = 0
  #for (loop in 0:20)
  #{
  Distance1A = 10   # 10
  Distance1B = 1000 # 1000
  DirRep$InvertedInside = 0
  DirRep$InvertedOutside = 0
  for (i in 1:nrow(DirRep))
  {
    # i=1
    first_start = DirRep$first_start[i]
    second_start = DirRep$second_start[i]
    # DirRep$InvertedInside (cases)
    InvTemp = InvRep[(InvRep$first_start > first_start &
                        InvRep$first_start < (first_start + Distance1A)) &
                       (
                         InvRep$second_start < second_start &
                           InvRep$second_start > (second_start - Distance1B)
                       ), ]
    if (nrow(InvTemp) > 0)   {
      DirRep$InvertedInside[i] = 1
    }
    # DirRep$InvertedOutside (controls)
    InvTemp = InvRep[(InvRep$first_start < first_start &
                        InvRep$first_start > (first_start - Distance1A)) &
                       (
                         InvRep$second_start > second_start &
                           InvRep$second_start < (second_start + Distance1B)
                       ), ]
    if (nrow(InvTemp) > 0)   {
      DirRep$InvertedOutside[i] = 1
    }
  }
  # number of DIID combinations is always < than IDDI with different Distance parameter!!!
  # using this difference we can touch the distance of the effect - walk with windows step by step...!
  # if this difference (iddi - diid) is stronger within the contact zone?
  # we can do the same for mammals!!!! and bacteria - longevity!!!
  
  DIID = nrow(DirRep[DirRep$InvertedInside == 1,])  # 260
  IDDI = nrow(DirRep[DirRep$InvertedOutside == 1,]) # 257
  table(DirRep$InvertedInside) # 1798; 976; 634; 296; 77; 24;
  table(DirRep$InvertedOutside) # 1859; 1099; 669; 326; 80; 36;

  DIIDInContactZone = nrow(DirRep[DirRep$ContactZone == 1 & DirRep$InvertedInside == 1,])
  IDDIInContactZone = nrow(DirRep[DirRep$ContactZone == 1 & DirRep$InvertedOutside == 1,])
  table(DirRep[DirRep$ContactZone == 1,]$InvertedInside)  # 4;
  table(DirRep[DirRep$ContactZone == 1,]$InvertedOutside) # 7
  
  DIIDOutOfContactZone = nrow(DirRep[DirRep$ContactZone == 0 & DirRep$InvertedInside == 1,])
  IDDIOutOfContactZone = nrow(DirRep[DirRep$ContactZone == 0 & DirRep$InvertedOutside == 1,])
  table(DirRep[DirRep$ContactZone == 0,]$InvertedInside)  # 25;
  table(DirRep[DirRep$ContactZone == 0,]$InvertedOutside) # 31
  
  X = cbind(c(DIIDInContactZone,IDDIInContactZone),c(DIIDOutOfContactZone,IDDIOutOfContactZone))
  F = fisher.test(X)
  FisherP.DiiD.IddI = F$p.value
  FisherOdds.DiiD.IddI = as.numeric(F$estimate)

## 6: derive PresenceOfDeletion for each direct repeat:
  Distance2 = 5 # 5+loop
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
table(DirRep$PresenceOfDeletion) # 2930/27
DirRep = DirRep[order(-DirRep$PresenceOfDeletion,-DirRep$InvertedInside),]
    
## 6: compare frequencies:
  
mean(DirRep[DirRep$InvertedInside == 1,]$PresenceOfDeletion)  
mean(DirRep[DirRep$InvertedInside == 0,]$PresenceOfDeletion)  
mean(DirRep[DirRep$InvertedOutside == 1,]$PresenceOfDeletion)  
mean(DirRep[DirRep$InvertedOutside == 0,]$PresenceOfDeletion)  

AB=data.frame(table(DirRep[DirRep$InvertedInside == 1,]$PresenceOfDeletion))  
CD=as.data.frame(table(DirRep[DirRep$InvertedInside == 0,]$PresenceOfDeletion))  
X = cbind(CD$Freq,AB$Freq) # (2/27) / (81/2847) 
F = fisher.test(X)
FisherP = F$p.value
FisherOdds = as.numeric(F$estimate)
FisherMatrixCDAB = paste(c(CD$Freq,AB$Freq),collapse = '_')
OneLine1 = data.frame(loop,Distance1A,Distance2,FisherP,FisherOdds,FisherMatrixCDAB,DIID,IDDI,DIIDInContactZone,IDDIInContactZone,DIIDOutOfContactZone,IDDIOutOfContactZone,FisherP.DiiD.IddI,FisherOdds.DiiD.IddI)
if (loop == 0) {Final = OneLine1}
if (loop >  0) {Final = rbind(Final,OneLine1)}
# }

## Fisher test within the contact zone is even stronger:

AB=data.frame(table(DirRep[DirRep$ContactZone == 1 & DirRep$InvertedInside == 1,]$PresenceOfDeletion))  
CD=as.data.frame(table(DirRep[DirRep$ContactZone == 1 & DirRep$InvertedInside == 0,]$PresenceOfDeletion))  
X = cbind(CD$Freq,AB$Freq) # (6/34) / (12/438) ; 
F = fisher.test(X) # 6.393955; 1.9421 p = 0.4012; (3/44) / (15/428) 

# if distance 1 == 50 and distance 2 = 5, DiiD have 0.0259 fraction of deletions but not DiiD - 0.0142
# if distance 1 == 30 and distance 2 = 5, DiiD have 0.04   fraction of deletions but not DiiD - 0.0144

### some proofs of assymetry:  better association with deletions if first distance is short
# 10	5	0.007159735	3.7005539 2677_20_253_7 # distance1A = 10, distance1B = 1000, distance2 = 5; 
X = cbind(c(2677,20),c(253,7)) # 2957 - 260
mosaicplot(X)    
fisher.test(X)

# the same but within the contact zone:

# 1000	5	0.06149588	2.6172263	2696_22_234_5 # distance1A = 1000, distance1B = 10, distance2 = 5; 
X = cbind(c(2696,22),c(234,5))
mosaicplot(X)    
fisher.test(X)

# the same but within the contact zone:


## reshufle direct and inverted repeats (types of repeats) within the major arc and count expected # of DIID
# observed number of DIID is 260 & IDDI is 257 
for (permut in 1:10000) # 10000
{
Perm = Rep
Perm$id_type = sample(Perm$id_type)
Perm$first_start = Perm$first_start + Perm$length/2
Perm$second_start = Perm$second_start + Perm$length/2
DirRep=Perm[Perm$id_type == 1,]
InvRep=Perm[Perm$id_type == 4,]
Distance1A = 10   # 10
Distance1B = 1000 # 1000
DirRep$InvertedInside = 0
DirRep$InvertedOutside = 0
for (i in 1:nrow(DirRep))
  { 
  first_start = DirRep$first_start[i]
  second_start = DirRep$second_start[i]
  InvTemp = InvRep[(InvRep$first_start > first_start & InvRep$first_start < (first_start+Distance1A)) & (InvRep$second_start < second_start & InvRep$second_start > (second_start - Distance1B)),] 
  if (nrow(InvTemp) > 0)   {DirRep$InvertedInside[i] = 1}
  InvTemp = InvRep[(InvRep$first_start < first_start & InvRep$first_start > (first_start-Distance1A)) & (InvRep$second_start > second_start & InvRep$second_start < (second_start + Distance1B)),] 
  if (nrow(InvTemp) > 0)   {DirRep$InvertedOutside[i] = 1}
  }
DIID = nrow(DirRep[DirRep$InvertedInside == 1,])
IDDI = nrow(DirRep[DirRep$InvertedOutside == 1,])
if (permut == 1) {DiidVec = DIID; IddiVec = IDDI}
if (permut > 1 ) {DiidVec = c(DiidVec,DIID); IddiVec = c(IddiVec,IDDI); }
}

summary(DiidVec)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 301.0   354.0   365.0   365.4   377.0   430.0 
summary(IddiVec)

median(DiidVec)
Pvalue = length(DiidVec[DiidVec <= 260])/length(DiidVec) # 0
RelativeDiff = (median(DiidVec) - 260)/median(DiidVec); RelativeDiff # 29% => expected median is 29% higher than observed 
hist(DiidVec, breaks = seq(250,450,1), col = 'gray')
abline(v=260, col = 'red', lwd = 2)

#### within the contact zone:
## reshufle direct and inverted repeats (types of repeats) within the major arc and count expected # of DIID
# observed number of DIID is 40 within the contact zone 

for (permut in 1:10000) # 10000
{
  Perm = Rep
  Perm = Perm[Perm$first_start >= 6000 & Perm$first_start <= 9000 & Perm$second_start >= 13000 & Perm$second_start <= 16000,]
  Perm$id_type = sample(Perm$id_type)
  Perm$first_start = Perm$first_start + Perm$length/2
  Perm$second_start = Perm$second_start + Perm$length/2
  DirRep=Perm[Perm$id_type == 1,]
  InvRep=Perm[Perm$id_type == 4,]
  Distance1A = 10   # 10
  Distance1B = 1000 # 1000
  DirRep$InvertedInside = 0
  DirRep$InvertedOutside = 0
  for (i in 1:nrow(DirRep))
  { 
    first_start = DirRep$first_start[i]
    second_start = DirRep$second_start[i]
    InvTemp = InvRep[(InvRep$first_start > first_start & InvRep$first_start < (first_start+Distance1A)) & (InvRep$second_start < second_start & InvRep$second_start > (second_start - Distance1B)),] 
    if (nrow(InvTemp) > 0)   {DirRep$InvertedInside[i] = 1}
    InvTemp = InvRep[(InvRep$first_start < first_start & InvRep$first_start > (first_start-Distance1A)) & (InvRep$second_start > second_start & InvRep$second_start < (second_start + Distance1B)),] 
    if (nrow(InvTemp) > 0)   {DirRep$InvertedOutside[i] = 1}
  }
  DIID = nrow(DirRep[DirRep$InvertedInside == 1,])
  IDDI = nrow(DirRep[DirRep$InvertedOutside == 1,])
  if (permut == 1) {DiidVec = DIID; IddiVec = IDDI}
  if (permut > 1 ) {DiidVec = c(DiidVec,DIID); IddiVec = c(IddiVec,IDDI); }
}

summary(DiidVec)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 301.0   354.0   365.0   365.4   377.0   430.0 
summary(IddiVec)

median(DiidVec) # 54
Pvalue = length(DiidVec[DiidVec <= 40])/length(DiidVec);  Pvalue # 0.0243
RelativeDiff = (median(DiidVec) - 40)/median(DiidVec); RelativeDiff # 0.26 => expected median is 29% higher than observed 
hist(DiidVec, breaks = seq(20,100,1), col = 'gray')
abline(v=40, col = 'red', lwd = 2)

dev.off()  
