rm(list=ls(all=TRUE)) # 19 May 2020

# pdf("../../Body/4Figures/MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.DiidForEachDD.R01.pdf")
  
##### 1: READ MITOBREAK AND KEEP ONLY MAJOR ARC DELETIONS:
breaks = read.table("../../Body/1Raw/MitoBreakDB_12122019.csv", sep = ',', header = TRUE)
breaks$X5..breakpoint = as.numeric(as.character(breaks$X5..breakpoint)); summary(breaks$X5..breakpoint)
breaks$X3..breakpoint = as.numeric(as.character(breaks$X3..breakpoint)); summary(breaks$X3..breakpoint)
breaks = breaks[!is.na(breaks$X3..breakpoint) & !is.na(breaks$X5..breakpoint),]
  
# although formally major arc is from 5782 through 16569 and till 110 we simplify it for now and take from 5782 till 16569
# OH: 110-441
# OL: 5721-5781
nrow(breaks); breaks = breaks[breaks$X5..breakpoint > 5781 & breaks$X3..breakpoint > 5781,]; nrow(breaks)
summary(breaks$X5..breakpoint)
summary(breaks$X3..breakpoint)

##### 2: read all degraded repeats of Victor (direct):
Rep = read.table("../../Body/1Raw/DegRepsHomo_sapiens 2018_07_02.csv", header = TRUE, sep = ';') 
Rep = Rep[Rep$id_type == 1,] # 1 = direct, 4 = inverted
nrow(Rep) # there are exactly 6304 files in Body/1Raw/VictorFoldings/ folders
# Rep$id is the file name in Body/1Raw/VictorFoldings/
Rep = Rep[Rep$first_start > 5781 & Rep$first_start < 16569 & Rep$second_start > 5781 & Rep$second_start < 16569,] 

DirRep = Rep

##### if not only direct repeats are important - just region of similarity.... the unit if analysis might be each 10 bp window? 
##### may be it is very restrictive to start from direct repeats? some sticky regions might be without direct repeats?
##### find them!!!???

##### 3: derive control zone
DirRep$ContactZone = 0
for (i in 1:nrow(DirRep))
{ # i=1
  if (DirRep$first_start[i] >= 6000 &  DirRep$first_start[i] <= 9000 & DirRep$second_start[i] >= 13000 &  DirRep$second_start[i] <= 16000)
  {DirRep$ContactZone[i] = 1}
}
table(DirRep$ContactZone)

##### 4: derive PresenceOfDeletion for each direct repeat:
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

DirRep$DistanceBetweenDD = DirRep$second_start - DirRep$first_start
summary(DirRep$DistanceBetweenDD)

### how to understand - which window to take? run opposite test and associate each deletion with repeats => majority of deletions should be associated with repeats.

# DirRep = DirRep[DirRep$DistanceBetweenDD > 2000,]

##### 5: read folding energies for each DirRep from Body/1Raw/VictorFoldings/...

# InputDir = "../../Body/1Raw/VictorFoldings/Dir1.G10.50bp.x10.50bp.C10.Dir2/result/"
# InputDir = "../../Body/1Raw/VictorFoldings/G10.Dir1.50bp.x10.50bp.Dir2.C10/result/"
InputDir = "../../Body/1Raw/VictorFoldings/Dir1.300bp.x10.300bp.Dir2/result/"
VecOfFiles = list.files(path = InputDir); length(VecOfFiles) 
Foldings = data.frame()
for (i in 1:length(VecOfFiles))
{ # i = 1   VecOfFiles[559]
  if (file.info(paste(InputDir,VecOfFiles[i],sep = ''))$size > 0)
    {
    Gibbs = read.table(paste(InputDir,VecOfFiles[i],sep = ''),sep = '\t')
    GibbsEnergy = gsub(".*\\(",'',Gibbs$V1[3]); GibbsEnergy = gsub("\\)",'',GibbsEnergy); 
    id = as.numeric(gsub(".csv",'',VecOfFiles[i]))
    OneLine = data.frame(id,as.numeric(GibbsEnergy))
    Foldings = rbind(Foldings,OneLine)
    }
}
names(Foldings)=c('id','GibbsEnergy')
summary(Foldings$GibbsEnergy)
Foldings$AbsGibbsEnergy = -1*Foldings$GibbsEnergy

###### 6 merge DirRep with Foldings
DirRep = merge(DirRep,Foldings, by = 'id')
summary(DirRep$AbsGibbsEnergy)

###### 7 statistics

# Poisson regression
summary(DirRep$PresenceOfDeletion)
table(DirRep$PresenceOfDeletion)
summary(glm(DirRep$PresenceOfDeletion ~ DirRep$AbsGibbsEnergy, family = "poisson"))

# NEXT!!!!!!!!!     DO NOT STICK TO DIRECT REPEATS ONLY - WALK THE WHOLE MAJOR ARC WITH 10 BP STEPS AND ESTIMATE  InvInv Gibbs Energy (probability of replication pausing)
# (16569 - 5782)/10 1000*1000/2
# NEXT!!!!!!!!!     DO NOT STICK TO DIRECT REPEATS ONLY - WALK THE WHOLE MAJOR ARC WITH 10 BP STEPS AND ESTIMATE  homology of direct (probability to realign)
