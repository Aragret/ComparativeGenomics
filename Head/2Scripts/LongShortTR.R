rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

library(gdata)

GenLength <- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header = TRUE, sep='\t')
tr = read.table('../../Body/2Derived/TRFinder.txt', sep='\t', header = TRUE)

GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

tr$FullLength = tr$End - tr$Start

tr$ConsensusLength = as.numeric(lapply(as.character(tr$Consensus), nchar))

summary(tr$ConsensusLength)

ShortTr = tr[tr$ConsensusLength < median(tr$ConsensusLength),]
LongTr = tr[tr$ConsensusLength >= median(tr$ConsensusLength),]

# calculate number of repeats

#### Short consensus length

length(unique(ShortTr$Species))
length(unique(LongTr$Species))

ShortRepNumber = c()
for(i in unique(ShortTr$Species)){
  a = nrow(ShortTr[ShortTr$Species == i,])
  ShortRepNumber = rbind(ShortRepNumber, c(as.character(i), a))
}

ShortRepNumber = as.data.frame(ShortRepNumber)
names(ShortRepNumber) = c('Species', 'RepNumber')

summary(ShortRepNumber$RepNumber)

ShortTrGl = merge(ShortRepNumber, GenLength, by = 'Species')

summary(ShortTrGl$RepNumber)

plot(ShortTrGl$RepNumber, ShortTrGl$GenerationLength_d)


#### Long consensus length

LongRepNumber = c()
for(i in unique(LongTr$Species)){
  a = nrow(LongTr[LongTr$Species == i,])
  LongRepNumber = rbind(LongRepNumber, c(as.character(i), a))
}

LongRepNumber = as.data.frame(LongRepNumber)
names(LongRepNumber) = c('Species', 'RepNumber')

summary(LongRepNumber$RepNumber)

LongTrGl = merge(LongRepNumber, GenLength, by = 'Species')

summary(LongTrGl$RepNumber)

plot(LongTrGl$RepNumber, LongTrGl$GenerationLength_d)

length(intersect(ShortTrGl$Species, LongTrGl$Species)) # 175 sp with both long and short repeats

cor.test(as.numeric(ShortTrGl$RepNumber), ShortTrGl$GenerationLength_d, method = 'spearman')
cor.test(as.numeric(LongTrGl$RepNumber), LongTrGl$GenerationLength_d, method = 'spearman')


data = merge(ShortTrGl, LongTrGl, by='Species')

setdiff(ShortTrGl$Species, LongTrGl$Species)
setdiff(LongTrGl$Species, ShortTrGl$Species)

OnlyShort = ShortTrGl[which(as.character(ShortTrGl$Species) %in% setdiff(ShortTrGl$Species, LongTrGl$Species)),]
OnlyLong = LongTrGl[which(as.character(LongTrGl$Species) %in% setdiff(LongTrGl$Species, ShortTrGl$Species)),]

plot(OnlyLong$RepNumber, OnlyLong$GenerationLength_d)
plot(OnlyShort$RepNumber, OnlyShort$GenerationLength_d)

for(i in 1:nrow(OnlyLong)){
  OnlyLong$ShortOrLong[i] = 1
}

for(i in 1:nrow(OnlyShort)){
  OnlyShort$ShortOrLong[i] = 0
}

ShortLongTr = rbind(OnlyShort, OnlyLong)

cor.test(as.numeric(OnlyShort$RepNumber), OnlyShort$GenerationLength_d, method = 'spearman')
cor.test(as.numeric(OnlyLong$RepNumber), OnlyLong$GenerationLength_d, method = 'spearman')


#######################################################################################################
### repeats length

VEC = unique(ShortTr$Species); length(VEC)
for (i in 1:length(VEC))
{  # i = 1
  species = VEC[i];
  TEMP = ShortTr[ShortTr$Species == species,]
  vec_all = c(1); vec_all = vec_all[-1];
  NumberOfTandemRepeats = nrow(TEMP);
  for (j in 1:nrow(TEMP))
  {
    vec = seq(TEMP$Start[j],TEMP$End[j],1);
    vec_all = c(vec_all,vec)
  }
  LengthOfTandemRepeats = length(vec_all); vec_all = unique(vec_all); LengthOfTandemRepeatsWithoutOverlaps = length(vec_all);
  result_line = data.frame(species, NumberOfTandemRepeats, LengthOfTandemRepeats, LengthOfTandemRepeatsWithoutOverlaps)
  if (i == 1) {ShortTrLength = result_line}
  if (i >  1) {ShortTrLength = rbind(ShortTrLength, result_line)}
}

VEC = unique(LongTr$Species); length(VEC)
for (i in 1:length(VEC))
{  # i = 1
  species = VEC[i];
  TEMP = LongTr[LongTr$Species == species,]
  vec_all = c(1); vec_all = vec_all[-1];
  NumberOfTandemRepeats = nrow(TEMP);
  for (j in 1:nrow(TEMP))
  {
    vec = seq(TEMP$Start[j],TEMP$End[j],1);
    vec_all = c(vec_all,vec)
  }
  LengthOfTandemRepeats = length(vec_all); vec_all = unique(vec_all); LengthOfTandemRepeatsWithoutOverlaps = length(vec_all);
  result_line = data.frame(species, NumberOfTandemRepeats, LengthOfTandemRepeats, LengthOfTandemRepeatsWithoutOverlaps)
  if (i == 1) {LongTrLength = result_line}
  if (i >  1) {LongTrLength = rbind(LongTrLength, result_line)}
}

### merge with GL

ShortTrLengthGL = merge(ShortTrLength, GenLength, by.x='species', by.y = 'Species')
LongTrLengthGL = merge(LongTrLength, GenLength, by.x='species', by.y = 'Species')

cor.test(ShortTrLengthGL$LengthOfTandemRepeats, ShortTrLengthGL$GenerationLength_d, method='spearman')
cor.test(LongTrLengthGL$LengthOfTandemRepeats, LongTrLengthGL$GenerationLength_d, method='spearman')

pdf('../../Body/4Figures/LongShortTR.pdf')
plot(LongTrLengthGL$LengthOfTandemRepeats, LongTrLengthGL$GenerationLength_d)
dev.off()

### without overlaps in species

data = merge(ShortTrLengthGL, LongTrLengthGL, by='species')

setdiff(ShortTrLengthGL$species, LongTrLengthGL$species)
setdiff(LongTrLengthGL$species, ShortTrLengthGL$species)

OnlyShortLength = ShortTrLengthGL[which(as.character(ShortTrLengthGL$species) %in% setdiff(ShortTrLengthGL$species, LongTrLengthGL$species)),]
OnlyLongLength = LongTrLengthGL[which(as.character(LongTrLengthGL$species) %in% setdiff(LongTrLengthGL$species, ShortTrLengthGL$species)),]

plot(OnlyLongLength$LengthOfTandemRepeats, OnlyLongLength$GenerationLength_d)
plot(OnlyShortLength$LengthOfTandemRepeats, OnlyShortLength$GenerationLength_d)

for(i in 1:nrow(OnlyLong)){
  OnlyLong$ShortOrLong[i] = 1
}

for(i in 1:nrow(OnlyShort)){
  OnlyShort$ShortOrLong[i] = 0
}

ShortLongTr = rbind(OnlyShort, OnlyLong)

cor.test(as.numeric(OnlyShortLength$LengthOfTandemRepeats), OnlyShortLength$GenerationLength_d, method = 'spearman')
cor.test(as.numeric(OnlyLongLength$LengthOfTandemRepeats), OnlyLongLength$GenerationLength_d, method = 'spearman')

plot(LongTrLengthGL$LengthOfTandemRepeats, LongTrLengthGL$GenerationLength_d)
