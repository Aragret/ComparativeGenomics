rm(list=ls(all=TRUE))

library(ggplot2)
library(ggpubr)
#library(mosaicplot) # install.packages('mosaicplot')

### 1: READ VICTOR's FILE 
Rep = read.table("../../Body/2Derived/Homo_sapiens_triangles.txt", sep = "\t", header = FALSE)
names(Rep)=c('DelStart','DelEnd','DelLength','MasterRepeat','AlternRepeats1','AlternRepeats2','AlternRepeats3','AlternRepeats4')
# MasterRepeat is a repeat with arms located close to given deletion breakpoints
# AlternativeRepeats1-4 are repeats where: first arm == arm1 of Master; first arm == arm2 of Master; second arm == arm1 of Master; second arm == arm2 of Master; 
# mb_del == реализованные делеции (существуют в MitoBreak); non_del - не реализованные == не существуют в MitoBreak

### 2: filter out deletions within major arc:
# поскольку координаты не такие простые (см ниже) - чтобы не париться можно взять все точки разрыва что больше чем 5781 и меньше чем 16569. 
# однако, позже можно подумать и взять во внимание хвостик (< 110) - добавить его к 16569
# OH: 110-441
# OL: 5721-5781
nrow(Rep)
Rep=Rep[Rep$DelStart > 5781 & Rep$DelStart < 16569 & Rep$DelEnd > 5781 & Rep$DelEnd < 16569,]
nrow(Rep)

### 3: concatenate all alternative deletions
Rep$AllAlternRepeats = paste(Rep$AlternRepeats1,Rep$AlternRepeats2,Rep$AlternRepeats3,Rep$AlternRepeats4,sep=',')
Rep$AllAlternRepeats = gsub("\\,\\[\\]\\,",",",Rep$AllAlternRepeats)
Rep$AllAlternRepeats = gsub("^\\[","",Rep$AllAlternRepeats)
Rep$AllAlternRepeats = gsub("\\]$","",Rep$AllAlternRepeats)

### 4: prepare dataset for analysis of realized and nonrealized deletions line by line
Rep$CenterOfRealizedRepeats = 0
Rep$CenterOfNonRealizedRepeats = 0
Rep$LengthOfRealizedRepeats = 0
Rep$LengthOfNonRealizedRepeats = 0
Rep$StartOfRealizedRepeats = 0
Rep$StartOfNonRealizedRepeats = 0
Rep$EndOfRealizedRepeats = 0
Rep$EndOfNonRealizedRepeats = 0

for (i in (1:nrow(Rep)))
{
# i = 1
# format of data to get a dataset of master and all alternative repeats for each deletion (for each line of the dataset)
temp =  Rep[i,]
AltRep = unlist(strsplit(temp$AllAlternRepeats, "\\]\\,\\[")); AltRep = AltRep[AltRep != ""]
AltRep = data.frame(AltRep); names(AltRep)=c('WholeLine')
if (nrow(AltRep) > 0)
  {
  AltRep$RepeatType = 'alternative'
  AltRep$WholeLine = as.character(AltRep$WholeLine)

  MasterRepeat = as.character(temp$MasterRepeat); MasterRepeat = gsub("^\\[","",MasterRepeat); MasterRepeat = gsub("\\]$","",MasterRepeat); 
  MasterRepeat = paste(MasterRepeat,'mb_del',sep=' ')
  MasterRepeat = data.frame(MasterRepeat); names(MasterRepeat)=c('WholeLine'); MasterRepeat$RepeatType = 'master'

  AllRep=rbind(MasterRepeat,AltRep)
  ReturnFifth = function(x)  {unlist(strsplit(x,' '))[5]}; AllRep$RealisedRepeat = apply(as.matrix(AllRep$WholeLine),1,FUN = ReturnFifth)
  ReturnFirst = function(x)  {as.numeric(unlist(strsplit(x,' '))[1])}; AllRep$RepStart = apply(as.matrix(AllRep$WholeLine),1,FUN = ReturnFirst)
  ReturnSecond = function(x)  {as.numeric(unlist(strsplit(x,' '))[2])}; AllRep$RepEnd = apply(as.matrix(AllRep$WholeLine),1,FUN = ReturnSecond)
  AllRep = AllRep[AllRep$RepStart > 5781 & AllRep$RepStart < 16569 & AllRep$RepEnd > 5781 & AllRep$RepEnd < 16569,]
  if (nrow(AllRep) > 0)
    {
    AllRep$Center = (AllRep$RepEnd - AllRep$RepStart)/2 + AllRep$RepStart
    Rep$CenterOfRealizedRepeats[i]  = mean(AllRep[AllRep$RealisedRepeat == 'mb_del',]$Center)
    Rep$CenterOfNonRealizedRepeats[i]= mean(AllRep[AllRep$RealisedRepeat == 'non_del',]$Center)
    AllRep$Length = AllRep$RepEnd - AllRep$RepStart
    Rep$LengthOfRealizedRepeats[i] = mean(AllRep[AllRep$RealisedRepeat == 'mb_del',]$Length)
    Rep$LengthOfNonRealizedRepeats[i] =  mean(AllRep[AllRep$RealisedRepeat == 'non_del',]$Length) 
    Rep$StartOfNonRealizedRepeats[i] =   mean(AllRep[AllRep$RealisedRepeat == 'non_del',]$RepStart)
    Rep$StartOfRealizedRepeats[i] =   mean(AllRep[AllRep$RealisedRepeat == 'mb_del',]$RepStart)
    Rep$EndOfNonRealizedRepeats[i] =   mean(AllRep[AllRep$RealisedRepeat == 'non_del',]$RepEnd)
    Rep$EndOfRealizedRepeats[i] =   mean(AllRep[AllRep$RealisedRepeat == 'mb_del',]$RepEnd)
    
    if (i == 1) {FinalAllRep = AllRep}
    if (i >  1) {FinalAllRep = rbind(FinalAllRep,AllRep)}
    }
  }
}
nrow(Rep) # 703 
Rep=Rep[Rep$CenterOfRealizedRepeats > 0,]; nrow(Rep) # 643
Rep=Rep[Rep$CenterOfNonRealizedRepeats >0,] # it means there is at least one extra (third) arm to analyze
nrow(Rep) # 643
Rep = Rep[!is.na(Rep$LengthOfRealizedRepeat),] # 618

##
pdf("../../Body/4Figures/RealizedVsNonrealizedDeletions.R.01.pdf")

## center is a bit higher in realized repeats
wilcox.test(Rep$CenterOfRealizedRepeats,Rep$CenterOfNonRealizedRepeats,paired = TRUE) # significant
t.test(Rep$CenterOfRealizedRepeats,Rep$CenterOfNonRealizedRepeats,paired = TRUE) # significant
summary(Rep$CenterOfRealizedRepeats)
summary(Rep$CenterOfNonRealizedRepeats)
boxplot(Rep$CenterOfRealizedRepeats,Rep$CenterOfNonRealizedRepeats, notch = TRUE, names=c('CenterOfRealizedRepeats','CenterOfNonRealizedRepeats'))

## length is longer in realized repeats
wilcox.test(Rep$LengthOfRealizedRepeats,Rep$LengthOfNonRealizedRepeats,paired = TRUE) # significant
t.test(Rep$LengthOfRealizedRepeats,Rep$LengthOfNonRealizedRepeats,paired = TRUE) # significant
summary(Rep$LengthOfRealizedRepeats)   # 5643  mean
summary(Rep$LengthOfNonRealizedRepats) # 3450  mean
summary(Rep$LengthOfRealizedRepeats - Rep$LengthOfNonRealizedRepeats) # 1998 mean
boxplot(Rep$LengthOfRealizedRepeats,Rep$LengthOfNonRealizedRepeats, notch = TRUE, names=c('LengthOfRealizedRepeats','LengthOfNonRealizedRepeats'))

## start and end in realized versus non-realized repeats 
wilcox.test(Rep$StartOfRealizedRepeats,Rep$StartOfNonRealizedRepeats,paired = TRUE) # significant
t.test(Rep$StartOfRealizedRepeats,Rep$StartOfNonRealizedRepeats,paired = TRUE) # significant
summary(Rep$StartOfRealizedRepeats) # mean 8678: 25%: 7401 - 75% 9764:
quantile(Rep$StartOfRealizedRepeats, 0.1) # 6465
quantile(Rep$StartOfRealizedRepeats, 0.9) # 10954

summary(Rep$StartOfNonRealizedRepeats) # mean 9394
summary(Rep$StartOfNonRealizedRepeats-Rep$StartOfRealizedRepeats) # mean 716 => non realised start later (~700 bp): 9394 - 8678

wilcox.test(Rep$EndOfRealizedRepeats,Rep$EndOfNonRealizedRepeats,paired = TRUE) # significant
t.test(Rep$EndOfRealizedRepeats,Rep$EndOfNonRealizedRepeats,paired = TRUE) # significant
summary(Rep$EndOfRealizedRepeats)    # 14321 mean
summary(Rep$EndOfNonRealizedRepeats) # 13039 mean
summary(Rep$EndOfNonRealizedRepeats-Rep$EndOfRealizedRepeats) # mean -1281=> non realised end earlier (~1300 bp): 13039 - 14321
quantile(Rep$EndOfRealizedRepeats, 0.1) # 13286
quantile(Rep$EndOfRealizedRepeats, 0.9) # 15863


## how to plot it in terms of X(start) and Y(end)? Just plot it? ALL?

plot(FinalAllRep[FinalAllRep$RealisedRepeat == 'non_del',]$RepStart,FinalAllRep[FinalAllRep$RealisedRepeat == 'non_del',]$RepEnd, pch = 16, col = "grey", xlim = c(5781,16569), ylim=c(16569,5781), xlab = '', ylab = '')
par(new=TRUE)
plot(FinalAllRep[FinalAllRep$RealisedRepeat == 'mb_del',]$RepStart,FinalAllRep[FinalAllRep$RealisedRepeat == 'mb_del',]$RepEnd, pch = 16, col = "red", xlim = c(5781,16569), ylim=c(16569,5781), xlab = 'Start', ylab = 'End')


sp = ggplot(FinalAllRep) +
  geom_point(aes(RepStart, RepEnd, col=RealisedRepeat)) +
  #scale_fill_manual(values=c("#404080", "#69b3a2")) + 
  scale_y_reverse() +
  theme_minimal() + xlab('Start') + ylab('End') + # scale_y_continuous(sec.axis = dup_axis(), breaks=c(15000, 12000, 9000, 6000)) +
  scale_x_continuous(breaks=c(15000, 12000, 9000, 6000), position = 'top') + 
  # theme(
  #   axis.text.x = element_blank(),
  #   axis.text.y = element_blank()) +
  scale_color_manual(values=c('red', 'grey'))


xplot <- ggplot(FinalAllRep, aes(RepStart, fill = RealisedRepeat)) +
  geom_histogram(alpha=0.4, position = 'dodge') + theme_minimal() +
  xlab('') + ylab('') + scale_x_continuous(breaks=c(6000, 9000, 12000, 15000)) +
  scale_y_continuous(breaks = c(100, 200, 300, 400, 500, 600)) +
  scale_fill_manual(values=c('red', 'grey'), name = "", labels = c("Realized repeats", "Non-realized repeats"))

yplot <- ggplot(FinalAllRep, aes(RepEnd, fill = RealisedRepeat)) +
  geom_histogram(alpha=0.4, position = 'dodge') + coord_flip() + scale_x_reverse(breaks=c(15000, 12000, 9000, 6000)) +
  theme_minimal() + xlab('') + ylab('') + 
  scale_fill_manual(values=c('red', 'grey'))

# Cleaning the plots

# Arranging the plot
ggarrange(xplot, NULL, sp, yplot, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)

wilcox.test(FinalAllRep[FinalAllRep$RealisedRepeat == 'mb_del',]$RepStart, FinalAllRep[FinalAllRep$RealisedRepeat == 'non_del',]$RepStart)
wilcox.test(FinalAllRep[FinalAllRep$RealisedRepeat == 'mb_del',]$RepEnd, FinalAllRep[FinalAllRep$RealisedRepeat == 'non_del',]$RepEnd)

FinalAllRep$RealisedRepeat = as.factor(FinalAllRep$RealisedRepeat)
summary(FinalAllRep$RealisedRepeat)

## higher chances of the repeat to be realised if it is in contact zone versus non-contact zone:
StartA = 6000 # 7000
StartB = 9000 # 10000
EndA = 13000
EndB = 16000
R.In = nrow(Rep[Rep$StartOfRealizedRepeats >= StartA & Rep$StartOfRealizedRepeats <= StartB & Rep$EndOfRealizedRepeats >= EndA & Rep$EndOfRealizedRepeats <= EndB,])
R.Out = nrow(Rep[(Rep$StartOfRealizedRepeats < StartA | Rep$StartOfRealizedRepeats > StartB) | (Rep$EndOfRealizedRepeats < EndA | Rep$EndOfRealizedRepeats > EndB),])
R.In + R.Out # 324 + 294 = 618 
NR.In = nrow(Rep[Rep$StartOfNonRealizedRepeats >= StartA & Rep$StartOfNonRealizedRepeats <= StartB & Rep$EndOfNonRealizedRepeats >= EndA & Rep$EndOfNonRealizedRepeats <= EndB,])
NR.Out = nrow(Rep[(Rep$StartOfNonRealizedRepeats < StartA | Rep$StartOfNonRealizedRepeats > StartB) | (Rep$EndOfNonRealizedRepeats < EndA | Rep$EndOfNonRealizedRepeats > EndB),])
NR.In + NR.Out # 79 + 598 = 618
X = cbind(c(R.In,R.Out),c(NR.In,NR.Out))
fisher.test(X) # odds 7.5, p-value < 2.2e-16
X = data.frame(X)
names(X) =c('R','NR')
row.names(X)=c('IN','OUT')
mosaicplot(X, main = '')


dev.off()
