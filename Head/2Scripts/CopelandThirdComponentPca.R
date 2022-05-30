rm(list=ls(all=TRUE))
library("ggplot2")   
library("reshape") # install.packages("reshape")                                       # Install reshape package
Pca = read.table("../../Body/1Raw/CopelandDeletionsPca/PCA_vectors_c0toc2_01.txt")

### 1: take only the third component (which is called component 2 in the file because they start from zero) = see "../../Body/1Raw/CopelandDeletionsPca/ReadMe.txt"
# and derive bin coordinates

Pca = Pca[12807:nrow(Pca),]
for (i in 1:nrow(Pca))
{ # i = 1
FirstBin = unlist(strsplit(Pca$V1[i],','))[1]
Pca$FirstBinStart[i] = as.numeric(unlist(strsplit(FirstBin,':'))[1])
Pca$FirstBinEnd[i] = as.numeric(unlist(strsplit(FirstBin,':'))[2])

SecondBin = unlist(strsplit(Pca$V1[i],','))[2]
Pca$SecondBinStart[i] = as.numeric(unlist(strsplit(SecondBin,':'))[1])
Pca$SecondBinEnd[i] = as.numeric(unlist(strsplit(SecondBin,':'))[2])
}

summary(Pca$SecondBinStart - Pca$SecondBinEnd) # 207 is a step!
summary(Pca$FirstBinStart - Pca$FirstBinEnd) # 207 is a step! 

Pca$FirstBinCenter = Pca$FirstBinStart + (Pca$FirstBinEnd - Pca$FirstBinStart)/2
Pca$SecondBinCenter = Pca$SecondBinStart + (Pca$SecondBinEnd - Pca$SecondBinStart)/2

### 2: plot PC3 and extract major arc. Major arc is from Ol (5721) till the end of mtDNA (16569) and a bit more (till Oh: 110)

par(mfrow=c(2,1))
plot(Pca$FirstBinCenter,Pca$V2, pch = 20, cex = 0.5)
plot(Pca$SecondBinCenter,Pca$V2, pch = 20, cex = 0.5)

Pca = Pca[Pca$FirstBinCenter > 5721 & Pca$SecondBinCenter > 5721,]
plot(Pca$FirstBinCenter,Pca$SecondBinCenter)

par(mfrow=c(2,1))
plot(Pca$FirstBinCenter,Pca$V2, pch = 20, cex = 0.5)
plot(Pca$SecondBinCenter,Pca$V2, pch = 20, cex = 0.5)


### plot 1 (rounded to 1kb cells)
str(Pca)
Pca$FirstBinCenterRound  = round(Pca$FirstBinCenter,-3) # till thousands
Pca$SecondBinCenterRound = round(Pca$SecondBinCenter,-3) # till thousands
Agg = aggregate(as.numeric(Pca$V2), by = list(Pca$FirstBinCenterRound,Pca$SecondBinCenterRound), FUN = mean)
names(Agg)=c('Start','End','Value')
head(Agg)
pdf("../../Body/4Figures/CopelandThirdComponentPca1.R.pdf")
ggp <- ggplot(Agg, aes(Start, End)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = Value))
ggp 
dev.off()

table(Agg$Start)
table(Agg$End)


### plot 2 (original cells) PAPER

Agg = aggregate(as.numeric(Pca$V2), by = list(Pca$FirstBinCenter,Pca$SecondBinCenter), FUN = mean)
names(Agg)=c('Start','End','Value')
head(Agg)
pdf("../../Body/4Figures/CopelandThirdComponentPca2.R.pdf")
ggp <- ggplot(Agg, aes(Start, End)) +                           # Create heatmap with ggplot2
  geom_tile(aes(fill = Value))
ggp 
dev.off()

ContactZone = Agg[Agg$Start >= 6000 & Agg$Start <= 9000 & Agg$End >= 13000 & Agg$Start <= 16000,]$Value
length(ContactZone)
wilcox.test(ContactZone,Agg$Value) # 4.48e-13
t.test(ContactZone,Agg$Value)      # 0.001571

table(Agg$Start)
table(Agg$End)
