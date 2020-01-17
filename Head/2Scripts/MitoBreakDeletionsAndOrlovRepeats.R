rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/MitoBreakDeletionsAndOrlovRepeats.R.01.pdf")
  
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
  
## 2: read Orlovs's direct perfect repeats
  
Rep = read.table("../../Body/1Raw/Homo_sapiens.input.out4out.SecondPart", header = TRUE, sep = '\t') # 767
Rep=Rep[Rep$RepName == 'Direct_repeat',] # 330
Rep$RepStart = as.numeric(as.character(Rep$RepStart)); Rep$RepEnd = as.numeric(as.character(Rep$RepEnd)); 
Rep = Rep[Rep$RepStart > 5781 & Rep$RepStart < 16569 & Rep$RepEnd > 5781 & Rep$RepEnd < 16569,] # 171 

## 3: plot repeats and breakpoints
par(mfrow=c(1,1))
plot(breaks$X5..breakpoint,breaks$X3..breakpoint,xlim = c(5781,16569), ylim=c(16569,5781), col = rgb(0.5,0.5,0.5,0.2), pch = 16, xlab = '5\'breakpoint',  ylab = '3\'breakpoint')
par(new= TRUE)
plot(Rep$RepEnd,Rep$RepStart, xlim = c(5781,16569), ylim=c(16569,5781), col = rgb(1,0.1,0.1,0.5), pch = 16, xlab = '5\'breakpoint',  ylab = '3\'breakpoint')
legend(12000, 8000, c('deletions','repeats'), col=c(rgb(0.5,0.5,0.5,0.5),rgb(1,0.1,0.1,0.5)), pch = 16)
# legend( text=)

dev.off()  
