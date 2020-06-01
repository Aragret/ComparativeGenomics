rm(list=ls(all=TRUE))

pdf("../../Body/4Figures/QuadruplexesAlongMtDnaDecription.R01.pdf")
  
GQ = read.table("../../Body/1Raw/QuadruplexFormingSequences/Homo_sapiens.genome.cut.gff", sep = '\t', header = FALSE)
GQ=GQ[,-c(1,2)]
names(GQ)=c('start','width','score','strand','type')
plot(GQ$start,GQ$score) # vertical lines mean many GQ with similar start and different scores 
  plot(GQ$start,GQ$score,xlim=c(6000,16560)) # vertical lines mean many GQ with similar start and different scores 


dev.off()  
