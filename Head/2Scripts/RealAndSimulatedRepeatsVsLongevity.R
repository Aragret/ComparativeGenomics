rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

Rep = read.table('../../Body/3Results/ecoreps.csv', sep = ';', header = TRUE)

pdf("../../Body/4Figures/RealAndSimulatedRepeatsVsLongevity.R.01.pdf")

VecOfP = c(); VecOfRho = c(); 
for (i in 8:ncol(Rep))
{ # i = 8
p = as.numeric(cor.test(Rep[,i],Rep$generlen, method = 'spearman')[3])
rho = as.numeric(cor.test(Rep[,i],Rep$generlen, method = 'spearman')[4])
VecOfP = c(VecOfP,p)
VecOfRho = c(VecOfRho,rho)
}
summary(VecOfP)
summary(VecOfRho)

plot(VecOfRho[2:101],-log10(VecOfP[2:101]), ylim=c(min(-log10(VecOfP)),max(-log10(VecOfP))), xlim = c(min(VecOfRho),max(VecOfRho)), pch = 16, col = rgb(0.1,0.1,0.1,0.5), xlab = '', ylab = ''); par(new=TRUE)
plot(VecOfRho[1],-log10(VecOfP[1]), ylim=c(min(-log10(VecOfP)),max(-log10(VecOfP))), xlim = c(min(VecOfRho),max(VecOfRho)), pch = 16, col = rgb(1,0,0,1), xlab = 'rho', ylab = '-log10(p value)');

hist(VecOfRho[2:101], breaks = 15, xlim = c(min(VecOfRho),max(VecOfRho)), col = 'grey')
abline(v = VecOfRho[1], col = 'red', lwd = 3)

### do the same with PICs and probably plot it (should look better):
par(mfrow=c(1,2))
plot(log2(Rep$generlen),log2(Rep$TLOfADRRreal))
plot(log2(Rep$generlen),log2(Rep$TLOfADRR1))

# dev.off()

#########################################################################################
######################### PICs

library(ape) # install.packages('ape') 

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

data = Rep[which(as.character(Rep$spece) %in% tree$tip.label),]
row.names(data) = data$spece

df_vec <- as.character(Rep$spece)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

tree2 <- drop.tip(tree, b)


### rho vs p 

VecOfP = c(); VecOfRho = c(); 
for (i in 8:ncol(Rep))
{ # i = 8
  p = as.numeric(cor.test(pic(log2(data[,i]), tree2), pic(log2(data$generlen), tree2), method = 'spearman')[3])
  rho = as.numeric(cor.test(pic(log2(data[,i]), tree2), pic(log2(data$generlen), tree2), method = 'spearman')[4])
  VecOfP = c(VecOfP,p)
  VecOfRho = c(VecOfRho,rho)
}
summary(VecOfP)
summary(VecOfRho)

par(mfrow=c(1,1))
plot(VecOfRho[2:101],-log10(VecOfP[2:101]), ylim=c(min(-log10(VecOfP)),max(-log10(VecOfP))), xlim = c(min(VecOfRho),max(VecOfRho)), pch = 16, col = rgb(0.1,0.1,0.1,0.5), xlab = '', ylab = ''); par(new=TRUE)
plot(VecOfRho[1],-log10(VecOfP[1]), ylim=c(min(-log10(VecOfP)),max(-log10(VecOfP))), xlim = c(min(VecOfRho),max(VecOfRho)), pch = 16, col = rgb(1,0,0,1), xlab = 'rho', ylab = '-log10(p value)');

par(mfrow=c(1,1))
hist(VecOfRho[2:101], breaks = 15, xlim = c(min(VecOfRho),max(VecOfRho)), col = 'grey')
abline(v = VecOfRho[1], col = 'red', lwd = 3)

##############

maxRho = match(min(VecOfRho), VecOfRho)

TempData = data[, c('generlen', 'TLOfADRRreal', 'TLOfADRR1')]
# cor.test(Rep$generlen, Rep$TLOfADRR1, method='spearman')
TempData[,1] = log2(TempData[,1]); TempData[,2] = log2(TempData[,2]); TempData[,3] = log2(TempData[,3]); 
contrasts <- as.data.frame(apply(TempData, 2, pic, tree2))
names(contrasts) = names(TempData)
cor.test(contrasts$generlen, contrasts$TLOfADRR1, method='spearman')
cor.test(contrasts$generlen, contrasts$TLOfADRRreal, method='spearman')

library(pacman)
p_load(tibble, dplyr, magrittr, purrr)
contrasts <- data %>% 
  select(generlen, TLOfADRRreal, TLOfADRR1) %>% 
  mutate_if(is.numeric, log2) %>% 
  map(pic, tree2)


summary(contrasts$generlen)
summary(pic(log2(data$generlen), tree2))

a = as.data.frame(contrasts)
# newdata <- subset(contrasts, !(contrasts$generlen > quantile(contrasts$generlen, probs=c(.03, .97))[2] | contrasts$generlen < quantile(contrasts$generlen, probs=c(.01, .99))[1]) ) 

newdata <- subset(a, !(a$generlen > quantile(a$generlen, probs=c(.03, .97))[2] | a$generlen < quantile(a$generlen, probs=c(.01, .99))[1]) ) 
newdata <- subset(newdata, !(newdata$TLOfADRRreal > quantile(newdata$TLOfADRRreal, probs=c(.03, .97))[2] | newdata$TLOfADRRreal < quantile(newdata$TLOfADRRreal, probs=c(.01, .99))[1]) ) 
newdata <- subset(newdata, !(newdata$TLOfADRR1 > quantile(newdata$TLOfADRR1, probs=c(.03, .97))[2] | newdata$TLOfADRR1 < quantile(newdata$TLOfADRR1, probs=c(.01, .99))[1]) ) 

par(mfrow=c(2,2))

# plot(contrasts$generlen, contrasts$TLOfADRRreal, col = rgb(0.1,0.1,0.1,0.5)); # cor.test(contrasts$generlen, contrasts$TLOfADRRreal, method = 'spearman', alternative = 'less') # nonsignificant
# plot(contrasts$generlen, contrasts$TLOfADRR1, col = rgb(0.1,0.1,0.1,0.5)); # cor.test(contrasts$generlen, contrasts$TLOfADRR1, method = 'spearman', alternative = 'less')       # marginally

plot(newdata$generlen, newdata$TLOfADRRreal, col = rgb(0.1,0.1,0.1,0.1), pch = 16, cex = 2); # cor.test(newdata$generlen, newdata$TLOfADRRreal, method = 'spearman', alternative = 'less') # nonsignificant
plot(newdata$generlen, newdata$TLOfADRR1, col = rgb(0.1,0.1,0.1,0.1), pch = 16, cex = 2,
     ylim = c(min(newdata$TLOfADRRreal), max(newdata$TLOfADRRreal))); cor.test(newdata$generlen, newdata$TLOfADRR1, method = 'spearman', alternative = 'less') # nonsignificant

dev.off()

