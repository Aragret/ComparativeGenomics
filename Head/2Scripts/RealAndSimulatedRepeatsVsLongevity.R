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

library(ape)

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

data = Rep[which(as.character(Rep$spece) %in% tree$tip.label),]
row.names(data) = data$spece

df_vec <- as.character(Rep$spece)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

tree2 <- drop.tip(tree, b)

TempData = data[, c('generlen', 'TLOfADRRreal', 'TLOfADRR1')]
contrasts <- as.data.frame(apply(TempData, 2, pic, tree2))
names(contrasts) = names(TempData)

par(mfrow=c(2,1))
plot(log(contrasts$generlen), log(contrasts$TLOfADRRreal))
plot(log(contrasts$generlen), log(contrasts$TLOfADRR1))

dev.off()
