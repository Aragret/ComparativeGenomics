library(ape)
library(gdata)

GenLength <- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")
GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

tree = read.tree('../../Body/2Derived/mtalign.aln.treefile.rooted')

str(tree$edge)

tree$edge[,2]

numberOfSpecies = length(tree$tip.label)

a = as.data.frame(tree$edge)
a = cbind(a, tree$edge.length)
externalBranches = a[a$V2 <= numberOfSpecies,]
externalBranches = cbind(externalBranches, tree$tip.label)

names(externalBranches) = c('Node', 'Tip', 'BranchLength', 'Species')

externalBranchesGl = merge(externalBranches, GenLength, by='Species')

cor.test(externalBranchesGl$GenerationLength_d, externalBranchesGl$BranchLength, method = 'spearman')
# rho -0.2042569
# pvalue 4.364e-08
