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


############## how to get sisters species

library(geiger)

max_node_number = max(tree$edge)
min_node_number = length(tree$tip.label) + 1

one_line = c()
for (i in min_node_number:max_node_number){
  descendants = tips(tree, i)
  if (length(descendants) == 2){
    one_line = rbind(one_line, c(i, descendants))
  }
}

sisters = as.data.frame(one_line)
names(sisters) = c('NodeNumber', 'Species1', 'Species2')


###########

sistersGL = merge(sisters, GenLength, by.x='Species1', by.y = "Species")
sistersGLfull = merge(sistersGL, GenLength, by.x='Species2', by.y = "Species")

branchLengthsVec = c(NULL)
for(i in 1:nrow(sistersGLfull)){
  # i = 1
  diff = sistersGLfull[i, 'GenerationLength_d.x'] - sistersGLfull[i, 'GenerationLength_d.y']
  if(diff > 0){
    longsp = as.character(sistersGLfull[i, 'Species1'])
    shortsp = as.character(sistersGLfull[i, 'Species2'])
  }
  if(diff < 0){
    longsp = as.character(sistersGLfull[i, 'Species2'])
    shortsp = as.character(sistersGLfull[i, 'Species1'])
  }
  branchLongsp = externalBranchesGl[externalBranchesGl$Species == longsp, 'BranchLength']
  branchShortsp = externalBranchesGl[externalBranchesGl$Species == shortsp, 'BranchLength']
  a = branchLongsp - branchShortsp
  branchLengthsVec = c(branchLengthsVec, a)
}

pdf('../../Body/4Figures/BranchLengthDiff.pdf')
summary(branchLengthsVec)
hist(branchLengthsVec, breaks = 50)

dev.off()
