rm(list=ls(all=TRUE))

library(ape)

#forTree = "((((((A, B), (C, D)), ((E, F), (G, H))), (((I, J), (K, L)), ((M, N), (O, P)))),(((Q, R), (S, T)),((U, V), (W, X)))),(Y,Z));"

forTree = "((((((Z, Y), (X, W)), ((V, U), (T, S))), (((R, Q), (P, O)), ((N, M), (L, K)))),(((J, I), (H, G)),((F, E), (D, C)))),(B,A));"


# forTree = "((((A, B), (C, D)), ((E, F), (G, H))), (((I, J), (K, L)), ((M, N), (O, P))));"

tree = read.tree(text=forTree)

plot(tree,no.margin=TRUE,edge.width=2)
plot(tree,no.margin=TRUE,edge.width=2); tiplabels()
plot(tree,no.margin=TRUE,edge.width=2); nodelabels()
plot(tree,no.margin=TRUE,edge.width=2); tiplabels(); nodelabels()

species = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z')
species = rev(species) # to follow the tree: A and B are ancestors and we start vectors from 

# species = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P')

# feature1 = c(1, 1.5, 3, 3.5, 7, 7.5, 8.5, 9, 140, 143, 150, 153, 180, 183, 190, 193, 1000, 1100, 1300, 1400, 2000, 2100, 2300, 2400)
# feature2 = c(1, 1.5, 3, 3.5, 7, 7.5, 8.5, 9, 140, 143, 150, 153, 180, 183, 190, 193, 1000, 1100, 1300, 1400, 2000, 2100, 2300, 2400)

GenLengthContinuosIncrease = rev(c(1:26))
GenLengthStepIncrease = rev(c(rep(1,2),rep(2,8),rep(3,8),rep(4,8)))

GcContentMtDnaContinuosIncrease  = rev(seq(0.01,0.26,0.01))
GcContentMtDnaStepIncrease  = rev(c(rep(0.1,2),rep(0.2,8),rep(0.3,8),rep(0.4,8)))

# correlations without PIC normalization: 
cor.test(GenLengthContinuosIncrease, GcContentMtDnaContinuosIncrease, method = 'spearman') # 1
cor.test(GenLengthStepIncrease, GcContentMtDnaContinuosIncrease, method = 'spearman') # 0.96
cor.test(GenLengthStepIncrease, GcContentMtDnaStepIncrease, method = 'spearman') # 1 

# PICs
data = as.data.frame(list(species, GenLengthStepIncrease, GcContentMtDnaStepIncrease))
# data = as.data.frame(list(species, GenLengthContinuosIncrease, GcContentMtDnaStepIncrease))
# data = as.data.frame(list(species, GenLengthContinuosIncrease, GcContentMtDnaContinuosIncrease))

names(data) = c('species', 'feature1', 'feature2')
treeComputeBranchLength = root(compute.brlen(tree), node = 27) # node of the root = number of species + 1

a = pic(data$feature1, treeComputeBranchLength); 
a                           # it starts from the root and goes to shallow nodes
length(a);                  # number of contrasts = number of internal nodes = number of species - 1
summary(a);
a = as.numeric(a)
b = pic(data$feature2, treeComputeBranchLength);
b                           # it starts from the root and goes to shallow nodes
length(b);                  # number of contrasts = number of internal nodes = number of species - 1
summary(b);
b = as.numeric(b)

### if we dont' want to analyze deep contrasts - how we choose nodes, close to external? ALINA!!!!
### if we want to visualize data (traits) using branch lengths - how to do it? ALINA!!!!

### ANALYSIS OF CONTRASTS
## if both cotrasts are zero - we have to delete them (!!), becuase esentialy it means that both features are not changing at all.

# rank correlation gives too good results!!!! Almost the whole vector consists of zeroes!!! How it is possible
cor.test(a, b)
cor.test(a, b, method = 'spearman')
cor.test(a, b, method = 'kendall') 

cor.test(a[3:25], b[3:25])
cor.test(a[3:25], b[3:25], method = 'spearman') 
cor.test(a[3:25], b[3:25], method = 'kendall') 

# lm?
A = lm(a~ 0 + b)
summary(A)

############## ALINA's ORIGINAL CODE 


#feature1 = c(140, 143, 150, 153, 180, 183, 190, 193, 1000, 1100, 1300, 1400, 2000, 2100, 2300, 2400)
#feature2 = c(1, 1.5, 3, 3.5, 7, 7.5, 8.5, 9, 140, 143, 150, 153, 180, 183, 190, 193)
# feature1 = c(143, 140, 150, 153, 180, 183, 190, 193, 1100, 1000, 1300, 1400, 2000, 2100, 2300, 2400)
# feature2 = c(1.5, 1, 3, 3.5, 7, 7.5, 8.5, 9, 143, 140, 150, 153, 180, 183, 190, 193)

cor.test(feature1, feature2) # 0.9633101
cor.test(feature1, feature2, method = 'spearman') # 1

data = as.data.frame(list(species, feature1, feature2))

names(data) = c('species', 'feature1', 'feature2')

treeComputeBranchLength = root(compute.brlen(tree), node = 17)

a = pic(feature1, treeComputeBranchLength)
b = pic(feature2, treeComputeBranchLength)

cor.test(a, b) # 0.8619552
cor.test(a, b, method = 'spearman') # 0.9990719

library(pacman)
p_load(tibble, dplyr, magrittr, purrr)
contrasts <- data %>% 
  select(feature1, feature2) %>% 
  #   mutate_if(is.numeric, log2) %>% 
  map(pic, treeComputeBranchLength)

summary(a)
summary(contrasts$feature1)

summary(b)
summary(contrasts$feature2)

# all contrasts are negative when a vector contains ascending numbers. After rearragment 
# results in PICs are changed (simple corr doesn't changed)