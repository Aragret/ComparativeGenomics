rm(list=ls(all=TRUE))

library(ape)

forTree = "((((A, B), (C, D)), ((E, F), (G, H))), (((I, J), (K, L)), ((M, N), (O, P))));"

tree = read.tree(text=forTree)

plot(tree,no.margin=TRUE,edge.width=2)
tiplabels()
nodelabels()

species = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P')
# feature1 = c(1, 1.5, 3, 3.5, 7, 7.5, 8.5, 9, 140, 143, 150, 153, 180, 183, 190, 193, 1000, 1100, 1300, 1400, 2000, 2100, 2300, 2400)
# feature2 = c(1, 1.5, 3, 3.5, 7, 7.5, 8.5, 9, 140, 143, 150, 153, 180, 183, 190, 193, 1000, 1100, 1300, 1400, 2000, 2100, 2300, 2400)

feature1 = c(140, 143, 150, 153, 180, 183, 190, 193, 1000, 1100, 1300, 1400, 2000, 2100, 2300, 2400)
feature2 = c(1, 1.5, 3, 3.5, 7, 7.5, 8.5, 9, 140, 143, 150, 153, 180, 183, 190, 193)

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
# results in PICs are chanched (simple corr doesn't changed)