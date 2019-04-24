library(ggplot2)

Dloops = read.table('../../Body/2Derived/dloops_control_regions.txt', header=TRUE, sep = '\t',
                    row.names = NULL)
names(Dloops) = c('Species', 'Feature_name', 'Feature_location', 'Strand', 'Sequence', 'Notes')
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
CHOR = CHOR[CHOR$TAXON != "AncientFish", ]
####### 

Dloops$Length = as.numeric(lapply(as.character(Dloops$Sequence), nchar))

summary(Dloops$Length)

pdf('../../Body/4Figures/DloopsLength.pdf')

hist(Dloops$Length, breaks = 50)

#######

data = merge(CHOR, Dloops[, c('Species', 'Length')])

for(taxon in unique(data$TAXON)){
  temp_data = data[data$TAXON == taxon,]
  a = cor.test(temp_data$GenomeLength, temp_data$Length)
  plot(temp_data$GenomeLength, temp_data$Length, main = c(taxon, a$estimate))
  print(c(taxon, cor.test(temp_data$GenomeLength, temp_data$Length)))
}

# cor.test(data$GenomeLength, data$Length) # r = 0.6997531

ggplot(data, aes(Length, GenomeLength, col=TAXON)) +
  geom_point() + xlab('D-loops length')


########## Dloop coverage

# data$DloopCoverage = data$Length / data$GenomeLength

# ggplot(data, aes(DloopCoverage, GenomeLength, col=TAXON)) +
#   geom_point() + xlab('D-loops coverage')

########### number of dloops
summary(Dloops$Feature_name)

length(unique(Dloops$Species))

not_unique = data.frame()

for(sp in unique(Dloops$Species)){
  if(nrow(Dloops[Dloops$Species == sp,]) > 1){
    not_unique = rbind(not_unique, Dloops[Dloops$Species == sp,])
  }
}

summary(not_unique$Length)

not_unique = merge(not_unique, CHOR[, c('Species', 'TAXON')])

summary(not_unique$TAXON)

write.table(not_unique, '../../Body/3Results/MultipleDloops.txt', sep='\t')

not_unique$MultipleDloop = 1

multDl = merge(Dloops, not_unique, all.x = TRUE)

multDl[is.na(multDl$MultipleDloop),]$MultipleDloop = 0
multDl = merge(multDl, data[, c('Species', 'GenomeLength')], by='Species', all.x = TRUE)

legend_title <- "Presence of multiple d-loops"
ggplot(multDl, aes(Length, GenomeLength, col=as.factor(MultipleDloop))) +
  geom_point() + xlab('Dloops Length') + 
  labs(colour=legend_title)

####################### 

AGG = aggregate(data$Length, by=list(data$Species), sum)

names(AGG) = c('Species', 'DloopsLength')

oneSpOneCr = merge(AGG, CHOR[, c('Species', 'TAXON', 'GenomeLength', 'ECO.Maximum.longevity..yrs.')], by='Species', all.y = FALSE)

ggplot(oneSpOneCr, aes(DloopsLength, GenomeLength, col=TAXON)) +
  geom_point() + xlab('Sum of multiple dloops')

dev.off()

##########################################################################################
##########################################################################################

a = merge(data, not_unique, all.x = TRUE)
a = a[, c('Species', 'TAXON', 'Length', 'GenomeLength', 'MultipleDloop')]

pdf('../../Body/4Figures/DloopLength.R.Figure1.pdf')

for(taxon in unique(data$TAXON)){
  # taxon = "Actinopterygii"
  temp_data = a[a$TAXON == taxon,]
  # a = cor.test(temp_data$GenomeLength, temp_data$Length)
  
  
  print(ggplot(temp_data, aes(Length, GenomeLength)) +
    geom_point(size = 3, alpha = 0.1) +
    ggtitle(taxon) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlim(min(data$Length), max(data$Length)) + ylim(min(data$GenomeLength), max(data$GenomeLength)) +
    geom_point(aes(Length, GenomeLength, size = MultipleDloop), alpha = 0.3))

}

dev.off()


########################################################################################
#### tables

for(taxon in unique(oneSpOneCr$TAXON)){
  TempData = oneSpOneCr[oneSpOneCr$TAXON == taxon, ]
  print(c(taxon, nrow(TempData)))
  print(summary(lm(TempData$GenomeLength ~ TempData$DloopsLength)))
}

##############################################################
############# add PICs

library(ape) # install.packages('ape') 

tree <- read.tree("../../Body/1Raw/mtalign.aln.treefile.rooted")

data = oneSpOneCr[which(as.character(oneSpOneCr$Species) %in% tree$tip.label),]
row.names(data) = data$Species

df_vec <- as.character(data$Species)
tree_vec <- tree$tip.label

a <- setdiff(df_vec, tree_vec)
b <- setdiff(tree_vec, df_vec)

tree2 <- drop.tip(tree, b)

#############################################################

for(taxon in unique(data$TAXON)){
  # taxon = 'Amphibia'
  TempData = data[data$TAXON == taxon, ]
  
  df_vec <- as.character(TempData$Species)
  tree_vec <- tree$tip.label
  
  a <- setdiff(df_vec, tree_vec)
  b <- setdiff(tree_vec, df_vec)
  
  tree2 <- drop.tip(tree, b)
  
  
  print(c(taxon, nrow(TempData)))
  print(summary(lm(pic(log2(TempData$GenomeLength), tree2) ~ pic(log2(TempData$DloopsLength), tree2), na.action=na.exclude)))
}

