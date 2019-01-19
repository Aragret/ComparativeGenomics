Dloops = read.table('../../Body/2Derived/dloops_control_regions.txt', header=TRUE, sep = '\t',
                    row.names = NULL)
names(Dloops) = c('Species', 'Feature_name', 'Feature_location', 'Strand', 'Sequence', 'Notes')
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')

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

dev.off()



