rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")

library(gdata)
library(Biostrings)
library(seqinr)

GenLength <- read.xls("../../Body/1Raw/GenerationLengthForMammals.xlsx")
CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header = TRUE, sep='\t')

GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

dupl = GenLength[GenLength$Species == "Neophocaena_phocaenoides",]
gl = c(mean(dupl$GenerationLength_d), 'Neophocaena_phocaenoides')
GenLength = GenLength[GenLength$Species != 'Neophocaena_phocaenoides',]
GenLength = rbind(GenLength, gl)
write.table(GenLength, '../../Body/2Derived/GenerationLength.txt', sep='\t', quote = F, row.names = F)

fastaFile <- readDNAStringSet("../../Body/2Derived/mammalia_genomes.fa")
Species = names(fastaFile)
Sequence = paste(fastaFile)
Mamm <- data.frame(Species, Sequence)

data = merge(Mamm, GenLength, by='Species')
data = merge(data, CHOR[, c('Species', 'REP.DirRepLength', 'GenomeLength', 'A', 'T', 'G', 'C', by='Species')])
data = data[, -10]

write.fasta(as.list(data$Sequence), data$Species, '../../Body/2Derived/MammalsGenomesWithEcology.fasta')
