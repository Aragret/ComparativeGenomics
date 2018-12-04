
library(seqinr)
library(gdata)

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
CodonTable = unzip('../../Body/2Derived/AllGenesCodons.zip')
CodonTable = read.table(CodonTable, header=TRUE, sep='\t')
GenLength = read.xls('../../Body/1Raw/GenerationLengthForMammals.xlsx')

GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

### all 13 proteins

CodonTable = CodonTable[CodonTable$Quality == 1,]
for(i in unique(CodonTable$Species)){
  if(nrow(CodonTable[CodonTable$Species == i,]) != 13){
    CodonTable = CodonTable[CodonTable$Species != i,]
  }
}

genes = as.character(unique(CodonTable$Gene))


### ecology

CHOR = CHOR[, c('Species', 'TAXON', 'ECO.Female.maturity..days.', 'ECO.Maximum.longevity..yrs.',
                'ECO.Body.mass..g.')]
CHOR = merge(CHOR, GenLength, by='Species', all.x=TRUE)
data = merge(CodonTable, CHOR, by='Species')

Final <- data[rowSums(is.na(data[, c(10,11,12,13)])) != 4,]

length(unique(Final$Species)) # 1300

Actinopterygii = Final[Final$TAXON == 'Actinopterygii',] # 217 sp.
Reptilia = Final[Final$TAXON == 'Reptilia',] # 93
Amphibia = Final[Final$TAXON == 'Amphibia',] # 32
Mammalia = Final[Final$TAXON == 'Mammalia',] # 731
Aves = Final[Final$TAXON == 'Aves',] # 180

#### write to fasta
write.fasta()
