
library(seqinr)
library(gdata)

CHOR = read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
CodonTable = unzip('../../Body/2Derived/AllGenesCodonUsageNoOverlap.txt.zip')
CodonTable = read.table(CodonTable, header=TRUE, sep='\t')
GenLength = read.xls('../../Body/1Raw/GenerationLengthForMammals.xlsx')

GenLength$Species = gsub(' ','_',GenLength$Scientific_name)
GenLength = GenLength[,c(14,16)]

### all 13 proteins

for(i in unique(CodonTable$Species)){
  if(nrow(CodonTable[CodonTable$Species == i,]) != 13){
    CodonTable = CodonTable[CodonTable$Species != i,]
  }
}

genes = as.character(unique(CodonTable$Gene))


### ecology

CHOR = CHOR[, c('Species', 'ECO.Female.maturity..days.', 'ECO.Maximum.longevity..yrs.',
                'ECO.Body.mass..g.')]
CHOR = merge(CHOR, GenLength, by='Species', all.x=TRUE)
data = merge(CodonTable, CHOR, by='Species')

Final <- data[rowSums(is.na(data[, c(79:82)])) != 4,]

length(unique(Final$Species)) # 1300

#### write to fasta
dir.create('../../Body/2Derived/ForAlign')

for(i in unique(Final$Gene)){
  for(j in unique(Final$Class)){
    GeneTable = Final[Final$Gene == i,]
    GeneTaxonTable = GeneTable[GeneTable$Class == j,]
    write.fasta(as.list(GeneTaxonTable$AminoNoOverlap), GeneTaxonTable$Species, sprintf('../../Body/2Derived/ForAlign/%s_%s.fna', i, j),
                open = 'w')
  }
}

