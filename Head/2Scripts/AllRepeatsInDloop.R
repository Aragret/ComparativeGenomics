rm(list=ls(all=TRUE))  # remove everything from R memory (old variables, datasets...) 

CHOR <- read.table('../../Body/2Derived/MitGenomics.txt', header=TRUE, sep='\t')
CHOR <- CHOR[CHOR$TAXON != 'AncientFish',]

dloops <- read.table('../../Body/2Derived/dloops_control_regions.txt', header=TRUE, sep='\t',
                     row.names = NULL)
names(dloops) <- c('Species', 'Feature_name', 'Feature_location', 'Strand', 'Sequence', 'Notes')
dloops$Length <- as.numeric(lapply(as.character(dloops$Sequence), nchar))

trInDloop <- read.table('../../Body/2Derived/TRinDloops.txt', header=TRUE, sep='\t')
orlovRepeatsInDloop <- read.table('../../Body/2Derived/OrlovRepeatsInDLoop.txt', header=TRUE, sep='\t')

