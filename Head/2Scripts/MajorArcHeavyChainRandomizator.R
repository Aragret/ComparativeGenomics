  rm(list=ls(all=TRUE))

  MajArc = read.table("../../Body/1Raw/heavy_major.fasta", sep = '\t', header = TRUE)
  for (i in 1:nrow(MajArc))
  { # i = 1
    if (i == 1) {VecMajArc = as.character(MajArc[i,])}
    if (i >  1) {VecMajArc = c(VecMajArc,as.character(MajArc[i,]))}
  }
  
  VecMajArc = unlist(strsplit(VecMajArc,''))
  length(VecMajArc)
  table(VecMajArc)
  #   A    C    G    T 
  # 3011 3226 1208 2556 
  
  RandMajArc = sample(c(rep('A',3011),rep('C',3226),rep('G',1208),rep('T',2556)))
  length(RandMajArc)
  for (i in 1:10000) {RandMajArc = sample(RandMajArc)}
  
  RandMajArc = data.frame(paste(RandMajArc, collapse = ''))
  names(RandMajArc) = c('>MajorArcHeavyChainRandomizator.StartFromNucContent.txt')
  write.table(RandMajArc,"../../Body/2Derived/MajorArcHeavyChainRandomizator.StartFromNucContent.txt", row.names = FALSE, quote = FALSE)
  