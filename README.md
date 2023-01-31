# ComparativeGenomics

Structure of folders:

```
Body/

  1Raw/         # raw data for reference, nothing to modify here
  
  2Derived/     # derived data = modified raw data for statistical analyses
  
  3Results/     # one table (one file) - one statistical analysis
  
  4Figures/     # only derived figures


Head/

  1KeyPapers/   # pdfs of key papers that inspired the work
  
  2Scripts/     # scripts and documentation
  
  3Conferences/ # our abstracts, posters, slides
  
  4MyPapers/    # our drafts, paper




============Scripts============

1) Tandem Repeats analysis
LongShortTR.R (IN:GenerationLengthForMammals.xlsx,MitGenomics.txt,TRFinder.txt; OUT:LongShortTR.pdf)

2) CentersOfDeletions
RealandRandomCentersOfDeletions.py
MitoBreak_vs_RandomInMajorArc.R

3) Realized Vs Nonrealized
RealizedVsNonrealizedDeletions.R
RealizedVsNonrealizedDeletions.tex

4) DIID
MitoBreakDeletionsAndInteractionOfDirectAndInvertedOrlovRepeats.R
MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.DiidForEachDD.R
MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.DiidForEachDD.RnaFolding.R
MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.R
MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.SearchForBestDiid.R
MitoBreakDeletionsAndInteractionOfDirectAndInvertedVictorRepeats.ShortInvInvINextToFirstDir.R
MitoBreakDeletionsAndTwoRnaFoldingsSimulatingInvIInvAndDirDir.R

Other:
BranchLengthVsGenerationLength.R
IN:GenerationLengthForMammals.xlsx
OUT:BranchLengthDiff.pdf

Clustering.ipynb
IN:MitoBreakDB_12122019.csv
OUT:cluster_spike_stats.csv, clusters.png

CopelandThirdComponentPca.R 
IN:PCA_vectors_c0toc2_01.txt
OUT:CopelandThirdComponentPca1.R.pdf,CopelandThirdComponentPca2.R.pdf

KnksGenTimelm.R
IN:MitGenomics.txt,KaKsRL_data.txt
OUT: Merged data

lin_model.R
IN:MitGenomics.txt,GenerationLengthForMammals.xlsx,mtalign.aln.treefile.rooted
OUT:lin_model.R.NucleotidesPairwiseCorrs.pdf,spearman

logRegDeletionProbability.R
IN:MitoBreakDB_12122019.csv,compare.square.200.38windows.txt,compare.square.200.38windows_2.txt,score_out_matrix; 
OUT:Bublik

MajorArcAlign.py
IN:1_mitoref_majorarc.fasta
OUT:output

MajorArcHeavyChainRandomizator.R
IN:heavy_major.fasta
OUT:MajorArcHeavyChainRandomizator.StartFromNucContent.txt

MajorArcParseAlign.py
IN:score_out
OUT:output

MammalsGenomesEcology.R
IN:GenerationLengthForMammals.xlsx,MitGenomics.txt,mammalia_genomes.fa
OUT:MammalsGenomesWithEcology.fasta,MammalsGenomesEcology.txt

MitoBreakDeletionsAndOrlovRepeats.R
IN:MitoBreakDB_12122019.csv,Reverse1_100bp.x10.Complement2_100bp/result/,QuadruplexFormingSequences/Homo_sapiens.genome.cut.gff
OUT:MitoBreakDeletionsAndTwoRnaFoldingsSimulatingInvIInvAndDirDir.R.pdf,InvInvFoldings.Complement1_100bp.100bp.Complement2_100bp.txt,DirDIrFoldings.Reverse1_100bp.x10.Complement2_100bp.txt

NaiveSimulation.RepeatsAsAFunctionOfNuclFractions.R
IN:RepsCountPseudoStep.csv,RepsCount.csv
OUT:NaiveSimulation.RepeatsAsAFunctionOfNuclFractions.R.01.pdf

QuadruplexesAlongMtDnaDecription.R
IN:Homo_sapiens.genome.cut.gff
OUT:QuadruplexesAlongMtDnaDecription.R01.pdf

RealAndSimulatedRepeatsVsLongevity.R
IN:ecoreps.csv,mtalign.aln.treefile.rooted
OUT:RealAndSimulatedRepeatsVsLongevity.R.01.pdf

SlipAndJump.R
IN:100x100.csv,Link_matrix_direct_major_activ_left.csv,MitoBreakDB_12122019.csv,Link_matrix1000_major.csv,Link_matrix100hydra_major.csv,Link_matrix_1000_invert_major_activ_left.csv,Link_matrix_invert_major_activ_left.modified.csv
OUT:SlipAndJump.R.01.pdf,heatmap_global_folding_sw100.pdf,violin_rep_folding_infsign_np.pdf,violin_rep_folding_infsign_p.pdf,SlipAndJump.HomologyAndRepeats.txt,SlipAndJump.R.02.pdf,heatmap_microhomology_AIC.pdf,heatmap_microhomology_AIC_wt_deledions.pdf

testPhylogeneticContrasts.R
OUT:spearman)
  
```
-------------------------------------------------------------------------

Center For Mitochondrial Functional Genomics, School of Life Science, 
Immanuel Kant Federal Baltic University, 
Kaliningrad, Russia.