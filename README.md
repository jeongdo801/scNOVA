# Epigenome_analysis_infer_transcriptome

Goal of the pipeline
<br/><br/>
Analysis pipeline (workflow management using Snakemake)
  - Preprocessing to remove duplicate and supplementary reads (samtools, biobambam)
  - Input for the inference: pre-processed bam files (bam)
<br/><br/>
PART1. Infer transcriptome with pseudo-bulk population
  - Feature1 : samtools merge for pseudo-bulk NO, Strand_seq_deeptool_Genes_for_CNN.pl
  - Feature2 : single-cell variance : Deeptool_matrix_CNN_C0.R
<br/><br/>
PART2. Infer transcriptome at the single-cell level
  - Strand_seq_deeptool_Genes_for_CNN.pl
  - R_ML_Features_sc.R
