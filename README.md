# Epigenome_analysis_infer_transcriptome

Goal of the pipeline
<br/><br/>
Analysis pipeline (workflow management using Snakemake)
  - Preprocessing to remove duplicate and supplementary reads (samtools, biobambam)
  - Input for the inference: pre-processed bam files (bam)
<br/><br/>
PART0. Extract single-cell and subclonal copy-number variations
<br/><br/>
PART1. Infer transcriptome with pseudo-bulk population
  - Feature1 : samtools merge for pseudo-bulk NO 
    - Create folders/subclone and copy single-cell libraries to the clonality, merge bam files in each folders, run Strand_seq_deeptool_Genes_for_CNN.pl for merged libraries
    - For normalization total mapped read needs to be calculated (Strand_seq_deeptool_chr_length.pl)
    - Make_ML_Features_BCLL01 (processing and normalization to make Feature1)
    - Output: Features_reshape_BCLL01_P1P2_C3_orientation_norm.txt
  - Feature2 : single-cell variance 
    - For each folders of subclone, run Strand_seq_deeptool_Genes_for_CNN.pl for single-cell libraries
    - Deeptool_matrix_CNN.R to calculate coefficient of variation of each bins
    - Make_ML_Features_sc_var_BCLL01 (processing and normaliization to make Feature2)
    - Output: Features_reshape_BCLL01_C3_Resid_orientation.txt
  - Combine five layers of feature sets 
    - Make_ML_Features_BCLL01_combine.R
    - For NO, copy-number normalization will be performed
<br/><br/>
PART2. Infer transcriptome at the single-cell level
  - Strand_seq_deeptool_Genes_for_CNN.pl
  - R_ML_Features_sc.R (This is for normalization)
  - Combine four layers of feature sets
<br/><br/>  
PART3. Infer transcriptome using CNN
  - Deeplearning_Nucleosome_with_mono_var_GC_CpG_RT_leave_Chr1_out_BCLL01_CNnorm_sc_wovar_ypred.py
