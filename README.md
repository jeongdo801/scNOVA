scNOVA
====================================

scNOVA : Single-Cell Nucleosome Occupancy and genetic Variation Analysis - summarized in a single [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


## Overview of this workflow
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

 
## System requirements
This workflow is mean to be run in a Unix-based operating system.


## Installation
1. unix based tools : SAMtools/1.3.1-foss-2016b, biobambam2/2.0.76-foss-2016b, deeptools/2.5.1-foss-2016b-Python-2.7.12
2. python packages : cuDNN, CUDA, TensorFlow, scikit-learn, matplotlib
3. R packages : DESeq2, matrixStats, pheatmap, gplots, umap, Rtsne, factoextra, pracma 


## Setup
1. Download this pipeline
2. Add your single-cell data (input_bam folder)
3. Add the subclonality information (input_user folder)
4. Change the project name in the Snakefile
5. Launch the run_pipeline.sh script

## References

For information on scTRIP see

> Sanderes *et al.*, 2019 (doi: https://doi.org/10.1038/s41587-019-0366-x)





