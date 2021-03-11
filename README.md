scNOVA
====================================

scNOVA : Single-Cell Nucleosome Occupancy and genetic Variation Analysis - summarized in a single [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


## Overview of this workflow
PART0. Pre-requirement step - preparation of single-cell genetic information:
Mosaicatcher: https://github.com/friendsofstrandseq/mosaicatcher-pipeline
<br/><br/>
PART1. Read preprocessing
<br/><br/>
PART2. Read counting to generate single-cell genebody NO
<br/><br/>
PART3. Infer expressed genes of subclones
<br/><br/>
PART4. DE analysis of subclones
<br/><br/>
PART5. (Optional) Infer Single-cell TF motif accessibility using chromVAR (20210108 updated), by default, Roadmap epigenomics DHS (Enhancers) will be used to define CREs
<br/><br/>
PART6. (Optional) Infer haplotype-resolved genebody NO (20210108 updated)
<br/><br/>
Main output (GM20509 example data)
1. Single-cell level NO table : result/GM20509_sort_geneid.txt (20210108 updated)
2. Infer expression probability for each clones : result_CNN/DNN_train80_output_ypred_clone1_annot.txt, result_CNN/DNN_train80_output_ypred_clone2_annot.txt
3. Infer differential expression table : result/Result_scNOVA_infer_expression_table.txt (20210108 updated)
4. Heatmap and UMAP visualization of significant hits : result_plots/Result_scNOVA_plots_GM20509.pdf
5. (Optional) Single-cell level TF motif deviation z-score : result/motif_dev_zscore_chromVAR_DHS_2kb_Enh_GM20509.txt (20210108 updated)
6. (Optional) Haplotype-resolved NO of genebody and CREs : result_haplo/Deeptool_Genebody_H1H2_sort.txt, result_haplo/Deeptool_DHS_2kb_H1H2_sort.txt (20210108 updated)
<br/><br/> 
## System requirements
This workflow is mean to be run in a Unix-based operating system.


## Installation
1. unix based tools : SAMtools/1.3.1-foss-2016b, biobambam2/2.0.76-foss-2016b, deeptools/2.5.1-foss-2016b-Python-2.7.12
2. python packages : cuDNN, CUDA, TensorFlow, scikit-learn, matplotlib
3. R packages : DESeq2, matrixStats, pheatmap, gplots, umap, Rtsne, factoextra, pracma, chromVAR, nabor, motifmatchr 


## Setup
1. **Download this pipeline**
	* git lfs install
	* git clone https://github.com/jeongdo801/scNOVA.git

2. **Preparation of input files**
	* Add your single-cell bam and index files (input_bam/*.bam)
	* Add key result files from mosaicatcher output in the input_user folder
		* input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE.txt
		* input_user/strandphaser_output.txt
	* Add the subclonality information (input_user/input_subclonality.txt)
	* Add the genes within copy number changed region to mask in the infer differential expression result, if it's not provided, genes will not be masked. (input_user/input_SV_affected_genes.txt) 

3. **Change the project name in the Snakefile**
4. **Launch the run_pipeline.sh script**

## References

For information on scTRIP see

> Sanders *et al.*, 2019 (doi: https://doi.org/10.1038/s41587-019-0366-x)





