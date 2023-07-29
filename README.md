scNOVA
====================================

scNOVA : Single-Cell Nucleosome Occupancy and genetic Variation Analysis - summarized in a single [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


## Overview of this workflow
PART0. Pre-requirement step - preparation of single-cell genetic information using scTRIP analysis (Sanders et al. 2020) <br>
Mosaicatcher: https://github.com/friendsofstrandseq/mosaicatcher-pipeline
<br/><br/>
PART1. Read preprocessing
<br/><br/>
PART2. Read counting to generate single-cell genebody NO
<br/><br/>
PART3. Infer expressed genes of subclones
<br/><br/>
PART4. DE analysis of subclones (default and alternative mode)
<br/><br/>
PART5. (Optional) Infer Single-cell TF motif accessibility using chromVAR (20210108 updated), by default, Roadmap epigenomics DHS (Enhancers) will be used to define CREs
<br/><br/>
PART6. (Optional) Infer haplotype-resolved genebody NO (20210108 updated)
<br/><br/>
Main output
1. Single-cell level NO table : `result/{SAMPLE}_sort_geneid.txt`
2. Infer expression probability for each clones : `result_CNN/DNN_train80_output_ypred_clone1_annot.txt`,<br> `result_CNN/DNN_train80_output_ypred_clone2_annot.txt`
3. Infer differential expression table : `result/Result_scNOVA_infer_expression_table.txt`, <br> `result/result_PLSDA_{SAMPLE}.txt` (20220121 updated)
4. Heatmap and UMAP visualization of significant hits : `result_plots/Result_scNOVA_plots_{SAMPLE}.pdf`,<br>`result_plots/Result_scNOVA_plots_{SAMPLE}_alternative_PLSDA.pdf`
5. (Optional) Single-cell level TF motif deviation z-score : `result/motif_dev_zscore_chromVAR_DHS_2kb_Enh_{SAMPLE}.txt`
6. (Optional) Haplotype-resolved NO of genebody and CREs : `result_haplo/Deeptool_Genebody_H1H2_sort.txt, result_haplo/Deeptool_DHS_2kb_H1H2_sort.txt`
<br/><br/> 
## System requirements
This workflow is mean to be run in a Unix-based operating system.


## Installation
As of now, unix based tools and python packages are expected to be available via modules. R packages are automatically downloaded and used via conda.
1. unix based tools : SAMtools/1.3.1-foss-2016b, biobambam2/2.0.76-foss-2016b, deeptools/2.5.1-foss-2016b-Python-2.7.12
2. python packages : cuDNN, CUDA, TensorFlow, scikit-learn, matplotlib
3. R packages: handled via conda install. You can choose not to use conda, in which case you will need the following packages: DESeq2, matrixStats, pheatmap, gplots, umap, Rtsne, factoextra, pracma, chromVAR, nabor, motifmatchr 

## Setup
1. **Download this pipeline**
	* git lfs install
	* git clone https://github.com/jeongdo801/scNOVA.git
        * install dependencies (see further below)
2. **Preparation of input files**
	* Add your single-cell bam and index files (input_bam/*.bam)
	* Add key result files from mosaicatcher output in the input_user folder
		* input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE.txt
		* input_user/strandphaser_output.txt
	* Add the subclonality information (input_user/input_subclonality.txt)
	* Add the genes within copy number changed region to mask in the infer differential expression result, if it's not provided, genes will not be masked. (input_user/input_SV_affected_genes.txt) 

3. **Change the project name in the Snakefile**
4. **Setup dependencies (--use-conda is recommended as follows)**
	* `snakemake -j 4 --use-conda --conda-create-envs-only`
5. **Launch the run_pipeline.sh script**

## Dependencies and conda
This tool reduces the burden of installing R dependencies and their correct versions by using conda. After following the installation and setup (steps 1-3), use `snakemake -j 4 --use-conda --conda-create-envs-only` to prepare the conda environment.

Note 1: If you do not have a lot of space in your `HOME` directory (i.e. `~`), make sure conda installs packages into a directory where you have more space. To do so, open `~/.condarc` (which already might have content if you have used conda before) and add the following lines to it: 
```
pkgs_dirs:
  - <path/to/directory_with_more_space>/envs/pkgs/
envs_dirs:
  - <path/to/directory_with_more_space>/envs/
```

Note 2: If you want to install all dependencies manually and not use conda, simply remove `--use-conda` from `run_pipeline.sh`, which will make `snakemake` ignore all `conda` statements.

### Rare conda issue
A note on running hundreds of jobs at the same time (cluster environment):

If running many concurrent jobs, a rare race condition can occur in which two R environments are activated at the same time. During this activation, a file is deleted for a short amount of time, which will throw an error if a separate activation is taking place at the same time. The issue is described [here](https://github.com/conda-forge/r-base-feedstock/issues/67) and can be worked around by applying [this fix](https://github.com/kpalin/r-base-feedstock/commit/9eda35bdc8ea2c2433cbc6b94c2e978b4d7cd8d4), which has not yet been merged into any branch unfortunately.

The error that gets thrown looks like this: `/path/to/pipeline/.snakemake/conda/<hash>/lib/R/bin/R: line 248: /path/to/pipeline/.snakemake/conda/<hash>/lib/R/etc/ldpaths: No such file or directory`

## Configuration for the special use
By default, to make CNN feature, scNOVA performs library size normalization followed by local copy number normalization.
However, for the samples with more dramatic karyotypic changes (e.g. changes in ploidy status), we provide a option to use copy number normalization before normalization by library size. To do so, users can change two lines in the `Snake.config.json`
```
{
    ...
    "count_norm" : "utils/count_norm_ploidy.R",
    ...
    "combine_features" : "utils/combine_features_ploidy.R",
    ...
}
```


## References
For detailed information on scNOVA see

> Jeong and Grimes *et al.*, 2022 (https://www.nature.com/articles/s41587-022-01551-4)


For information on scTRIP see

> Sanders *et al.*, 2020 (doi: https://doi.org/10.1038/s41587-019-0366-x)





