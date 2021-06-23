import os
import os.path

###############################################################################
 # SETTINGS                                                                    #
 #                                                                             #
 # Set the sample name                                                         #
 # Update 20201028                                                                            #
SAMPLE_NAME = "GM20509"
CLONE_NAME = ["clone1", "clone2"]
 #                                                                             #
 #                                                                             #
 # By default, the pipeline will expect file to be in a subfolder called       #
 # 'bam' and to be names *.bam and *.bai                                       #
 #                                                                             #
BAMFILE, = glob_wildcards("input_bam/{cell}.bam")
BAM_SC, = glob_wildcards("input_bam/{single_cells}.sort.mdup.bam")
#                                                                             #
 #                                                                             #
abbreviate_names = False
 #                                                                             #
 ###############################################################################

configfile: "Snake.config.json"

rule all:
    input:
        expand("input_user/sv_calls_all_print.txt"),
        expand("result/Features_reshape_{c}_orientation_CN_correct0.txt", c = CLONE_NAME),
        expand("input_user/sv_calls_all_print_CREs.txt"),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        expand("result/{s}.tab", s = SAMPLE_NAME),
        expand("result/{s}.npz", s = SAMPLE_NAME),
        expand("result/{s}_sort.txt", s = SAMPLE_NAME),
        expand("result/{s}_sort_geneid.txt", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}.tab", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}.npz", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sc.tab", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sc.npz", s = SAMPLE_NAME),
        expand("result/Deeptool_chr_length_{s}.tab", s = SAMPLE_NAME),
        expand("result/Deeptool_chr_length_{s}.npz", s = SAMPLE_NAME),
        expand("result/Deeptool_chr_length_{s}_sc.tab", s = SAMPLE_NAME),
        expand("result/Deeptool_chr_length_{s}_sc.npz", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sort.txt", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sort_lab.txt", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sort_lab_final.txt", s = SAMPLE_NAME),
        expand("result/Features_reshape_{s}_{c}_orientation_norm_qc.pdf", s = SAMPLE_NAME, c = CLONE_NAME),
        expand("result/Features_reshape_{c}_orientation_norm.txt", c = CLONE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sc_sort.txt", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sc_sort_lab.txt", s = SAMPLE_NAME),
        expand("result/Deeptool_Genes_for_CNN_{s}_sc_sort_lab_final.txt", s = SAMPLE_NAME),
        expand("result/Features_reshape_{s}_{c}_Resid_orientation_qc.pdf", s = SAMPLE_NAME, c = CLONE_NAME),
        expand("result/Features_reshape_{c}_Resid_orientation.txt", c = CLONE_NAME),
        expand("result/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_{c}.txt", c = CLONE_NAME),
        expand("result/Expression_all_{c}.txt", c = CLONE_NAME),
        expand("result/Features_reshape_all_TSS_matrix_woM_all_RT_{c}.txt", c = CLONE_NAME),
        expand("result_CNN/DNN_train80_output_ypred_{c}.csv", c = CLONE_NAME),
        expand("result_CNN/DNN_train40_output_ypred_{c}.csv", c = CLONE_NAME),
        expand("result_CNN/DNN_train20_output_ypred_{c}.csv", c = CLONE_NAME),
        expand("result_CNN/DNN_train5_output_ypred_{c}.csv", c = CLONE_NAME),
        expand("result_CNN/DNN_train80_output_ypred_{c}_annot.txt", c = CLONE_NAME),
        expand("result_CNN/DNN_train40_output_ypred_{c}_annot.txt", c = CLONE_NAME),
        expand("result_CNN/DNN_train20_output_ypred_{c}_annot.txt", c = CLONE_NAME),
        expand("result_CNN/DNN_train5_output_ypred_{c}_annot.txt", c = CLONE_NAME),
        expand("result_plots/Result_scNOVA_plots_{s}.pdf", s = SAMPLE_NAME),
        expand("result/{s}_CREs_2kb.tab", s = SAMPLE_NAME),
        expand("result/{s}_CREs_2kb.npz", s = SAMPLE_NAME),
        expand("result/{s}_CREs_2kb_sort.txt", s = SAMPLE_NAME),
        expand("result/{s}_CREs_2kb_sort_num.txt", s = SAMPLE_NAME),
        expand("result/{s}_CREs_2kb_sort_num_sort_for_chromVAR.txt", s = SAMPLE_NAME),
#        expand("result/SC_CN_regions_all_hg38_v2_resize_2kb_sc.txt"),
        expand("result/motif_dev_zscore_chromVAR_DHS_2kb_Enh_{s}.txt", s = SAMPLE_NAME),
        expand("result/motif_dev_chromVAR_DHS_2kb_Enh_{s}.Rda", s = SAMPLE_NAME),
        expand("result/motif_variability_chromVAR_DHS_2kb_Enh_{s}.txt", s = SAMPLE_NAME),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W1.bam", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W2.bam", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C1.bam", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C2.bam", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam.bai", cell=BAMFILE),
        expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam.bai", cell=BAMFILE),
        expand("nucleosome_sampleA/result.H1.bam"),
        expand("nucleosome_sampleB/result.H2.bam"),
        expand("input_user/strandphaser_output_copy.txt"),
        expand("result_haplo/Deeptool_DHS_2kb_H1H2.tab"),
        expand("result_haplo/Deeptool_DHS_2kb_H1H2.npz"),
        expand("result_haplo/Deeptool_DHS_2kb_H1H2_sort.txt"),
        expand("result_haplo/Deeptool_Genebody_H1H2.tab"),
        expand("result_haplo/Deeptool_Genebody_H1H2.npz"),
        expand("result_haplo/Deeptool_Genebody_H1H2_sort.txt"),
#        expand("result_features_sc/CNN_preprocessing_{s}_sc.pdf", s = SAMPLE_NAME),
#        expand("result_features_sc_CN/SC_CN_DHS_for_CNN_{s}.txt", s = SAMPLE_NAME),
#        expand("result_features_sc/Features_reshape_all_TSS_matrix_woM_all_sc.txt"),
 #       expand("result_features_sc/{s}_wovar_exp.txt", s = SAMPLE_NAME),

 #
 # PART 0
 # Preparation of input
 #

rule generate_CN_for_CNN:
    input:
        subclone = "input_user/input_subclonality.txt",
        sv_calls_all = "input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE.txt",
        Deeptool_result_final = "utils/Deeptool_Genes_for_CNN_merge_sort_lab_final.txt",
        CNN_features_annot = "utils/bin_Genes_for_CNN_reshape_annot.txt",
    output:
        sv_calls_all_print = "input_user/sv_calls_all_print.txt",
        CN_result_data1 = expand("result/Features_reshape_{clone}_orientation_CN_correct0.txt", clone = ["clone1", "clone2"]),
    params:
        generate_CN_for_CNN = config["generate_CN_for_CNN"]
    log:
        "log/generate_CN_for_CNN.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.generate_CN_for_CNN} {input.subclone} {input.sv_calls_all} {input.Deeptool_result_final} {input.CNN_features_annot} {output.sv_calls_all_print} 
        """



rule generate_CN_for_chromVAR:
    input:
        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",
        TES_matrix = "utils/Strand_seq_matrix_TES_for_SVM.txt",
        Genebody_matrix  = "utils/Strand_seq_matrix_Genebody_for_SVM.txt",
        DHS_matrix_resize = "utils/regions_all_hg38_v2_resize_2kb_sort_num_sort_for_chromVAR.bed",
        subclone = "input_user/input_subclonality.txt",
        sv_calls_all = "input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE.txt",
    output:
        sv_calls_all_print = "input_user/sv_calls_all_print_CREs.txt",
    params:
        generate_CN_for_chromVAR = config["generate_CN_for_chromVAR"]
    log:
        "log/generate_CN_for_chromVAR.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.generate_CN_for_chromVAR} {input.TSS_matrix} {input.TES_matrix} {input.Genebody_matrix} {input.DHS_matrix_resize} {input.subclone} {input.sv_calls_all} {output.sv_calls_all_print} 
        """


 #
 # PART I
 # Read preprocessing
 #
	
rule remove_low_quality_reads:
    input:
        bam = "input_bam/{cell}.bam"
    output:
        bam_pre = "bam/{cell}.sc_pre_mono.bam",
        bam_header = "bam/{cell}.header_test.sam"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
		samtools view -H {input} > {output.bam_header} 
		samtools view -F 2304 {input.bam} | awk -f utils/awk_1st.awk | cat {output.bam_header} - | samtools view -Sb - > {output.bam_pre}	
        """

rule sort_bam:
    input:
        "bam/{cell}.sc_pre_mono.bam"
    output:
        "bam/{cell}.sc_pre_mono_sort_for_mark.bam"
    threads:
        2
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools sort -@ {threads} -O BAM -o {output} {input}
        """

rule index_num1:
    input:
        "bam/{cell}.sc_pre_mono_sort_for_mark.bam"
    output:
        "bam/{cell}.sc_pre_mono_sort_for_mark.bam.bai"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools index {input}
        """	
	
rule remove_dup:
    input:
        bam="bam/{cell}.sc_pre_mono_sort_for_mark.bam"
    output:
        bam_uniq="bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam",
        bam_metrix="bam/{cell}.sc_pre_mono.metrix_dup.txt"
    shell:
        """
        module load biobambam2/2.0.76-foss-2016b
        bammarkduplicates markthreads=2 I={input.bam} O={output.bam_uniq} M={output.bam_metrix} index=1 rmdup=1
        """

rule index_num2:
    input:
        "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam"
    output:
        "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai"
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools index {input}
        """

 #
 # PART II
 # Read counting
 #

rule count_reads:
    input:
        bam = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        bai = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai", cell=BAMFILE)
    output:
        tab = "result/" + SAMPLE_NAME + ".tab",
        npz = "result/" + SAMPLE_NAME + ".npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/bin_Genebody_all.bed --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """


rule count_sort_by_coordinate:
    input:
        "result/" + SAMPLE_NAME + ".tab"
    output:
        "result/" + SAMPLE_NAME + "_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_sort_annotate_geneid:
    input:
        count_table = "result/" + SAMPLE_NAME + "_sort.txt",
        GB_matrix = "utils/Strand_seq_matrix_Genebody_for_SCDE.txt"
    output:
        "result/" + SAMPLE_NAME + "_sort_geneid.txt"
    params:
        count_sort_annotate_geneid = config["count_sort_annotate_geneid"]
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.count_sort_annotate_geneid} {input.count_table} {input.GB_matrix} {output}  
        """


#
# PART III
# Infer expressed genes of subclones
#

rule merge_bam_clones:
    input:
        bam = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        bai = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai", cell=BAMFILE)
    output:
        bam = expand("bam_merge/{clone}.merge.bam", clone=CLONE_NAME),
        bai = expand("bam_merge/{clone}.merge.bam.bai", clone=CLONE_NAME),
    shell:
        """
        perl utils/merge_bam_clones.pl
        """


rule count_reads_for_DNN:
    input:
        bam = expand("bam_merge/{clone}.merge.bam", clone=CLONE_NAME),
        bai = expand("bam_merge/{clone}.merge.bam.bai", clone=CLONE_NAME)
    output:
        tab = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + ".tab",
        npz = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + ".npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/bin_Genes_for_CNN_sort.txt --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """


rule count_reads_for_DNN_sc:
    input:
        bam = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        bai = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai", cell=BAMFILE)
    output:
        tab = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc.tab",
        npz = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc.npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/bin_Genes_for_CNN_sort.txt --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """


rule count_reads_chr_length:
    input:
        bam = expand("bam_merge/{clone}.merge.bam", clone=CLONE_NAME),
        bai = expand("bam_merge/{clone}.merge.bam.bai", clone=CLONE_NAME)
    output:
        tab = "result/Deeptool_chr_length_" + SAMPLE_NAME + ".tab",
        npz = "result/Deeptool_chr_length_" + SAMPLE_NAME + ".npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/bin_chr_length.bed --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """


rule count_reads_chr_length_sc:
    input:
        bam = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        bai = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai", cell=BAMFILE)
    output:
        tab = "result/Deeptool_chr_length_" + SAMPLE_NAME + "_sc.tab",
        npz = "result/Deeptool_chr_length_" + SAMPLE_NAME + "_sc.npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/bin_chr_length.bed --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """


rule count_reads_for_DNN_sort:
    input:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + ".tab"
    output:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_reads_for_DNN_sort_lab:
    input:
        count_reads_sort = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sort.txt",
        Ref_bed = "utils/bin_Genes_for_CNN_num_sort.txt",
    output:
        count_reads_sort_label = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sort_lab.txt"
    params:
        count_sort_label = config["count_sort_label"]
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.count_sort_label} {input.count_reads_sort} {input.Ref_bed} {output.count_reads_sort_label}
        """


rule count_reads_for_DNN_sort_label_sort:
    input:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sort_lab.txt"
    output:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sort_lab_final.txt"
    shell:
        """
        sort -k4,4n -t$'\t' {input} > {output}
		"""



rule count_reads_for_DNN_normalization:
    input:
        count_reads_chr_length = "result/Deeptool_chr_length_" + SAMPLE_NAME + ".tab",
        count_reads_sort_label = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sort_lab.txt",
        CNN_features_annot = "utils/bin_Genes_for_CNN_reshape_annot.txt",
        table_CpG = "utils/Features_reshape_CpG_orientation.txt",
        table_GC = "utils/Features_reshape_GC_orientation.txt",
        table_size = "utils/Features_reshape_size_orientation.txt",
        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",
        FPKM = "utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
        CN_result_data1 = "result/Features_reshape_{clone}_orientation_CN_correct0.txt",
    output:
        plot = "result/Features_reshape_" + SAMPLE_NAME + "_{clone}_orientation_norm_qc.pdf",
        table_mononuc_norm_data1 = "result/Features_reshape_{clone}_orientation_norm.txt",
    params:
        count_norm = config["count_norm"]
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.count_norm} {input.count_reads_chr_length} {input.count_reads_sort_label} {input.CNN_features_annot} {input.table_CpG} {input.table_GC} {input.table_size} {input.TSS_matrix} {input.FPKM} {input.CN_result_data1} {output.plot} {output.table_mononuc_norm_data1}
        """

rule count_reads_for_DNN_sc_sort:
    input:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc.tab"
    output:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


rule count_reads_for_DNN_sc_sort_lab:
    input:
        count_reads_sort = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort.txt",
        Ref_bed = "utils/bin_Genes_for_CNN_num_sort.txt",
    output:
        count_reads_sort_label = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort_lab.txt"
    params:
        count_sort_label = config["count_sort_label"]
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.count_sort_label} {input.count_reads_sort} {input.Ref_bed} {output.count_reads_sort_label}
        """


rule count_reads_for_DNN_sc_sort_label_sort:
    input:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort_lab.txt"
    output:
        "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort_lab_final.txt"
    shell:
        """
        sort -k4,4n -t$'\t' {input} > {output}
        """



rule generate_feature_sc_var:
    input:
        subclone_list = "input_user/input_subclonality.txt",
        count_reads_sc_sort = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort_lab_final.txt",
        Ref_bed_annot = "utils/bin_Genes_for_CNN_num_sort_ann_sort_GC_ensemble.txt",
        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",    
        CNN_features_annot = "utils/bin_Genes_for_CNN_reshape_annot.txt",
        FPKM = "utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
        CN_result_data1 = "result/Features_reshape_{clone}_orientation_CN_correct0.txt",
    output:
        plot = "result/Features_reshape_" + SAMPLE_NAME + "_{clone}_Resid_orientation_qc.pdf",
        table_mononuc_var_data1 = "result/Features_reshape_{clone}_Resid_orientation.txt",
    params:
        feature_sc_var = config["feature_sc_var"]
    log:
        "log/generate_feature_sc_var_{clone}.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.feature_sc_var} {input.subclone_list} {input.count_reads_sc_sort} {input.Ref_bed_annot} {input.TSS_matrix} {input.CNN_features_annot} {input.FPKM} {input.CN_result_data1} {output.plot} {output.table_mononuc_var_data1} > {log} 2>&1
       """


rule combine_features:
    input:
        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",
        table_GC_imput = "utils/Features_reshape_GC_orientation_impute.txt",
        table_CpG_imput = "utils/Features_reshape_CpG_orientation_impute.txt",
        table_RT = "utils/Features_reshape_RT_orientation.txt",
        table_mononuc_norm_data1 = "result/Features_reshape_{clone}_orientation_norm.txt",
        CN_result_data1 = "result/Features_reshape_{clone}_orientation_CN_correct0.txt",
        table_mononuc_var_data1 = "result/Features_reshape_{clone}_Resid_orientation.txt",
        FPKM = "utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
    output:
        features = "result/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_{clone}.txt",
        exp = "result/Expression_all_{clone}.txt",
        TSS_annot = "result/Features_reshape_all_TSS_matrix_woM_all_RT_{clone}.txt",
    params:
        combine_features = config["combine_features"]
    log:
        "log/combine_features_{clone}.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.combine_features} {input.TSS_matrix} {input.table_GC_imput} {input.table_CpG_imput} {input.table_RT} {input.table_mononuc_norm_data1} {input.CN_result_data1} {input.table_mononuc_var_data1} {input.FPKM} {output.features} {output.exp} {output.TSS_annot}
        """


rule infer_expressed_genes:
    input:
        features = "result/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_{clone}.txt",
        TSS_annot = "result/Features_reshape_all_TSS_matrix_woM_all_RT_{clone}.txt",
    output:
        train80 = "result_CNN/DNN_train80_output_ypred_{clone}.csv",
        train40 = "result_CNN/DNN_train40_output_ypred_{clone}.csv",
        train20 = "result_CNN/DNN_train20_output_ypred_{clone}.csv",
        train5 = "result_CNN/DNN_train5_output_ypred_{clone}.csv",
    log:
        "log/infer_expressed_genes_{clone}.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        module load foss/2019b
        module load Python/3.7.4-GCCcore-8.3.0
        module load cuDNN/7.6.4.38-gcccuda-2019b
        module load CUDA/10.1.243-GCC-8.3.0
        module load TensorFlow/1.15.0-fosscuda-2019b-Python-3.7.4
        module load scikit-learn/0.21.3-foss-2019b-Python-3.7.4
        module load matplotlib/3.1.1-foss-2019b-Python-3.7.4
        python utils/Deeplearning_Nucleosome_predict_train_RPE.py {input.features} {input.TSS_annot} {output.train80} {output.train40} {output.train20} {output.train5}
        """


rule annot_expressed_genes:
    input:
        TSS_annot = "result/Features_reshape_all_TSS_matrix_woM_all_RT_{clone}.txt",
        train80 = "result_CNN/DNN_train80_output_ypred_{clone}.csv",
        train40 = "result_CNN/DNN_train40_output_ypred_{clone}.csv",
        train20 = "result_CNN/DNN_train20_output_ypred_{clone}.csv",
        train5 = "result_CNN/DNN_train5_output_ypred_{clone}.csv",
    output:
        train80_annot = "result_CNN/DNN_train80_output_ypred_{clone}_annot.txt",
        train40_annot = "result_CNN/DNN_train40_output_ypred_{clone}_annot.txt",
        train20_annot = "result_CNN/DNN_train20_output_ypred_{clone}_annot.txt",
        train5_annot = "result_CNN/DNN_train5_output_ypred_{clone}_annot.txt",
    params:
        annot_expressed = config["annot_expressed"]
    log:
        "log/annot_expressed_genes_{clone}.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.annot_expressed} {input.TSS_annot} {input.train80} {input.train40} {input.train20} {input.train5} {output.train80_annot} {output.train40_annot} {output.train20_annot} {output.train5_annot}
        """

#
# PART IV
# DE analysis of epigenome of subclones
#

rule infer_differential_gene_expression:
    input:
        Genebody_NO = "result/" + SAMPLE_NAME + "_sort.txt",
        clonality = "input_user/input_subclonality.txt",
        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",
        GB_matrix = "utils/Strand_seq_matrix_Genebody_for_SCDE.txt",
        CNN_result1 = "result_CNN/DNN_train80_output_ypred_clone1_annot.txt",
        CNN_result2 = "result_CNN/DNN_train80_output_ypred_clone2_annot.txt",
        input_matrix = "input_user/input_SV_affected_genes.txt",
    output:
        "result_plots/Result_scNOVA_plots_" + SAMPLE_NAME + ".pdf",
    params:
        infer_diff_gene_expression = config["infer_diff_gene_expression"]
    log:
        "log/infer_diff_gene_expression.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.infer_diff_gene_expression} {input.Genebody_NO} {input.clonality} {input.TSS_matrix} {input.GB_matrix} {input.CNN_result1} {input.CNN_result2} {input.input_matrix} {output}
        """


#
# PART V
# Infer Single-cell TF motif accessibility using chromVAR
#
rule count_reads_CREs:
    input:
        bam = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam", cell=BAMFILE),
        bai = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.bai", cell=BAMFILE)
    output:
        tab = "result/" + SAMPLE_NAME + "_CREs_2kb.tab",
        npz = "result/" + SAMPLE_NAME + "_CREs_2kb.npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/regions_all_hg38_v2_resize_2kb_sort.bed --bamfiles {input.bam} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """

rule count_sort_by_coordinate_CREs:
    input:
        "result/" + SAMPLE_NAME + "_CREs_2kb.tab"
    output:
        "result/" + SAMPLE_NAME + "_CREs_2kb_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """

rule count_sort_annotate_chrid_CREs:
    input:
        "result/" + SAMPLE_NAME + "_CREs_2kb_sort.txt"
    output:
        "result/" + SAMPLE_NAME + "_CREs_2kb_sort_num.txt"
    params:
        count_sort_annotate_chrid_CREs = config["count_sort_annotate_chrid_CREs"]
    log:
        "log/count_sort_annotate_chrid_CREs.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.count_sort_annotate_chrid_CREs} {input} {output} 
        """

rule count_sort_annotate_chrid_CREs_sort:
    input:
        "result/" + SAMPLE_NAME + "_CREs_2kb_sort_num.txt"
    output:
        "result/" + SAMPLE_NAME + "_CREs_2kb_sort_num_sort_for_chromVAR.txt"
    shell:
        """
        sort -k1,1n -k2,2n -k3,3n -t$'\t' {input} > {output}
        """

#rule generate_CREs_CN:
#    input:
#        Result_enh_prom = "utils/Result_enh_prom_sort_new_num_sort_for_chromVAR.txt",
#        DHS_matrix_resize = "utils/regions_all_hg38_v2_resize_2kb_sort_num_sort_for_chromVAR.bed",
#        sv_calls  = "input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE.txt",
#    output:
#        "result/SC_CN_regions_all_hg38_v2_resize_2kb_sc.txt",
#    params:
#        generate_CREs_CN = config["generate_CREs_CN"]
#    log:
#        "log/generate_CREs_CN.log"
#    shell:
#        """
#        Rscript {params.generate_CREs_CN} {input.Result_enh_prom} {input.DHS_matrix_resize} {input.sv_calls} {output}
#        """

rule NO_chromVAR:
    input:
        DHS_annot_resize = "utils/Result_enh_prom_sort_new_num_sort_for_chromVAR.txt",
        DHS_matrix_resize = "utils/regions_all_hg38_v2_resize_2kb_sort_num_sort_for_chromVAR.bed",
        SC_CN_DHS = "input_user/sv_calls_all_print_CREs.txt",
        data1_resize = "result/" + SAMPLE_NAME + "_CREs_2kb_sort_num_sort_for_chromVAR.txt",
        class_label = "input_user/input_subclonality.txt",
    output:
        motif_zscore = "result/motif_dev_zscore_chromVAR_DHS_2kb_Enh_" + SAMPLE_NAME + ".txt",
        motif_dev = "result/motif_dev_chromVAR_DHS_2kb_Enh_" + SAMPLE_NAME + ".Rda",
        motif_variability = "result/motif_variability_chromVAR_DHS_2kb_Enh_" + SAMPLE_NAME + ".txt",
    params:
        NO_chromVAR = config["NO_chromVAR"]
    log:
        "log/NO_chromVAR.log"
    conda: "envs/scNOVA.yaml"
    shell:
        """
        Rscript {params.NO_chromVAR} {input.DHS_annot_resize} {input.DHS_matrix_resize} {input.SC_CN_DHS} {input.data1_resize} {input.class_label} {output.motif_zscore} {output.motif_dev} {output.motif_variability} 
        """




#
# PART VI
# Infer haplotype-resolved genebody NO
#
rule split_bam_WC:
    input:
        "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam"
    output:
        bam_header = "bam/{cell}.header_WC.sam",
        bam_C1 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C1.bam",
        bam_C2 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C2.bam",
        bam_W1 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W1.bam",
        bam_W2 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W2.bam",
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools view -H {input} > {output.bam_header}
        samtools view -f 99 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_C1}
        samtools view -f 147 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_C2}
        samtools view -f 83 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_W1}
        samtools view -f 163 {input} | cat {output.bam_header} - | samtools view -Sb - > {output.bam_W2}
        """

rule split_bam_WC_merge:
    input:
        bam_C1 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C1.bam",
        bam_C2 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C2.bam",
        bam_W1 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W1.bam",
        bam_W2 = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W2.bam",
    output:
        bam_C = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam",
        bam_W = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam",
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools merge {output.bam_C} {input.bam_C1} {input.bam_C2}
        samtools merge {output.bam_W} {input.bam_W1} {input.bam_W2}
        """

rule split_bam_WC_index:
    input:
        bam_C = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam",
        bam_W = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam",
    output:
        bam_C_ind = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam.bai",
        bam_W_ind = "bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam.bai",
    shell:
        """
        module load SAMtools/1.3.1-foss-2016b
        samtools index {input.bam_C}
        samtools index {input.bam_W}
        """

rule perl_split_sc:
    input:
        strandphaser_output = "input_user/strandphaser_output.txt",
        bam_C_ind = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.C.bam.bai", cell=BAMFILE),
        bam_W_ind = expand("bam/{cell}.sc_pre_mono_sort_for_mark_uniq.bam.W.bam.bai", cell=BAMFILE),
    output:
        nucleosome_sampleA = "nucleosome_sampleA/result.H1.bam",
        nucleosome_sampleB = "nucleosome_sampleB/result.H2.bam",
        perl_split_output = "input_user/strandphaser_output_copy.txt",
    shell:
        """
        perl utils/perl_test_all_snake.pl
        """

rule count_reads_CREs_haplo:
    input:
        bam1 = "nucleosome_sampleA/result.H1.bam",
        bam2 = "nucleosome_sampleB/result.H2.bam",
    output:
        tab = "result_haplo/Deeptool_DHS_2kb_H1H2.tab",
        npz = "result_haplo/Deeptool_DHS_2kb_H1H2.npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/regions_all_hg38_v2_resize_2kb_sort.bed --bamfiles {input.bam1} {input.bam2} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """

rule count_reads_CREs_haplo_sort_by_coordinate:
    input:
        "result_haplo/Deeptool_DHS_2kb_H1H2.tab"
    output:
        "result_haplo/Deeptool_DHS_2kb_H1H2_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """

rule count_reads_genebody_haplo:
    input:
        bam1 = "nucleosome_sampleA/result.H1.bam",
        bam2 = "nucleosome_sampleB/result.H2.bam",
    output:
        tab = "result_haplo/Deeptool_Genebody_H1H2.tab",
        npz = "result_haplo/Deeptool_Genebody_H1H2.npz",
    shell:
        """
        module load deeptools/2.5.1-foss-2016b-Python-2.7.12
        multiBamSummary BED-file --BED utils/bin_Genebody_all.bed --bamfiles {input.bam1} {input.bam2} \
            --extendReads --outRawCounts {output.tab} -out {output.npz}
        """

rule count_reads_genebody_haplo_sort_by_coordinate_genebody:
    input:
        "result_haplo/Deeptool_Genebody_H1H2.tab"
    output:
        "result_haplo/Deeptool_Genebody_H1H2_sort.txt"
    shell:
        """
        sort -k1,1 -k2,2n -k3,3n -t$'\t' {input} > {output}
        """


#
# PART V
# Single-cell level inference of expression
#

#rule generate_feature_sc:
#    input:
#        Deeptool_result_final = "result/Deeptool_Genes_for_CNN_" + SAMPLE_NAME + "_sc_sort_lab_final.txt",
#        CNN_features_annot = "utils/bin_Genes_for_CNN_reshape_annot.txt",
#        table_CpG = "utils/Features_reshape_CpG_orientation.txt",
#        table_GC = "utils/Features_reshape_GC_orientation.txt",
#        table_size = "utils/Features_reshape_size_orientation.txt",
#        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",
#        Deeptool_mapped = "result/Deeptool_chr_length_" + SAMPLE_NAME + "_sc.tab",
#        FPKM = "utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
#    output:
#        "result_features_sc/CNN_preprocessing_" + SAMPLE_NAME + "_sc.pdf",
#    params:
#        generate_feature_sc = config["generate_feature_sc"]
#    log:
#        "log/generate_feature_sc.log"
#    shell:
#        """
#        Rscript {params.generate_feature_sc} {input.Deeptool_result_final} {input.CNN_features_annot} {input.table_CpG} {input.table_GC} {input.table_size} {input.TSS_matrix} {input.Deeptool_mapped} {input.FPKM} {output}
#        """


#rule generate_feature_sc_CN:
#    input:
#        Deeptool_result_final_merge = "utils/Deeptool_Genes_for_CNN_merge_sort_lab_final.txt",
#        sv_calls  = "input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE.txt",
#        CNN_features_annot = "utils/bin_Genes_for_CNN_reshape_annot.txt",
#    output:
#        "result_features_sc_CN/SC_CN_DHS_for_CNN_" + SAMPLE_NAME + ".txt",
#    params:
#        generate_feature_sc_CN = config["generate_feature_sc_CN"]
#    log:
#        "log/generate_feature_sc_CN.log"
#    shell:
#        """
#        Rscript {params.generate_feature_sc_CN} {input.Deeptool_result_final_merge} {input.sv_calls} {input.CNN_features_annot} {output} 
#        """

#rule combine_feature_sc:
#    input:
#        TSS_matrix = "utils/Strand_seq_matrix_TSS_for_SVM.txt",
#        subclonality = "input_user/input_subclonality.txt",
#        table_GC_imput = "utils/Features_reshape_GC_orientation_impute.txt",
#        table_CpG_imput = "utils/Features_reshape_CpG_orientation_impute.txt",
#        table_RT = "utils/Features_reshape_RT_orientation.txt",
#        FPKM = "utils/FPKM_sort_LCL_RPE_19770_renamed.txt",
#    output:
#        "result_features_sc/Features_reshape_all_TSS_matrix_woM_all_sc.txt",
#    params:
#        combine_feature_sc = config["combine_feature_sc"]
#    log:
#        "log/combine_feature_sc.log"
#    shell:
#        """
#        Rscript {params.combine_feature_sc} {input.TSS_matrix} {input.subclonality} {input.table_GC_imput} {input.table_CpG_imput} {input.table_RT} {input.FPKM} {output} 
#        """


#rule merge_feature_sc:
#    input:
#        feature_sc = expand("result_features_sc/Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_{single_cells}_wovar_exp.txt", single_cells=BAM_SC)
#    output:
#        "result_features_sc/" + SAMPLE_NAME + "_wovar_exp.txt"
#    shell:
#        """
#        cat {input} > {output} 
#        """

