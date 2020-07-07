## 1) Load feature sets


##-------------------------------------------------------------------------------------
##Features: Benchmarking : Nucleosome occupancy (mono) 2K around TSS, NDR (-150 to 50bp)
##-------------------------------------------------------------------------------------
setwd('/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_SCDE')


TSS_matrix <- read.table("Strand_seq_matrix_TSS_for_SVM.txt", header=TRUE, sep ='\t')
NDR_matrix <- read.table("Strand_seq_matrix_NDR_for_SCDE.txt", header=TRUE, sep ='\t')



##-------------------------------------------------------------------------------------
##Features: Sequence : GC%, CpG%, RT
##-------------------------------------------------------------------------------------
setwd('/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_SCDE')
table_GC_imput <- read.table("Features_reshape_GC_orientation_impute.txt", header=F, sep ='\t', comment.char = "")
table_CpG_imput <- read.table("Features_reshape_CpG_orientation_impute.txt", header=F, sep ='\t', comment.char = "")
table_RT <- read.table("Features_reshape_RT_orientation.txt", header=F, sep ='\t', comment.char = "")

##-------------------------------------------------------------------------------------
##Features: Nucleosome occupancy 150 bins and copy-number normalization
##-------------------------------------------------------------------------------------
setwd('/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_SCDE')
table_mononuc_norm_data1 <- read.table("Features_reshape_BCLL01_P1P2_C0_orientation_norm.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_norm_data2 <- read.table("Features_reshape_BCLL01_P1P2_C1_orientation_norm.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_norm_data3 <- read.table("Features_reshape_BCLL01_P1P2_C2_orientation_norm.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_norm_data4 <- read.table("Features_reshape_BCLL01_P1P2_C3_orientation_norm.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_norm_data5 <- read.table("Features_reshape_BCLL01_P1P2_C4_orientation_norm.txt", header=F, sep ='\t', comment.char = "")

##Normalization by copy number (19757 X 150 copy number matrix)
setwd('/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_SCDE')
CN_result_data1 <- read.table("Features_reshape_BCLL01_P1P2_C0_orientation_CN_correct0.txt", sep = '\t', header=F)
CN_result_data2 <- read.table("Features_reshape_BCLL01_P1P2_C1_orientation_CN_correct0.txt", sep = '\t', header=F)
CN_result_data3 <- read.table("Features_reshape_BCLL01_P1P2_C2_orientation_CN_correct0.txt", sep = '\t', header=F)
CN_result_data4 <- read.table("Features_reshape_BCLL01_P1P2_C3_orientation_CN_correct0.txt", sep = '\t', header=F)
CN_result_data5 <- read.table("Features_reshape_BCLL01_P1P2_C4_orientation_CN_correct0.txt", sep = '\t', header=F)

table_mononuc_norm_data1_cn <- table_mononuc_norm_data1/CN_result_data1
table_mononuc_norm_data2_cn <- table_mononuc_norm_data2/CN_result_data2
table_mononuc_norm_data3_cn <- table_mononuc_norm_data3/CN_result_data3
table_mononuc_norm_data4_cn <- table_mononuc_norm_data4/CN_result_data4
table_mononuc_norm_data5_cn <- table_mononuc_norm_data5/CN_result_data5

##-------------------------------------------------------------------------------------
##Features: Nucleosome occupancy Residual of the CV square 150 bins
##-------------------------------------------------------------------------------------
setwd('/Users/jeong/Documents/Strand_Seq/Deeptool/deeptool_SCDE')
table_mononuc_var_data1 <- read.table("Features_reshape_BCLL01_P1P2_C0_Resid_orientation.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_var_data2 <- read.table("Features_reshape_BCLL01_P1P2_C1_Resid_orientation.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_var_data3 <- read.table("Features_reshape_BCLL01_P1P2_C2_Resid_orientation.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_var_data4 <- read.table("Features_reshape_BCLL01_P1P2_C3_Resid_orientation.txt", header=F, sep ='\t', comment.char = "")
table_mononuc_var_data5 <- read.table("Features_reshape_BCLL01_P1P2_C4_Resid_orientation.txt", header=F, sep ='\t', comment.char = "")


## 2) Load gene expression information (Expressed / Unexpressed)

setwd('/Users/jeong/Documents/Strand_Seq/RNA_seq/htseq_DEanalysis_2019_renamed')
FPKM <- read.table("FPKM_sort_LCL_RPE_19770_renamed.txt", header=T, sep ='\t', comment.char = "")
TSS_matrix_woM <- TSS_matrix[TSS_matrix[,2]!="chrM",]
FPKM_woM <- FPKM[TSS_matrix[,2]!="chrM",]
FPKM_woM_LCL <- cbind(as.matrix(rowMeans(FPKM_woM[,1:9])), as.matrix(rowMeans(FPKM_woM[,1:9])), as.matrix(rowMeans(FPKM_woM[,1:9])), as.matrix(rowMeans(FPKM_woM[,1:9])))
Expression_label <- matrix(0, nrow(FPKM_woM_LCL), 4)
for (i in 1:4){
  Expression_label[FPKM_woM_LCL[,i]>1,i] <- 1
}



## 3) Generate feature sets and target vector (Expressed = 1 / Unexpressed = 0)


Features_label <- as.matrix(Expression_label[,1])
TSS_matrix_woM_all <- TSS_matrix_woM
table_RT_all <- table_RT
for (i in 1:4){
  Features_label <- rbind(Features_label,as.matrix(Expression_label[,(i+1)]))
  TSS_matrix_woM_all <- rbind(TSS_matrix_woM_all, TSS_matrix_woM)
  table_RT_all <- rbind(table_RT_all, as.matrix(table_RT))
}


##Alternative way to make input format for the transposed features (150X3)
Features_t1 <- cbind(table_mononuc_norm_data1_cn[,1], table_mononuc_var_data1[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
Features_t2 <- cbind(table_mononuc_norm_data2_cn[,1], table_mononuc_var_data2[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
Features_t3 <- cbind(table_mononuc_norm_data3_cn[,1], table_mononuc_var_data3[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
Features_t4 <- cbind(table_mononuc_norm_data4_cn[,1], table_mononuc_var_data4[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
Features_t5 <- cbind(table_mononuc_norm_data5_cn[,1], table_mononuc_var_data5[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
#Features_t6 <- cbind(table_mononuc_norm_data6_cn[,1], table_mononuc_var_data6[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
#Features_t7 <- cbind(table_mononuc_norm_data7_cn[,1], table_mononuc_var_data7[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
#Features_t8 <- cbind(table_mononuc_norm_data8_cn[,1], table_mononuc_var_data8[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
#Features_t9 <- cbind(table_mononuc_norm_data9_cn[,1], table_mononuc_var_data9[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
#Features_t10 <- cbind(table_mononuc_norm_data10_cn[,1], table_mononuc_var_data10[,1], table_GC_imput[,1], table_CpG_imput[,1], table_RT[,1]/100)
for (i in 2:150){
  Features_t1 <- cbind(Features_t1, cbind(table_mononuc_norm_data1_cn[,i], table_mononuc_var_data1[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  Features_t2 <- cbind(Features_t2, cbind(table_mononuc_norm_data2_cn[,i], table_mononuc_var_data2[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  Features_t3 <- cbind(Features_t3, cbind(table_mononuc_norm_data3_cn[,i], table_mononuc_var_data3[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  Features_t4 <- cbind(Features_t4, cbind(table_mononuc_norm_data4_cn[,i], table_mononuc_var_data4[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  Features_t5 <- cbind(Features_t5, cbind(table_mononuc_norm_data5_cn[,i], table_mononuc_var_data5[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  #Features_t6 <- cbind(Features_t6, cbind(table_mononuc_norm_data6_cn[,i], table_mononuc_var_data6[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  #Features_t7 <- cbind(Features_t7, cbind(table_mononuc_norm_data7_cn[,i], table_mononuc_var_data7[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  #Features_t8 <- cbind(Features_t8, cbind(table_mononuc_norm_data8_cn[,i], table_mononuc_var_data8[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  #Features_t9 <- cbind(Features_t9, cbind(table_mononuc_norm_data9_cn[,i], table_mononuc_var_data9[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
  #Features_t10 <- cbind(Features_t10, cbind(table_mononuc_norm_data10_cn[,i], table_mononuc_var_data10[,i], table_GC_imput[,i], table_CpG_imput[,i], table_RT[,i]/100))
}
Features_tall <- rbind(Features_t1, Features_t2, Features_t3, Features_t4, Features_t5)
Features_both_sub <- Features_tall


standard_svm <- data.frame(Features_both_sub)
standard_svm$Type<-rep(0, nrow(standard_svm))
standard_svm[Features_label=="0",ncol(Features_both_sub)+1] <- 0
standard_svm[Features_label=="1",ncol(Features_both_sub)+1] <- 1
standard_svm$Type <- as.factor(standard_svm$Type)

standard_svm_RT <- standard_svm[is.na(rowSums(table_RT_all))==0,]
TSS_matrix_woM_all_RT <- TSS_matrix_woM_all[is.na(rowSums(table_RT_all))==0,]


##This is to practice leave one chromosome validation
write.table(standard_svm_RT[,1:750], "Features_reshape_all_orientation_norm_var_GC_CpG_RT_T_comb3_BCLL01_P1P2.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
write.table(standard_svm_RT[,751], "Expression_all_BCLL01_P1P2.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
write.table(TSS_matrix_woM_all_RT[,c(2:5, 60:61)], "Features_reshape_all_TSS_matrix_woM_all_RT_BCLL01_P1P2.txt", row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
