rm(list=ls())
load("../../../../out/Writeup11b/Writeup11b_adams_T_preprocessed.RData")
load("../../../../out/Writeup11b/Writeup11b_habermann_T_preprocessed.RData")

library(Seurat)
library(DESeq2)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Disease_Identity %in% c("Control", "IPF"))] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
tab_vec <- table(adams$Subject_Identity)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Subject_Identity %in% subj_keep)] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
table(adams$Subject_Identity, adams$Disease_Identity)

adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]

keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$Diagnosis %in% c("Control", "IPF"))] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)
tab_vec <- table(habermann$Sample_Name)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$Sample_Name %in% subj_keep)] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)
table(habermann$Sample_Name, habermann$Diagnosis)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann <- Seurat::NormalizeData(habermann,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
habermann <- Seurat::FindVariableFeatures(habermann,
                                          selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
length(intersect(habermann_df_genes, adams_df_genes))
length(unique(c(habermann_df_genes, adams_df_genes)))
length(intersect(Seurat::VariableFeatures(habermann), Seurat::VariableFeatures(adams)))
length(unique(c(Seurat::VariableFeatures(habermann), Seurat::VariableFeatures(adams))))

hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

all_genes <- unique(c(Seurat::VariableFeatures(habermann),
                      Seurat::VariableFeatures(adams),
                      habermann_df_genes,
                      adams_df_genes,
                      hk_genes,
                      cycling_genes))
all_available_genes <- intersect(rownames(adams), rownames(habermann))
all_genes <- intersect(all_available_genes, all_genes)

adams[["RNA"]]@var.features <- all_genes
habermann[["RNA"]]@var.features <- all_genes

#######################

case_individuals <- unique(adams@meta.data[which(adams$Disease_Identity == "IPF"),"Subject_Identity"])
control_individuals <- unique(adams@meta.data[which(adams$Disease_Identity == "Control"),"Subject_Identity"])
mat <- as.matrix(adams[["RNA"]]@counts[adams[["RNA"]]@var.features,])
pseudomat_case <- sapply(case_individuals, function(indiv){
  idx <- which(adams$Subject_Identity == indiv)
  rowSums(mat[,idx])
})
pseudomat_control <- sapply(control_individuals, function(indiv){
  idx <- which(adams$Subject_Identity == indiv)
  rowSums(mat[,idx])
})
pseudomat <- cbind(pseudomat_case, pseudomat_control)
covariate <- as.data.frame(t(sapply(c(case_individuals, control_individuals), function(indiv){
  idx <- which(adams$Subject_Identity == indiv)[1]
  data.frame(Disease_Identity = adams$Disease_Identity[idx],
             Gender = adams$Gender[idx],
             Age = adams$Age[idx],
             Tobacco = adams$Tobacco[idx])
})))
covariate$Disease_Identity <- as.factor(unlist(covariate$Disease_Identity))
covariate$Gender <- as.factor(unlist(covariate$Gender))
covariate$Age <- scale(as.numeric(covariate$Age))
covariate$Tobacco <- as.factor(unlist(covariate$Tobacco))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudomat,
                                      colData = covariate,
                                      design = ~ Gender + Age + Tobacco + Disease_Identity)
dds <- DESeq2::DESeq(dds)
adams_res <- DESeq2::results(dds)

################

case_individuals <- unique(habermann@meta.data[which(habermann$Diagnosis == "IPF"),"Sample_Name"])
control_individuals <- unique(habermann@meta.data[which(habermann$Diagnosis == "Control"),"Sample_Name"])
mat <- as.matrix(habermann[["RNA"]]@counts[habermann[["RNA"]]@var.features,])
pseudomat_case <- sapply(case_individuals, function(indiv){
  idx <- which(habermann$Sample_Name == indiv)
  rowSums(mat[,idx])
})
pseudomat_control <- sapply(control_individuals, function(indiv){
  idx <- which(habermann$Sample_Name == indiv)
  rowSums(mat[,idx])
})
pseudomat <- cbind(pseudomat_case, pseudomat_control)
covariate <- as.data.frame(t(sapply(c(case_individuals, control_individuals), function(indiv){
  idx <- which(habermann$Sample_Name == indiv)[1]
  data.frame(Diagnosis = habermann$Diagnosis[idx],
             Gender = habermann$Gender[idx],
             Age = habermann$Age[idx],
             Tobacco = habermann$Tobacco[idx])
})))
covariate$Diagnosis <- as.factor(unlist(covariate$Diagnosis))
covariate$Gender <- as.factor(unlist(covariate$Gender))
covariate$Age <- scale(as.numeric(covariate$Age))
covariate$Tobacco <- as.factor(unlist(covariate$Tobacco))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = pseudomat,
                                      colData = covariate,
                                      design = ~ Gender + Age + Tobacco + Diagnosis)
dds <- DESeq2::DESeq(dds)
habermann_res <- DESeq2::results(dds)

#######################

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
adams_df_genes_others <- unique(df_mat$gene[which(df_mat$cellType %in% c("B", "Macrophage", "Macrophage Alveolar", "NK"))])
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X
file_vec <- c("Macrophages_disease_vs_control_.csv", "Monocytes_disease_vs_control_.csv",
              "B_Cells_disease_vs_control_.csv", "NK_Cells_disease_vs_control_.csv")
habermann_df_genes_others <- unique(unlist(lapply(file_vec, function(file_suffix){
  df_mat <- read.csv(paste0("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/", file_suffix),
                     sep = ",")
  df_mat$X
})))
de_genes <- unique(c(adams_df_genes, habermann_df_genes))
other_genes <- unique(c(adams_df_genes_others, habermann_df_genes_others))

hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_list <- list(de_genes,
                  setdiff(other_genes, de_genes),
                  setdiff(unique(c(hk_genes, cycling_genes)), c(other_genes, de_genes)))
names(gene_list) <- c("Published DE gene", "Other interest genes", "Housekeeping gene")

#######################

adams_teststat_vec <- adams_res[,"stat"]
adams_teststat_vec[is.na(adams_teststat_vec)] <- 0
habermann_teststat_vec <- habermann_res[,"stat"]
habermann_teststat_vec[is.na(habermann_teststat_vec)] <- 0
names(adams_teststat_vec) <- rownames(adams_res)
names(habermann_teststat_vec) <- rownames(habermann_res)
all(names(adams_teststat_vec) == names(habermann_teststat_vec))

xlim <- quantile(habermann_teststat_vec, probs = c(0.01, 0.99))*1.1
ylim <- quantile(adams_teststat_vec, probs = c(0.01, 0.99))*1.1

col_template_vec <- c(2,4,3)
idx_list <- lapply(gene_list, function(gene_vec){
  which(names(habermann_teststat_vec) %in% gene_vec)
})
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(habermann_teststat_vec))
for(i in 1:length(idx_list)){
  col_vec[idx_list[[i]]] <- col_template_vec[i]
}

png("../../../../out/fig/Writeup11b/habermann_adams_teststatistic_DESeq2.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(habermann_teststat_vec, adams_teststat_vec,
     xlab = "Habermann", ylab = "Adams", pch = 16, col = col_vec,
     xlim = xlim, ylim = ylim,
     main = "T-cells\nUsing DESeq2 on Pseudo-bulk")
graphics.off()

png("../../../../out/fig/Writeup11b/habermann_adams_teststatistic_DESeq2_separate.png",
    height = 1200, width = 3600, res = 300, units = "px")
par(mfrow = c(1,3))
for(i in 1:length(idx_list)){
  plot(habermann_teststat_vec[idx_list[[i]]], adams_teststat_vec[idx_list[[i]]],
       xlab = "Habermann", ylab = "Adams", pch = 16, col = col_template_vec[i],
       xlim = xlim, ylim = ylim,
       main = "T-cells\nUsing DESeq2 on Pseudo-bulk")
  lines(c(-15,15), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-15,15), lwd = 2, lty = 2)
}
graphics.off()



