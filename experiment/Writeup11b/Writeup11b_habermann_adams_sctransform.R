rm(list=ls())
load("../../../../out/Writeup11b/Writeup11b_adams_T_preprocessed.RData")
load("../../../../out/Writeup11b/Writeup11b_habermann_T_preprocessed.RData")

library(Seurat)
library(eSVD2)
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

set.seed(10)
adams <- Seurat::SCTransform(adams,
                             residual.features = adams[["RNA"]]@var.features)
mat <- as.matrix(Matrix::t(adams[["SCT"]]@scale.data[adams[["SCT"]]@var.features,]))
covariate_dat <- as.data.frame(adams@meta.data[,c("percent.mt", "Disease_Identity", "Subject_Identity",
                                                  "Gender", "Age", "Tobacco")])
covariate_dat$Disease_Identity <- as.factor(covariate_dat$Disease_Identity)
covariate_dat$Subject_Identity <- as.factor(covariate_dat$Subject_Identity)
covariate_dat$Gender <- as.factor(covariate_dat$Gender)
covariate_dat$Tobacco <- as.factor(covariate_dat$Tobacco)
for(j in 1:ncol(mat)){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  tmp_df <- cbind(mat[,j], covariate_dat)
  colnames(tmp_df)[1] <- "y"
  lmer_res <- lme4::lmer(y ~ Disease_Identity + percent.mt + Gender + Age + Tobacco + (1 | Subject_Identity), data = tmp_df)
  pred_vec <- stats::predict(lmer_res)
  mat[,j] <- mat[,j] - pred_vec
}

case_idx <- which(adams$Disease_Identity == "IPF")
control_idx <- which(adams$Disease_Identity == "Control")
adams_wilcoxon_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  res <- stats::wilcox.test(x = mat[case_idx,j],
                            y = mat[control_idx,j])
  if(mean(mat[case_idx,j]) < mean(mat[control_idx,j])){
    res$statistic
  } else {
    -res$statistic
  }
})
names(adams_wilcoxon_vec) <- colnames(mat)

##

set.seed(10)
habermann <- Seurat::SCTransform(habermann,
                                 residual.features = habermann[["RNA"]]@var.features)
mat <- as.matrix(Matrix::t(habermann[["SCT"]]@scale.data[habermann[["SCT"]]@var.features,]))
covariate_dat <- as.data.frame(habermann@meta.data[,c("Diagnosis", "Sample_Name",
                                        "percent.mt", "Gender", "Age", "Tobacco")])
covariate_df <- data.frame(covariate_dat)
covariate_df$Gender <- as.factor(covariate_df$Gender)
covariate_df$Diagnosis <- as.factor(covariate_df$Diagnosis)
covariate_df$Tobacco <- as.factor(covariate_df$Tobacco)
covariate_df$Sample_Name <- as.factor(covariate_df$Sample_Name)
for(j in 1:ncol(mat)){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  tmp_df <- cbind(mat[,j], covariate_dat)
  colnames(tmp_df)[1] <- "y"
  lmer_res <- lme4::lmer(y ~ Diagnosis + percent.mt + Gender + Age + Tobacco + (1 | Sample_Name), data = tmp_df)
  pred_vec <- stats::predict(lmer_res)
  mat[,j] <- mat[,j] - pred_vec
}

case_idx <- which(habermann$Diagnosis == "IPF")
control_idx <- which(habermann$Diagnosis == "Control")
habermann_wilcoxon_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  res <- stats::wilcox.test(x = mat[case_idx,j],
                            y = mat[control_idx,j])
  if(mean(mat[case_idx,j]) < mean(mat[control_idx,j])){
    res$statistic
  } else {
    -res$statistic
  }
})
names(habermann_wilcoxon_vec) <- colnames(mat)


adams_wilcoxon_vec2 <- rep(0, length(all_genes))
for(i in 1:length(adams_wilcoxon_vec2)){
  gene_name <- all_genes[i]
  idx <- which(names(adams_wilcoxon_vec) == gene_name)
  if(length(idx) == 1) {
    adams_wilcoxon_vec2[i] <- sign(adams_wilcoxon_vec[idx])*log(abs(adams_wilcoxon_vec[idx]))
  }
}
names(adams_wilcoxon_vec2) <- names(adams_wilcoxon_vec)

habermann_wilcoxon_vec2 <- rep(0, length(all_genes))
for(i in 1:length(habermann_wilcoxon_vec2)){
  gene_name <- all_genes[i]
  idx <- which(names(habermann_wilcoxon_vec) == gene_name)
  if(length(idx) == 1) {
    habermann_wilcoxon_vec2[i] <- sign(habermann_wilcoxon_vec[idx])*log(abs(habermann_wilcoxon_vec[idx]))
  }
}
names(habermann_wilcoxon_vec2) <- names(habermann_wilcoxon_vec)

#############################

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

xlim <- quantile(habermann_wilcoxon_vec2, probs = c(0.01, 0.99))*1.1
ylim <- quantile(adams_wilcoxon_vec2, probs = c(0.01, 0.99))*1.1

col_template_vec <- c(2,4,3)
idx_list <- lapply(gene_list, function(gene_vec){
  which(names(habermann_wilcoxon_vec2) %in% gene_vec)
})
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(habermann_wilcoxon_vec2))
for(i in 1:length(idx_list)){
  col_vec[idx_list[[i]]] <- col_template_vec[i]
}


png("../../../../out/fig/Writeup11b/habermann_adams_teststatistic_sctransformWilcoxon.png",
    height = 2500, width = 2500, res = 300, units = "px")
plot(habermann_wilcoxon_vec2, adams_wilcoxon_vec2,
     xlab = "Habermann", ylab = "Adams", pch = 16, col = col_vec,
     xlim = xlim, ylim = ylim,
     main = "T-cells\nUsing SCTransform and Wilcoxon")
graphics.off()

png("../../../../out/fig/Writeup11b/habermann_adams_teststatistic_sctransformWilcoxon_separate.png",
    height = 1200, width = 3600, res = 300, units = "px")
par(mfrow = c(1,3))
for(i in 1:length(idx_list)){
  plot(habermann_wilcoxon_vec2[idx_list[[i]]], adams_wilcoxon_vec2[idx_list[[i]]],
       xlab = "Habermann", ylab = "Adams", pch = 16, col = col_template_vec[i],
       xlim = xlim, ylim = ylim,
       main = "T-cells\nUsing SCTransform and Wilcoxon")
  lines(c(-15,15), rep(0,2), lwd = 2, lty = 2)
  lines(rep(0,2), c(-15,15), lwd = 2, lty = 2)
}
graphics.off()

