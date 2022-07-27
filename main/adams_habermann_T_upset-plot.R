rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/adams_T_sctransform.RData")
sctransform_adams <- de_result
load("../../../out/main/habermann_T_sctransform.RData")
sctransform_habermann <- de_result

load("../../../out/main/habermann_T_esvd.RData")
eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj$teststat_vec <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T)
metadata <- habermann@meta.data
metadata[,"Sample_Name"] <- as.factor(metadata[,"Sample_Name"])
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Sample_Name",
                                           metadata = metadata,
                                           verbose = 1)
eSVD_obj_habermann <- eSVD_obj

load("../../../out/main/adams_T_esvd.RData")
eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj$teststat_vec <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T)
metadata <- adams@meta.data
metadata[,"Subject_Identity"] <- as.factor(metadata[,"Subject_Identity"])
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "Subject_Identity",
                                           metadata = metadata,
                                           verbose = 1)
eSVD_obj_adams <- eSVD_obj

############

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
de_genes_others <- unique(c(adams_df_genes_others, habermann_df_genes_others))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
de_genes_others <- setdiff(de_genes_others, de_genes)
cycling_genes <- setdiff(cycling_genes, c(de_genes_others, de_genes))
hk_genes <- setdiff(hk_genes, c(cycling_genes, de_genes_others, de_genes))

#####################

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_adams,
                             metadata = adams@meta.data,
                             covariate_individual = "Subject_Identity")
teststat_vec <- eSVD_obj_adams$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_adams <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_adams <- locfdr::locfdr(gaussian_teststat_adams, plot = 0)
fdr_vec <- locfdr_res_adams$fdr
names(fdr_vec) <- names(gaussian_teststat_adams)
null_mean <- locfdr_res_adams$fp0["mlest", "delta"]
null_sd <- locfdr_res_adams$fp0["mlest", "sigma"]
pvalue_vec <- sapply(gaussian_teststat_adams, function(x){
  if(x < null_mean) {
    stats::pnorm(x, mean = null_mean, sd = null_sd)*2
  } else {
    (1-stats::pnorm(x, mean = null_mean, sd = null_sd))*2
  }
})
logpvalue_vec <- -log10(pvalue_vec)
idx_adams <- order(logpvalue_vec, decreasing = T)[1:length(unique(c(adams_df_genes, habermann_df_genes)))]
eSVD_adams_de <- names(teststat_vec)[idx_adams]

####

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_habermann,
                             metadata = habermann@meta.data,
                             covariate_individual = "Sample_Name")
teststat_vec <- eSVD_obj_habermann$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_habermann <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_habermann <- locfdr::locfdr(gaussian_teststat_habermann, plot = 0)
fdr_vec <- locfdr_res_habermann$fdr
names(fdr_vec) <- names(gaussian_teststat_habermann)
null_mean <- locfdr_res_habermann$fp0["mlest", "delta"]
null_sd <- locfdr_res_habermann$fp0["mlest", "sigma"]
pvalue_vec <- sapply(gaussian_teststat_habermann, function(x){
  if(x < null_mean) {
    stats::pnorm(x, mean = null_mean, sd = null_sd)*2
  } else {
    (1-stats::pnorm(x, mean = null_mean, sd = null_sd))*2
  }
})
logpvalue_vec <- -log10(pvalue_vec)
idx_habermann <- order(logpvalue_vec, decreasing = T)[1:length(unique(c(adams_df_genes, habermann_df_genes)))]
eSVD_habermann_de <- names(teststat_vec)[idx_habermann]

###############################

sctransform_adams_degenes <- rownames(sctransform_adams)[order(sctransform_adams[,"p_val"], decreasing = F)[1:length(unique(c(adams_df_genes, habermann_df_genes)))]]
sctransform_habermann_degenes <- rownames(sctransform_habermann)[order(sctransform_habermann[,"p_val"], decreasing = F)[1:length(unique(c(adams_df_genes, habermann_df_genes)))]]

de_vec <- unique(c(adams_df_genes, habermann_df_genes))

length(de_vec)
length(adams_df_genes)
length(habermann_df_genes)
length(intersect(adams_df_genes, habermann_df_genes))

length(eSVD_adams_de)
length(intersect(eSVD_adams_de, habermann_df_genes))
length(intersect(eSVD_adams_de, adams_df_genes))
length(intersect(eSVD_adams_de, de_vec))
length(intersect(eSVD_adams_de, hk_genes))

length(eSVD_habermann_de)
length(intersect(eSVD_habermann_de, adams_df_genes))
length(intersect(eSVD_habermann_de, habermann_df_genes))
length(intersect(eSVD_habermann_de, de_vec))
length(intersect(eSVD_habermann_de, hk_genes))

length(intersect(eSVD_adams_de, eSVD_habermann_de))

length(sctransform_adams_degenes)
length(intersect(sctransform_adams_degenes, habermann_df_genes))
length(intersect(sctransform_adams_degenes, adams_df_genes))
length(intersect(sctransform_adams_degenes, de_vec))
length(intersect(sctransform_adams_degenes, hk_genes))

length(sctransform_habermann_degenes)
length(intersect(sctransform_habermann_degenes, adams_df_genes))
length(intersect(sctransform_habermann_degenes, habermann_df_genes))
length(intersect(sctransform_habermann_degenes, de_vec))
length(intersect(sctransform_habermann_degenes, hk_genes))

length(intersect(sctransform_adams_degenes, sctransform_habermann_degenes))



