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

habermann[["RNA"]]@var.features <- all_genes

#############################

mat <- as.matrix(Matrix::t(habermann[["RNA"]]@counts[habermann[["RNA"]]@var.features,]))
covariate_dat <- habermann@meta.data[,c("Diagnosis", "Sample_Name",
                                        "percent.mt", "Gender", "Age", "Tobacco")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Gender"] <- as.factor(covariate_df[,"Gender"])
covariate_df[,"Diagnosis"] <- factor(covariate_df[,"Diagnosis"], levels = c("Control", "IPF"))
covariate_df[,"Tobacco"] <- as.factor(covariate_df[,"Tobacco"])
covariate_df[,"Sample_Name"] <- as.factor(covariate_df[,"Sample_Name"])

covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        mixed_effect_variables = "Sample_Name")

case_control_variable <- "Diagnosis_IPF"
library_size_variable <- "Log_UMI"
mixed_effect_variables <- c(colnames(covariates)[grep("^Sample_Name", colnames(covariates))])

#############################
print("Initializing")
time_start1 <- Sys.time()
eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                    bool_intercept = T,
                                    covariates = covariates,
                                    case_control_variable = case_control_variable,
                                    k = 30,
                                    lambda = 0.1,
                                    mixed_effect_variables = mixed_effect_variables,
                                    offset_variables = "Log_UMI",
                                    verbose = 1)
time_end1 <- Sys.time()

pval_thres <- min(exp(quantile(eSVD_obj$initial_Reg$log_pval, probs = 0.05)), 1e-5)
# pval_thres <- 0.001
eSVD_obj <- eSVD2:::apply_initial_threshold(eSVD_obj = eSVD_obj,
                                            pval_thres = pval_thres,
                                            verbose = 1)

offset_variables <- setdiff(colnames(eSVD_obj$covariates), case_control_variable)
print("First fit")
time_start2 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = offset_variables,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_First",
                             fit_previous = "fit_Init")
time_end2 <- Sys.time()

print("Second fit")
time_start3 <- Sys.time()
eSVD_obj <- eSVD2:::opt_esvd(input_obj = eSVD_obj,
                             l2pen = 0.1,
                             max_iter = 100,
                             offset_variables = NULL,
                             tol = 1e-6,
                             verbose = 1,
                             fit_name = "fit_Second",
                             fit_previous = "fit_First")
time_end3 <- Sys.time()

print("Nuisance estimation")
time_start4 <- Sys.time()
eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                      verbose = 1)
time_end4 <- Sys.time()

save(date_of_run, session_info, habermann, covariate_df,
     eSVD_obj,
     time_start1, time_end1, time_start2, time_end2,
     time_start3, time_end3, time_start4, time_end4,
     file = "../../../../out/Writeup11d/Writeup11d_habermann_T_esvd.RData")
