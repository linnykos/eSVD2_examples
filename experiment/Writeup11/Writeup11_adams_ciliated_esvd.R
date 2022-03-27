rm(list=ls())
load("../../../../out/Writeup11/Writeup11_adams_epithelial_preprocessed.RData")
load("../../../../out/Writeup11/Writeup11_habermann_epithelial_preprocessed.RData")

library(Seurat)
library(eSVD2)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(0, ncol(adams))
keep_vec[intersect(which(adams$Manuscript_Identity == "Ciliated"),
                   which(adams$Disease_Identity %in% c("Control", "IPF")))] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)
tab_vec <- table(adams$Subject_Identity)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(adams))
keep_vec[which(adams$Subject_Identity %in% subj_keep)] <- 1
adams$keep <- keep_vec
adams <- subset(adams, keep == 1)

adams <- Seurat::NormalizeData(adams,
                               normalization.method = "LogNormalize", scale.factor = 10000)
adams <- Seurat::FindVariableFeatures(adams,
                                      selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "Ciliated")]

keep_vec <- rep(0, ncol(habermann))
keep_vec[intersect(which(habermann$celltype == "Ciliated"),
                   which(habermann$Diagnosis %in% c("Control", "IPF")))] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)
tab_vec <- table(habermann$Sample_Name)
subj_keep <- names(tab_vec)[which(tab_vec >= 50)]
keep_vec <- rep(0, ncol(habermann))
keep_vec[which(habermann$Sample_Name %in% subj_keep)] <- 1
habermann$keep <- keep_vec
habermann <- subset(habermann, keep == 1)

Seurat::DefaultAssay(habermann) <- "RNA"
habermann <- Seurat::NormalizeData(habermann,
                               normalization.method = "LogNormalize", scale.factor = 10000)
habermann <- Seurat::FindVariableFeatures(habermann,
                                      selection.method = "vst", nfeatures = 5000)

df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/Ciliated_disease_vs_control_.csv",
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

#############################

mat <- as.matrix(Matrix::t(adams[["RNA"]]@counts[adams[["RNA"]]@var.features,]))
covariate_dat <- adams@meta.data[,c("percent.mt", "Disease_Identity", "Subject_Identity",
                                    "Gender", "Age", "Ethnicity", "Tobacco")]
tmp <- as.character(covariate_dat$Ethnicity)
tmp[tmp != "white"] <- "nonwhite"
covariate_dat$Ethnicity <- as.factor(tmp)
covariate_df <- data.frame(covariate_dat)
covariate_df[,"Gender"] <- as.factor(covariate_df[,"Gender"])
covariate_df[,"Ethnicity"] <- as.factor(covariate_df[,"Ethnicity"])
covariate_df[,"Disease_Identity"] <- factor(covariate_df[,"Disease_Identity"], levels = c("Control", "IPF"))
covariate_df[,"Tobacco"] <- as.factor(covariate_df[,"Tobacco"])
covariate_df[,"Subject_Identity"] <- as.factor(covariate_df[,"Subject_Identity"])

covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        mixed_effect_variables = c("Subject_Identity"))

#####################
mixed_effect_variables <- colnames(covariates)[grep("^Subject_Identity_", colnames(covariates))]

time_start1 <- Sys.time()
esvd_init <- eSVD2:::initialize_esvd(dat = mat,
                                     covariates = covariates,
                                     case_control_variable = "Disease_Identity_IPF",
                                     k = 30,
                                     lambda = 0.1,
                                     mixed_effect_variables = mixed_effect_variables,
                                     offset_variables = "Log_UMI",
                                     verbose = 2,
                                     tmp_path = "../../../../out/Writeup11/Writeup11_adams_ciliated_esvd_tmp.RData")
time_end1 <- Sys.time()

save(date_of_run, session_info, adams, covariate_df,
     esvd_init, time_start1, time_end1,
     file = "../../../../out/Writeup11/Writeup11_adams_ciliated_esvd_tmp.RData")

#############

print("Starting first eSVD fit")

case_control_variable <- "Disease_Identity_IPF"
offset_var <- setdiff(colnames(esvd_init$covariates), case_control_variable)
offset_mat <- tcrossprod(esvd_init$covariates[,offset_var], esvd_init$b_mat[,offset_var])
covariate_init <- esvd_init$covariates[,case_control_variable,drop = F]
b_init <- esvd_init$b_mat[,case_control_variable,drop = F]

time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(esvd_init$x_mat,
                            esvd_init$y_mat,
                            mat,
                            family = "poisson",
                            nuisance_param_vec = NA,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = b_init,
                            covariates = covariate_init,
                            offset_vec = NULL,
                            offset_mat = offset_mat,
                            global_estimate = F,
                            l2pen = 0.1,
                            max_iter = 50,
                            run_cpp = F,
                            reparameterize = F,
                            reestimate_nuisance = F,
                            verbose = 1)
time_end2 <- Sys.time()

save(date_of_run, session_info, adams, covariate_df,
     esvd_init, time_start1, time_end1,
     esvd_res, time_start2, time_end2,
     file = "../../../../out/Writeup11/Writeup11_adams_ciliated_esvd_tmp.RData")

##################

print("Starting final fit, where library size coef can change")

covariates <- esvd_init$covariates
b_mat <- esvd_init$b_mat
b_mat[,case_control_variable] <- esvd_res$b_mat[,case_control_variable]

time_start3 <- Sys.time()
esvd_res_full <- eSVD2::opt_esvd(esvd_res$x_mat,
                                 esvd_res$y_mat,
                                 mat,
                                 family = "poisson",
                                 nuisance_param_vec = NA,
                                 library_size_vec = 1,
                                 method = "newton",
                                 b_init = b_mat,
                                 covariates = covariates,
                                 offset_vec = rep(0, nrow(mat)),
                                 offset_mat = NULL,
                                 global_estimate = F,
                                 l2pen = 0.1,
                                 max_iter = 50,
                                 run_cpp = T,
                                 reparameterize = F,
                                 reestimate_nuisance = F,
                                 verbose = 1)
time_end3 <- Sys.time()

save(date_of_run, session_info, adams, covariate_df,
     esvd_init, time_start1, time_end1,
     esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     file = "../../../../out/Writeup11/Writeup11_adams_ciliated_esvd_tmp.RData")

###########

print("Starting nuisance parameter estimation")
nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,case_control_variable,drop = F],
                       esvd_res_full$b_mat[,case_control_variable,drop = F])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res_full$covariates[,offset_var],
  esvd_res_full$b_mat[,offset_var]
))

time_start4 <- Sys.time()
nuisance_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  val <- tryCatch(eSVD2:::gamma_rate(x = mat[,j],
                                     mu = mean_mat_nolib[,j],
                                     s = library_mat[,j]),
                  error = function(c) 0)
  val
})
time_end4 <- Sys.time()

save(date_of_run, session_info, adams, covariate_df,
     esvd_init, time_start1, time_end1,
     esvd_res, time_start2, time_end2,
     esvd_res_full, time_start3, time_end3,
     nuisance_vec, time_start4, time_end4,
     file = "../../../../out/Writeup11/Writeup11_adams_ciliated_esvd_tmp.RData")

print("Finished")
