rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed.RData")

library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
covariate_dat <- sns@meta.data[,c("percent.mt", "individual", "region", "age", "sex",
                                  "RNA.Integrity.Number", "post.mortem.hours",
                                  "diagnosis", "Seqbatch")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"individual"] <- as.factor(covariate_df[,"individual"])
covariate_df[,"region"] <- as.factor(covariate_df[,"region"])
covariate_df[,"diagnosis"] <- factor(covariate_df[,"diagnosis"], levels = c("Control", "ASD"))
covariate_df[,"sex"] <- as.factor(covariate_df[,"sex"])
covariate_df[,"Seqbatch"] <- as.factor(covariate_df[,"Seqbatch"])
covariates <- format_covariates(dat = mat,
                                covariate_df = covariate_df,
                                mixed_effect_variables = c("individual", "Seqbatch"))

#####################

load("../../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]

case_control_variable = "diagnosis_ASD"
lambda = 0.01
offset_variables = c("Log_UMI")
offset_vec <- Matrix::rowSums(covariates[,offset_variables,drop=F])
p_val_thres = 0.05
verbose = 1
dat <- mat

covariates_nooffset <- covariates[,which(!colnames(covariates) %in% offset_variables)]

gene_idx <- which(colnames(dat) %in% de_gene_specific)
b_mat <- sapply(gene_idx, function(j){
  .lrt_coefficient(case_control_variable = case_control_variable,
                   covariates = covariates_nooffset,
                   include_intercept = T,
                   lambda = lambda,
                   offset_vec = offset_vec,
                   p_val_thres = p_val_thres,
                   vec = dat[,j],
                   verbose = 2,
                   verbose_gene_name = colnames(dat)[j])
})


gene_idx <- which(colnames(dat) %in% hk_genes)
b_mat2 <- sapply(gene_idx, function(j){
  .lrt_coefficient(case_control_variable = case_control_variable,
                   covariates = covariates_nooffset,
                   include_intercept = T,
                   lambda = lambda,
                   offset_vec = offset_vec,
                   p_val_thres = p_val_thres,
                   vec = dat[,j],
                   verbose = 2,
                   verbose_gene_name = colnames(dat)[j])
})

######################

vec <- dat[,"CDH2"]
# vec <- dat[,"LRRC47"]
tmp_cov <- covariates_nooffset[,c("Intercept", "percent.mt", "age", "RNA.Integrity.Number", "post.mortem.hours", "region_PFC", "sex_M", "diagnosis_ASD")]
df <- data.frame(cbind(vec, tmp_cov))
colnames(df)[1] <- "y"
glm1 <- stats::glm(y ~ . - 1, family = stats::poisson,
                   data = df,
                   offset = offset_vec)
tmp_cov2 <- tmp_cov[,which(colnames(tmp_cov) != "diagnosis_ASD")]
df2 <- df[,which(colnames(df) != "diagnosis_ASD")]
glm2 <- stats::glm(y ~ . - 1, family = stats::poisson,
                   data = df2,
                   offset = offset_vec)
anova(glm2, glm1, test = "Chi")

mean_vec1 <- exp(tmp_cov %*% stats::coef(glm1) + offset_vec)
log_vec <- log(vec/mean_vec1); log_vec[vec == 0] <- 0
deviance1 <- 2*sum(vec*log_vec - (vec - mean_vec1))

mean_vec2 <- exp(tmp_cov2 %*% stats::coef(glm2) + offset_vec)
log_vec <- log(vec/mean_vec2); log_vec[vec == 0] <- 0
deviance2 <- 2*sum(vec*log_vec - (vec - mean_vec2))
1-stats::pchisq(deviance2 - deviance1, df = 1)

