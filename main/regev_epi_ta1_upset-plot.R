rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regev_epi_ta1-noninflamed_sctransform.RData")
sctransform_noninflamed <- de_result
load("../../../out/main/regev_epi_ta1-inflamed_sctransform.RData")
sctransform_inflamed <- de_result

load("../../../out/main/regevEpi_ta1-inflamed_esvd.RData")
regevEpi_inflamed <- regevEpi
eSVD_obj_inflamed <- eSVD_obj
load("../../../out/main/regevEpi_ta1-noninflamed_esvd.RData")
regevEpi_noninflamed <- regevEpi
eSVD_obj_noninflamed <- eSVD_obj

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Non-Inflamed vs. He"))
sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Health"))
sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Non-In"))
noninf_de_genes <- sheet1[sheet1$ident == "TA 1","gene"]
inf_de_genes <- sheet2[sheet2$ident == "TA 1","gene"]
other_de_genes <- sheet3[sheet3$ident == "TA 1","gene"]
other_de_genes <- setdiff(other_de_genes, c(noninf_de_genes, inf_de_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
cycling_genes <- setdiff(cycling_genes, c(other_de_genes, noninf_de_genes, inf_de_genes))
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]

gene_names <- names(eSVD_obj$teststat_vec)
cycling_idx <- which(gene_names %in% cycling_genes)
inf_de_idx <- which(gene_names %in% inf_de_genes)
other_idx <- which(gene_names %in% other_de_genes)
noninf_de_idx <- which(gene_names %in% noninf_de_genes)
hk_idx <- which(gene_names %in% hk_genes)

other_idx <- setdiff(other_idx, c(inf_de_idx, noninf_de_idx))
cycling_idx <- setdiff(cycling_idx, c(inf_de_idx, noninf_de_idx, other_idx))
hk_idx <- setdiff(hk_idx, c(inf_de_idx, noninf_de_idx, other_idx, cycling_idx))

############################

# first for inflamed
eSVD_obj_inflamed$fit_Second$posterior_mean_mat <- NULL
eSVD_obj_inflamed$fit_Second$posterior_var_mat <- NULL
eSVD_obj_inflamed$teststat_vec <- NULL
eSVD_obj_inflamed <- eSVD2:::compute_posterior(input_obj = eSVD_obj_inflamed,
                                               bool_adjust_covariates = F,
                                               alpha_max = NULL,
                                               bool_covariates_as_library = T)
metadata <- regevEpi_inflamed@meta.data
metadata[,"Sample"] <- as.factor(metadata[,"Sample"])
eSVD_obj_inflamed <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj_inflamed,
                                                    covariate_individual = "Sample",
                                                    metadata = metadata,
                                                    verbose = 1)

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_inflamed,
                             metadata = regevEpi_inflamed@meta.data,
                             covariate_individual = "Sample")
teststat_vec <- eSVD_obj_inflamed$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_inflamed <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_inflamed <- locfdr::locfdr(gaussian_teststat_inflamed, plot = 0)
fdr_vec <- locfdr_res_inflamed$fdr
names(fdr_vec) <- names(gaussian_teststat_inflamed)
null_mean <- locfdr_res_inflamed$fp0["mlest", "delta"]
null_sd <- locfdr_res_inflamed$fp0["mlest", "sigma"]
logpvalue_vec_inflamed <- sapply(gaussian_teststat_inflamed, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec_inflamed <- -(logpvalue_vec_inflamed/log10(exp(1)) + log10(2))
idx_inflamed <- order(logpvalue_vec_inflamed, decreasing = T)[1:length(unique(c(inf_de_idx,noninf_de_genes)))]
esvd_inflamed_degenes <- names(teststat_vec)[idx_inflamed]

##

# next for non-inflamed
eSVD_obj_noninflamed$fit_Second$posterior_mean_mat <- NULL
eSVD_obj_noninflamed$fit_Second$posterior_var_mat <- NULL
eSVD_obj_noninflamed$teststat_vec <- NULL
eSVD_obj_noninflamed <- eSVD2:::compute_posterior(input_obj = eSVD_obj_noninflamed,
                                                  bool_adjust_covariates = F,
                                                  alpha_max = NULL,
                                                  bool_covariates_as_library = T)
metadata <- regevEpi@meta.data
metadata[,"Sample"] <- as.factor(metadata[,"Sample"])
eSVD_obj_noninflamed <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj_noninflamed,
                                                       covariate_individual = "Sample",
                                                       metadata = metadata,
                                                       verbose = 1)

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_noninflamed,
                             metadata = regevEpi_noninflamed@meta.data,
                             covariate_individual = "Sample")
teststat_vec <- eSVD_obj_noninflamed$teststat_vec
p <- length(teststat_vec)
gaussian_teststat_noninflamed <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res_noninflamed <- locfdr::locfdr(gaussian_teststat_noninflamed, plot = 0)
fdr_vec <- locfdr_res_noninflamed$fdr
names(fdr_vec) <- names(gaussian_teststat_noninflamed)
null_mean <- locfdr_res_noninflamed$fp0["mlest", "delta"]
null_sd <- locfdr_res_noninflamed$fp0["mlest", "sigma"]
logpvalue_vec_noninflamed <- sapply(gaussian_teststat_noninflamed, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec_noninflamed <- -(logpvalue_vec_noninflamed/log10(exp(1)) + log10(2))
idx_noninflamed <- order(logpvalue_vec_noninflamed, decreasing = T)[1:length(unique(c(inf_de_idx,noninf_de_genes)))]
esvd_noninflamed_degenes <- names(teststat_vec)[idx_noninflamed]

######################

sctransform_inflamed_degenes <- rownames(sctransform_inflamed)[order(sctransform_inflamed[,"p_val"], decreasing = F)[1:length(unique(c(inf_de_idx,noninf_de_genes)))]]
sctransform_noninflamed_degenes <- rownames(sctransform_noninflamed)[order(sctransform_noninflamed[,"p_val"], decreasing = F)[1:length(unique(c(inf_de_idx,noninf_de_genes)))]]

de_vec <- unique(c(noninf_de_genes, inf_de_genes))

length(noninf_de_genes)
length(inf_de_genes)
length(intersect(noninf_de_genes, inf_de_genes))

length(esvd_noninflamed_degenes)
length(intersect(esvd_noninflamed_degenes, noninf_de_genes))
length(intersect(esvd_noninflamed_degenes, de_vec))

length(esvd_inflamed_degenes)
length(intersect(esvd_noninflamed_degenes, inf_de_genes))
length(intersect(esvd_inflamed_degenes, de_vec))

length(intersect(esvd_inflamed_degenes, esvd_noninflamed_degenes))


length(sctransform_inflamed_degenes)
length(intersect(sctransform_inflamed_degenes, inf_de_genes))
length(intersect(sctransform_inflamed_degenes, de_vec))

length(sctransform_noninflamed_degenes)
length(intersect(sctransform_noninflamed_degenes, noninf_de_genes))
length(intersect(sctransform_noninflamed_degenes, de_vec))

length(intersect(sctransform_noninflamed_degenes, sctransform_inflamed_degenes))

