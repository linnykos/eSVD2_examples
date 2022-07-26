rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta2-inflamed_esvd.RData")
regevEpi_inflamed <- regevEpi
eSVD_obj_inflamed <- eSVD_obj
load("../../../out/main/regevEpi_ta2-noninflamed_esvd.RData")
regevEpi_noninflamed <- regevEpi
eSVD_obj_noninflamed <- eSVD_obj

# next for non-inflamed
df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj_noninflamed,
                             metadata = regevEpi_noninflamed@meta.data,
                             covariate_individual = "Sample")
teststat_vec <- eSVD_obj_noninflamed$teststat_vec
p <- length(teststat_vec)
gaussian_teststat <- sapply(1:p, function(j){
  qnorm(pt(teststat_vec[j], df = df_vec[j]))
})

locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
fdr_vec <- locfdr_res$fdr
names(fdr_vec) <- names(gaussian_teststat)
null_mean <- locfdr_res$fp0["mlest", "delta"]
null_sd <- locfdr_res$fp0["mlest", "sigma"]
pvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    stats::pnorm(x, mean = null_mean, sd = null_sd)*2
  } else {
    (1-stats::pnorm(x, mean = null_mean, sd = null_sd))*2
  }
})
logpvalue_vec_noninflamed <- -log10(pvalue_vec)
idx_noninflamed <- which(fdr_vec < 0.1)

############################




