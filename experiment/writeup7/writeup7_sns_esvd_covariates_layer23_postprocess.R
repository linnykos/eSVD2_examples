rm(list=ls())

library(Seurat)

load("../../../../out/writeup7/writeup7_sns_esvd_covariates_layer23_2.RData")

png("../../../../out/fig/writeup7/sns_esvd_covariates_layer23_asd_hist.png",
    width = 1800, height = 1500, units = "px", res = 300)
hist(esvd_res$b_mat[,3], main = "Human brain (SNS, with covariates)\neSVD, ASD coef. histogram",
     col = "gray", xlab = "ASD coefficient", breaks = 50)
graphics.off()
