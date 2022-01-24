rm(list=ls())
load("../../../../out/writeup8g/writeup8g_sns_layer23_esvd.RData")
source("../writeup8g/fano_nuisance.R")
library(Rmpfr)

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

library_idx <- which(colnames(esvd_res_full$covariates) == "Log_UMI")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
mean_mat_nolib <- pmin(mean_mat_nolib, 1e4)

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})

nuisance_vec <- rep(NA, ncol(mat))
for(j in 1:ncol(mat)){
  print(j)
  if(j %% 100 == 0) {
    save.image("../../../../out/writeup8g/writeup8g_sns_layer23_esvd_postprocess_rmpfr_tmp.RData")
  }
  nuisance_vec[j] <- calculate_fano_parameter(y = mat[,j],
                                              mu = mean_mat_nolib[,j],
                                              sf = library_mat[,j])
}

save.image("../../../../out/writeup8g/writeup8g_sns_layer23_esvd_postprocess_rmpfr.RData")
