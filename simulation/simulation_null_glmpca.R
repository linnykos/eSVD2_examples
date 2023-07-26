rm(list=ls())
library(Seurat)
load("../eSVD2_examples/simulation/simulation_null.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- seurat_obj[["RNA"]]@counts
set.seed(10)
K <- 5
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "poi",
                             ctl = list(verbose = T),
                             minibatch = "stochastic",
                             X = covariate[,c("Log_UMI", "Sex", "Age")])

nat_mat1 <- tcrossprod(as.matrix(glmpca_res$factors),
                       as.matrix(glmpca_res$loadings))
nat_mat2 <- tcrossprod(as.matrix(glmpca_res$X[,"(Intercept)"]),
                       as.matrix(glmpca_res$coefX[,"(Intercept)"]))
nat_mat <- nat_mat1 + nat_mat2
nat_mat <- sweep(nat_mat, MARGIN = 1, STATS = glmpca_res$offsets, FUN = "+")

pred_mat <- exp(nat_mat)
case_idx <- which(covariate[,"CC"] == 1)
control_idx <- which(covariate[,"CC"] == 0)

p <- ncol(pred_mat)
pvalue_vec <- sapply(1:p, function(j){
  zz <- stats::wilcox.test(x = pred_mat[case_idx,j],
                           y = pred_mat[control_idx,j])
  zz$p.value
})

plot(sort(pvalue_vec[-c(1:10)]),
     seq(0,1,length.out = length(pvalue_vec[-c(1:10)])), asp = T)
lines(c(0,1), c(0,1), col = 2, lty = 2)

save(glmpca_res,
     pvalue_vec,
     date_of_run, session_info,
     file = "../eSVD2_examples/simulation/simulation_null_glmpca.RData")
