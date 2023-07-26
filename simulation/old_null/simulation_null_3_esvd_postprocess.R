rm(list=ls())
library(Seurat)
library(eSVD2)
load("../eSVD2_examples/simulation/simulation_null_3_esvd.RData")

1-length(eSVD_obj$dat@x)/prod(dim(eSVD_obj$dat)) # percent of 0's

# see how well we estimated the case-control
plot(z_mat[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"], asp = T)
plot(z_mat[,"Sex"], eSVD_obj$fit_Second$z_mat[,"Sex"], asp = T)
plot(z_mat[,"Age"], eSVD_obj$fit_Second$z_mat[,"Age"], asp = T)

# see how well we estimated the signal
mat1 <- tcrossprod(x_mat, y_mat) + tcrossprod(covariate[,"CC"], z_mat[,"CC"])
mat2 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
set.seed(10)
subset <- sample(1:prod(dim(mat1)), 1e5)
plot(as.numeric(mat1[subset]), as.numeric(mat2[subset]), asp = T,
     pch = 16, col = rgb(0.5,0.5,0.5,0.2))

# see how well we estimated the full signal
mat2 <- tcrossprod(eSVD_obj$fit_Second$x_mat, eSVD_obj$fit_Second$y_mat) + tcrossprod(eSVD_obj$covariates, eSVD_obj$fit_Second$z_mat)
set.seed(10)
subset <- sample(1:prod(dim(mat1)), 1e5)
plot(as.numeric(eSVD_obj$dat[subset]), as.numeric(exp(mat2[subset])), asp = T,
     pch = 16, col = rgb(0.5,0.5,0.5,0.2))

