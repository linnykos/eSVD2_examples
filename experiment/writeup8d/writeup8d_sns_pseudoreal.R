rm(list=ls())
library(Seurat)

load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd.RData")
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("esvd_res2")]
rm(list = ls_vec)

n <- nrow(esvd_res2$x_mat)
p <- nrow(esvd_res2$y_mat)
autism_idx <- which(colnames(esvd_res2$covariates) == "diagnosis_ASD")

#################3
# fix the simulation data
esvd_res2$b_mat[,autism_idx] <- 0
esvd_res2$covariates[,"age"] <- esvd_res2$covariates[,"age"]/10
esvd_res2$covariates[,"RNA.Integrity.Number"] <- esvd_res2$covariates[,"RNA.Integrity.Number"]/2
esvd_res2$covariates[,"post.mortem.hours"] <- esvd_res2$covariates[,"post.mortem.hours"]/10
esvd_res2$covariates[,"nFeature_RNA"] <- esvd_res2$covariates[,"nFeature_RNA"]/1000

y_max <- quantile(abs(esvd_res2$y_mat), prob = 0.99)
esvd_res2$y_mat[which(abs(esvd_res2$y_mat) >= y_max)] <- y_max*sign(esvd_res2$y_mat[which(abs(esvd_res2$y_mat) >= y_max)])

x_max <- quantile(abs(esvd_res2$x_mat), prob = 0.99)
esvd_res2$x_mat[which(abs(esvd_res2$x_mat) >= x_max)] <- x_max*sign(esvd_res2$x_mat[which(abs(esvd_res2$x_mat) >= x_max)])

b_max <- quantile(abs(esvd_res2$b_mat), prob = 0.95)
esvd_res2$b_mat[which(abs(esvd_res2$b_mat) >= b_max)] <- b_max*sign(esvd_res2$b_mat[which(abs(esvd_res2$b_mat) >= b_max)])

x_row_max <- matrixStats::rowMaxs(abs(esvd_res2$x_mat))
idx <- which(x_row_max >= 2.5)
esvd_res2$x_mat[idx,] <- esvd_res2$x_mat[idx,]/2

y_row_max <- matrixStats::rowMaxs(abs(esvd_res2$y_mat))
idx <- which(y_row_max >= 1)
esvd_res2$y_mat[idx,] <- esvd_res2$y_mat[idx,]/2

b_row_max <- matrixStats::rowMaxs(abs(esvd_res2$b_mat))
idx <- which(b_row_max >= 3)
esvd_res2$b_mat[idx,] <- esvd_res2$b_mat[idx,]/5

esvd_res2$b_mat[,c(3:6)] <- esvd_res2$b_mat[,c(3:6)]/2

################

nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(esvd_res2$covariates, esvd_res2$b_mat)
nat_mat[nat_mat >= 5] <- 5

set.seed(10)
autism_gene_idx <- sample(1:p, size = round(p/100))
multiplier <- 5
for(j in autism_gene_idx){
  vec <- nat_mat[,j]
  vec_autism <- vec[which(esvd_res2$covariates[,"diagnosis_ASD"] > 0.5)]
  vec_control <- vec[which(esvd_res2$covariates[,"diagnosis_ASD"] < 0.5)]

  target <- multiplier * max(mean(exp(vec_control)), mean(exp(vec_autism)))
  current <- mean(exp(vec_autism))
  value <-  log(target/current)
  esvd_res2$b_mat[j,autism_idx] <- value
}

true_esvd_model <- esvd_res2
nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(esvd_res2$covariates, esvd_res2$b_mat)
nat_mat[nat_mat >= 5] <- 5
quantile(nat_mat[which(esvd_res2$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx])
quantile(nat_mat[which(esvd_res2$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx])

#################

# now simulate data
mat <- nat_mat
for(j in 1:ncol(mat)){
  set.seed(j)
  mat[,j] <- stats::rnbinom(n,
                            size = esvd_res2$nuisance_param_vec[j],
                            mu = exp(nat_mat[,j]))
}
mat[mat >= 200] <- 200
tmp <- mat[which(esvd_res2$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx]
length(which(tmp == 0))/prod(dim(tmp))
quantile(tmp[tmp > 0])

tmp <- mat[which(esvd_res2$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx]
length(which(tmp == 0))/prod(dim(tmp))
quantile(tmp[tmp > 0])

######################
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("true_esvd_model", "mat", "autism_gene_idx")]
rm(list = ls_vec)

# now fit

# initialization
K <- 10
n <- nrow(mat)
p <- ncol(mat)
covariates <- esvd_res2$covariates
esvd_res_truth <- esvd_res2

time_start1 <- Sys.time()
init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "neg_binom2",
                                   covariates = covariates,
                                   column_set_to_one = "Log-UMI",
                                   offset_vec = rep(0, nrow(mat)),
                                   verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

quantile(init_res$b_mat[autism_gene_idx,autism_idx])
quantile(init_res$b_mat[-autism_gene_idx,autism_idx])

###################3

print("Estimating NB via eSVD, round 1")
time_start2 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat,
                            family = "neg_binom2",
                            nuisance_param_vec = init_res$nuisance_param_vec,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = T,
                            global_estimate = T,
                            reparameterize = T,
                            bool_run_cpp = T,
                            max_iter = 50,
                            l2pen = 0.1,
                            verbose = 1)
time_end2 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

quantile(esvd_res$b_mat[autism_gene_idx,autism_idx])
quantile(esvd_res$b_mat[-autism_gene_idx,autism_idx])

print("Estimating NB via eSVD, round 2")
time_start3 <- Sys.time()
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(esvd_res$x_mat,
                             esvd_res$y_mat, mat,
                             family = "neg_binom2",
                             nuisance_param_vec = esvd_res$nuisance_param_vec,
                             library_size_vec = 1,
                             method = "newton",
                             b_init = esvd_res$b_mat,
                             covariates = esvd_res$covariates,
                             offset_vec = esvd_res$offset_vec,
                             reestimate_nuisance = T,
                             global_estimate = F,
                             reparameterize = T,
                             bool_run_cpp = T,
                             max_iter = 50,
                             tol = 1e-8,
                             l2pen = 0.1,
                             verbose = 1)
time_end3 <- Sys.time()
save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

quantile(esvd_res2$b_mat[autism_gene_idx,autism_idx])
quantile(esvd_res2$b_mat[-autism_gene_idx,autism_idx])




