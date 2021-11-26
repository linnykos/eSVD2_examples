rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson4.RData")
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]

#########

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

mat[mat == 0.25] <- 0
nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  # MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
})
quantile(nuisance_param_vec)
nuisance_param_vec <- pmin(nuisance_param_vec, 1e5)

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
library_mat <- pmin(library_mat, 50000)
AplusR <- sweep(mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + 1
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)
tmp <- posterior_mean_mat/mean_mat
quantile(tmp)

##############

gene_idx <- which(colnames(mat) == "SAT2")
# gene_idx <- which(colnames(mat) == "DEXI")
par(mfrow = c(1,2))
plot(mean_mat[,gene_idx],
     jitter(mat[,gene_idx]), asp = T,
     col = rgb(0,0,0,0.05))
plot(posterior_mean_mat[,gene_idx],
     jitter(mat[,gene_idx]), asp = T,
     col = rgb(0,0,0,0.05))
quantile(posterior_mean_mat[,gene_idx]/mean_mat[,gene_idx])

tmp <- cbind(mat[,gene_idx], mean_mat[,gene_idx])
eSVD2:::.compute_principal_angle(tmp)

###############

png("../../out/fig/writeup8e/tmp.png",
    height = 2000, width = 1500, units = "px", res = 300)
image(t(esvd_res_full$x_mat))
graphics.off()

png("../../out/fig/writeup8e/tmp2.png",
    height = 2000, width = 1500, units = "px", res = 300)
image(t(esvd_res_full$y_mat))
graphics.off()
