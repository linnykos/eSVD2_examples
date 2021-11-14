rm(list=ls())
library(Seurat)

load("../../../../out/writeup8d/writeup8d_sns_layer23_esvd.RData")
true_esvd <- esvd_res2
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("true_esvd")]
rm(list = ls_vec)

n <- nrow(true_esvd$x_mat)
p <- nrow(true_esvd$y_mat)
colnames(true_esvd$b_mat) <- colnames(true_esvd$covariates)
autism_idx <- which(colnames(true_esvd$covariates) == "diagnosis_ASD")
library_idx <- which(colnames(true_esvd$covariates) == "Log-UMI")

#################3
# fix the simulation data
true_esvd$b_mat[,autism_idx] <- 0
true_esvd$b_mat[,library_idx] <- 1

y_max <- quantile(abs(true_esvd$y_mat), prob = 0.99)
true_esvd$y_mat[which(abs(true_esvd$y_mat) >= y_max)] <- y_max*sign(true_esvd$y_mat[which(abs(true_esvd$y_mat) >= y_max)])

x_max <- quantile(abs(true_esvd$x_mat), prob = 0.99)
true_esvd$x_mat[which(abs(true_esvd$x_mat) >= x_max)] <- x_max*sign(true_esvd$x_mat[which(abs(true_esvd$x_mat) >= x_max)])

b_max <- quantile(abs(true_esvd$b_mat), prob = 0.95)
true_esvd$b_mat[which(abs(true_esvd$b_mat) >= b_max)] <- b_max*sign(true_esvd$b_mat[which(abs(true_esvd$b_mat) >= b_max)])

x_row_max <- matrixStats::rowMaxs(abs(true_esvd$x_mat))
idx <- which(x_row_max >= 2.5)
true_esvd$x_mat[idx,] <- true_esvd$x_mat[idx,]/2

y_row_max <- matrixStats::rowMaxs(abs(true_esvd$y_mat))
idx <- which(y_row_max >= 1)
true_esvd$y_mat[idx,] <- true_esvd$y_mat[idx,]/2

b_row_max <- matrixStats::rowMaxs(abs(true_esvd$b_mat))
idx <- which(b_row_max >= 3)
true_esvd$b_mat[idx,] <- true_esvd$b_mat[idx,]/5

true_esvd$b_mat[,c(3:6)] <- true_esvd$b_mat[,c(3:6)]/2
true_esvd$b_mat[,"age"] <- true_esvd$b_mat[,"age"]/10
true_esvd$b_mat[,"RNA.Integrity.Number"] <- true_esvd$b_mat[,"RNA.Integrity.Number"]/2
true_esvd$b_mat[,"post.mortem.hours"] <- true_esvd$b_mat[,"post.mortem.hours"]/10
true_esvd$b_mat[,"nFeature_RNA"] <- true_esvd$b_mat[,"nFeature_RNA"]/1000

################

nat_mat1 <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat)
nat_mat2 <- tcrossprod(true_esvd$covariates[,-library_idx], true_esvd$b_mat[,-library_idx])
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)
mean_mat[mean_mat >= 100] <- 100

set.seed(10)
autism_gene_idx <- sample(1:p, size = round(p/50))
multiplier_vec <- rnorm(length(autism_gene_idx), mean = 3, sd = 1)
for(j in 1:length(autism_gene_idx)){
  vec <- nat_mat[,autism_gene_idx[j]]
  vec_autism <- vec[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5)]
  vec_control <- vec[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5)]

  target <- multiplier_vec[j] * max(mean(exp(vec_control)), mean(exp(vec_autism)))
  current <- mean(exp(vec_autism))
  value <-  log(target/current)
  true_esvd$b_mat[autism_gene_idx[j],autism_idx] <- value
}

nat_mat1 <- tcrossprod(true_esvd$x_mat, true_esvd$y_mat)
nat_mat2 <- tcrossprod(true_esvd$covariates[,-library_idx], true_esvd$b_mat[,-library_idx])
nat_mat <- nat_mat1 + nat_mat2
nat_mat[nat_mat >= log(100)] <- log(100)
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx])
quantile(nat_mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx])

#################

# now simulate data
# first simulate the gamma
lambda_mat <- nat_mat
for(j in 1:ncol(lambda_mat)){
  set.seed(j)
  lambda_mat[,j] <- stats::rgamma(n,
                           shape = true_esvd$nuisance_param_vec[j],
                           scale = exp(nat_mat[,j])/true_esvd$nuisance_param_vec[j])
}

# scale up lambda according to the library size
s_vec <- p*exp(true_esvd$covariates[,library_idx])/matrixStats::rowSums2(lambda_mat)

mat <- lambda_mat
tol <- 1e-3
for(i in 1:nrow(mat)){
  set.seed(i)
  mat[i,] <- stats::rpois(p, lambda = s_vec[i]*lambda_mat[i,] + tol)
}
true_esvd$covariates[,library_idx] <- log(matrixStats::rowMeans2(mat))

tmp <- mat[which(true_esvd$covariates[,"diagnosis_ASD"] > 0.5), autism_gene_idx]
length(which(tmp == 0))/prod(dim(tmp))
quantile(tmp[tmp > 0])

tmp <- mat[which(true_esvd$covariates[,"diagnosis_ASD"] < 0.5), autism_gene_idx]
length(which(tmp == 0))/prod(dim(tmp))
quantile(tmp[tmp > 0])

######################
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("true_esvd", "mat", "autism_gene_idx",
                                "nat_mat", "lambda_mat" , "gamma_mat", "s_vec")]
rm(list = ls_vec)
save.image("../../../../out/writeup8d/writeup8d_sns_pseudoreal_data.RData")


true_mean <- eSVD2:::.mult_vec_mat(s_vec, lambda_mat)
max_value <- 30
png("../../../../out/fig/writeup8d/sns_pseudoreal_true_scatterplot.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                            mean_mat = true_mean,
                            size_vec = true_esvd$nuisance_param_vec,
                            quantile_shoulder = 0.5,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Pseudoreal: True",
                            included_col = rgb(0.5, 0.5, 0.5, 0.5),
                            excluded_col = rgb(0.5, 0, 0, 0.5),
                            include_percentage_in_main = T,
                            verbose = T)
graphics.off()


######################################

# now fit

# initialization
K <- 10
n <- nrow(mat)
p <- ncol(mat)
covariates <- true_esvd$covariates

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

