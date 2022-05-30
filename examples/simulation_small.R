rm(list=ls())
library(eSVD2)

set.seed(123)
num_indiv <- 20
n_per_indiv <- 100
p <- 150
k <- 5
n <- num_indiv * n_per_indiv
x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
covariates <- cbind(
  c(rep(0, n/2), rep(1, n/2)),
  c(rep(0, n/4), rep(1, n/4), rep(0, n/4), rep(1, n/4)),
  sapply(1:n, function(i){rnorm(1, mean = i/n*3, sd = 0.5)})
)
colnames(covariates) <- c("case_control", "gender", "Log_UMI")
for(i in 1:num_indiv){
  tmp <- rep(0, n)
  tmp[((i-1)*n_per_indiv+1):(i*n_per_indiv)] <- 1
  covariates <- cbind(covariates, tmp)
  colnames(covariates)[ncol(covariates)] <- paste0("individual_", i)
}
z_mat <- cbind(
  c(rep(0, .9*p), rep(.5, 2/3*(.1*p)), rep(1, 1/3*(.1*p))),
  rnorm(p),
  rep(1,p)
)
for(i in 1:num_indiv){
  tmp <- rnorm(p, mean = 0, sd = 1)
  z_mat <- cbind(z_mat, tmp)
}
colnames(z_mat) <-  colnames(covariates)
case_control_variable <- colnames(covariates)[1]
case_control_idx <- which(colnames(z_mat) == case_control_variable)
intercept <- .5
nat_mat_nolib <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx]) + intercept
library_mat <- exp(tcrossprod(covariates[,"Log_UMI"], z_mat[,"Log_UMI"]))
nuisance_vec <- rep(c(5, 1, 1/5), times = p/3)

# Simulate data
gamma_mat <- matrix(NA, nrow = n, ncol = p)
dat <- matrix(NA, nrow = n, ncol = p)
for(i in 1:n){
  for(j in 1:p){
    gamma_mat[i,j] <- stats::rgamma(n = 1,
                                    shape = nuisance_vec[j]*exp(nat_mat_nolib[i,j]),
                                    rate = nuisance_vec[j])
    dat[i,j] <- stats::rpois(n = 1, lambda = library_mat[i,j] * gamma_mat[i,j])
  }
}
dat <- pmin(dat, 200)
covariates[,"Log_UMI"] <- scale(log(Matrix::rowSums(dat)))
quantile(dat)
image(dat)
dat <- Matrix::Matrix(dat, sparse = T)
rownames(dat) <- paste0("c", 1:n)
colnames(dat) <- paste0("g", 1:p)
metadata <- data.frame(individual = factor(rep(1:num_indiv, each = n_per_indiv)))
rownames(metadata) <- rownames(dat)

# fit eSVD
eSVD_obj <- eSVD2::initialize_esvd(dat = dat,
                                   covariates = covariates,
                                   case_control_variable = case_control_variable,
                                   k = 5,
                                   lambda = 0.1,
                                   mixed_effect_variables = colnames(covariates)[grep("individual", colnames(covariates))],
                                   offset_variables = NULL,
                                   verbose = 1)
eSVD_obj <- eSVD2::apply_initial_threshold(eSVD_obj = eSVD_obj,
                                           pval_thres = 0.1)
eSVD_obj <- eSVD2::opt_esvd(input_obj = eSVD_obj,
                            max_iter = 50,
                            verbose = 1)
eSVD_obj <- eSVD2::estimate_nuisance(input_obj = eSVD_obj,
                                     verbose = 1)
eSVD_obj <- eSVD2::compute_posterior(input_obj = eSVD_obj,
                                     library_size_variable = "Log_UMI")
eSVD_obj <- eSVD2::compute_test_statistic(input_obj = eSVD_obj,
                                          covariate_individual = "individual",
                                          metadata = metadata)

plot(eSVD_obj$teststat_vec,
     pch = 16,
     col = ifelse(z_mat[,"case_control"] > 1e-6, 2, 1))
plot(eSVD_obj$fit_First$z_mat[,"case_control"],
     pch = 16,
     col = ifelse(z_mat[,"case_control"] > 1e-6, 2, 1))

###############

dat2 <- as.matrix(dat)
wilcoxon_vec <- sapply(1:p, function(j){
  stats::wilcox.test(
    x = dat2[1:(n/2),j],
    y = dat2[((n/2)+1):n,j]
  )$p.value
})
plot(wilcoxon_vec)

