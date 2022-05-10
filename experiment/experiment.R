rm(list=ls())
set.seed(123)
n <- 100
p <- 150
k <- 5
x_mat <- matrix(abs(rnorm(n * k))*.5, nrow = n, ncol = k)
y_mat <- matrix(abs(rnorm(p * k))*.5, nrow = p, ncol = k)
covariates <- cbind(c(rep(0, n/2), rep(1, n/2)),
                    matrix(abs(rnorm(n * 3, mean = 1, sd = 0.1)), nrow = n, ncol = 3))
colnames(covariates) <- paste0("covariate_", 1:4)
z_mat <- cbind(c(rep(0, p/2), rep(2, p/2)), rep(1,p), rep(1,p), rep(1,p))
colnames(z_mat) <-  colnames(covariates)
case_control_variable <- "covariate_1"
case_control_idx <- which(colnames(z_mat) == case_control_variable)
nat_mat_nolib <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates[,case_control_idx], z_mat[,case_control_idx])
library_mat <- exp(tcrossprod(covariates[,-case_control_idx], z_mat[,-case_control_idx]))
nuisance_vec <- rep(c(5, 1, 1/5), times = 50)

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
dat <- Matrix::Matrix(dat, sparse = T)
rownames(dat) <- paste0("c", 1:n)
colnames(dat) <- paste0("g", 1:p)
metadata <- data.frame(individual = factor(rep(1:4, each = n/4)))
rownames(metadata) <- rownames(dat)

posterior_mean_mat <- gamma_mat
posterior_var_mat <- sweep(gamma_mat, MARGIN = 2, STATS = nuisance_vec, FUN = "/")

case_individuals = c("1", "2")
control_individuals = c("3", "4")
covariate_individual = "individual"
metadata = metadata

#####

p <- ncol(posterior_mean_mat)

if(verbose >= 1) print("Computing individual-level statistics")
tmp <- .determine_individual_indices(case_individuals = case_individuals,
                                     control_individuals = control_individuals,
                                     covariate_individual = covariate_individual,
                                     metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                       n = nrow(posterior_mean_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

case_row_idx <- 1:length(case_individuals)
control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = F])
case_gaussian_var <- .compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = F],
  avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = F]
)
control_gaussian_var <- .compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = F],
  avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = F]
)

n1 <- length(case_individuals)
n2 <- length(control_individuals)
teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
  (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
names(teststat_vec) <- colnames(posterior_mean_mat)

teststat_vec
