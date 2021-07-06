rm(list=ls())
set.seed(5)
canon_vec <- 1:50
nat_mat <- sapply(log(canon_vec)+1, function(x){rep(x, 100)})
dat <- generate_data(nat_mat, family = "poisson")

nuisance_param_vec = NA
k = 2
family = "poisson"
library_size_vec = NA
covariates = matrix(1, nrow = nrow(dat), ncol = 1)
config = initialization_options()
verbose = 0

#####

stopifnot(is.character(family))
if(family != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))
if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
  nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
}

# estimate library sizes if asked
n <- nrow(dat); p <- ncol(dat)
family <- .string_to_distr_funcs(family)
library_size_vec <- .parse_library_size(dat, library_size_vec = library_size_vec)
rescaled_dat <- t(sapply(1:nrow(dat), function(i){
  dat[i,]/library_size_vec[i]
}))

# determine initial matrix taking into account to missing values and library size
# [note to self: this probably could be something a lot simpler]
rescaled_dat <- .matrix_completion(rescaled_dat, k = k)
init_res <- .determine_initial_matrix(rescaled_dat, k = k, family = family,
                                      nuisance_param_vec = nuisance_param_vec,
                                      max_val = config$max_val,
                                      tol = config$tol)
nat_mat <- init_res$nat_mat; domain <- init_res$domain

if(all(is.null(covariates)))
{
  b_mat <- NULL
} else {
  r <- ncol(covariates)
  # do regression
  tmp <- lapply(1:p, function(j){
    .regress_covariates(nat_mat[,j], covariates)
  })

  nat_mat <- sapply(1:p, function(j){tmp[[j]]$residual})
  b_mat <- sapply(1:p, function(j){tmp[[j]]$coef})
  baseline <- tcrossprod(covariates, b_mat)
}

# project inital matrix into space of low-rank matrices
nat_mat2 <- .initialize_nat_mat(nat_mat, k = k, baseline = baseline,
                               domain = domain, config = config)

Matrix::rankMatrix(nat_mat2)
quantile(nat_mat)
# res <- .factorize_matrix(nat_mat, k = k, equal_covariance = T)


###########3

svd_res <- .svd_truncated(nat_mat, K = k, symmetric = F, rescale  = F,
                          mean_vec = NULL, sd_vec = NULL, K_full_rank = F)
x_mat <- .mult_mat_vec(svd_res$u, sqrt(svd_res$d))
y_mat <- .mult_mat_vec(svd_res$v, sqrt(svd_res$d))
