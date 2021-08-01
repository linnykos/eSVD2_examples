rm(list=ls())

dat <- mat
k = K
family = "neg_binom"
nuisance_param_vec = nuisance_vec
library_size_vec = 1
covariates = covariates
config = eSVD2::initialization_options()
verbose = 1

##########

stopifnot(is.character(family))
if(family != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))
if(all(!is.na(nuisance_param_vec)) & length(nuisance_param_vec) == 1) {
  nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
}

n <- nrow(dat); p <- ncol(dat)
family <- eSVD2:::.string_to_distr_funcs(family)
library_size_vec <- eSVD2:::.parse_library_size(dat, library_size_vec = library_size_vec)
rescaled_dat <- t(sapply(1:nrow(dat), function(i){
  dat[i,]/library_size_vec[i]
}))

rescaled_dat <- eSVD2:::.matrix_completion(rescaled_dat, k = k)
init_res <- eSVD2:::.determine_initial_matrix(rescaled_dat, k = k, family = family,
                                      nuisance_param_vec = nuisance_param_vec,
                                      max_val = config$max_val,
                                      tol = config$tol)
nat_mat <- init_res$nat_mat; domain <- init_res$domain

tmp <- lapply(1:p, function(j){
  eSVD2:::.regress_covariates(nat_mat[,j], covariates)
})

r <- ncol(covariates)
nat_mat <- sapply(1:p, function(j){tmp[[j]]$residual})
b_mat <- do.call(rbind, (lapply(1:p, function(j){tmp[[j]]$coef})))
baseline <- tcrossprod(covariates, b_mat)
k2 <- k + r #[[note to self: check that this is correct]]

# project inital matrix into space of low-rank matrices
nat_mat <- eSVD2:::.initialize_nat_mat(nat_mat, k = k2, baseline = baseline,
                               domain = domain, config = config)

res <-  eSVD2:::.factorize_matrix(nat_mat, k = k, equal_covariance = T)
res <-  eSVD2:::.fix_rank_defficiency(res$x_mat, res$y_mat, domain = domain)

#####################

zz <- tcrossprod(res$x_mat, res$y_mat) + tcrossprod(covariates, b_mat)
idx <- which(zz >= 0, arr.ind = T)

