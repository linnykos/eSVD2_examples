rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_10_esvd.RData")

dispersion_est_vec <- eSVD_obj$fit_Second$nuisance_vec
posterior_mat <- eSVD_obj$fit_Second$posterior_mean_mat
p <- ncol(posterior_mat)
denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                           eSVD_obj$fit_Second$y_mat) +
  tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
denoised_mat <- exp(denoised_mat)
sparsity_vec <- sapply(1:p, function(j){
  length(which(eSVD_obj$dat[,j] == 0))/nrow(eSVD_obj$dat)
})

cor_vec <- sapply(1:p, function(j){
  stats::cor(denoised_mat[,j], posterior_mat[,j])
})
quantile(cor_vec)

plot(jitter(1/dispersion_est_vec), cor_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.5))
plot(sparsity_vec, cor_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.5))

plot(res2$pvalue_vec[-c(1:10)])
plot(multtest_res$pvalue_vec[-c(1:10)])

##########################

p <- ncol(eSVD_obj$dat)
col_vec <- rep(2, p)
col_vec[eSVD_obj$covariates[,"CC"] == 1] <- 3
j <- target_null_genes[1]
mean_vec <- eSVD_obj$fit_Second$posterior_mean_mat[,j]
mean_vec <- pmin(mean_vec, quantile(mean_vec, probs = 0.95))
var_vec <- eSVD_obj$fit_Second$posterior_var_mat[,j]
var_vec <- pmin(var_vec, quantile(var_vec, probs = 0.95))
plot(jitter(mean_vec), jitter(sqrt(var_vec)),
     asp = F, pch = 16, col = col_vec)
plot(mean_vec,
     asp = F, pch = 16, col = col_vec)

plot(denoised_mat[,j], pch = 16, col = col_vec,
     ylim = quantile(denoised_mat[,j], probs = c(0,0.95)))

dat <- eSVD_obj$dat
plot(dat[,950])
plot(dat[,950], ylim = c(0,5))

#############################

bool_adjust_covariates = F
alpha_max = 200
bool_covariates_as_library = T
bool_stabilize_underdispersion = T
library_min = 0.1
pseudocount = 1

input_obj <- eSVD_obj$dat
input_obj <- pmax(input_obj, 1)
case_control_variable <- .get_object(eSVD_obj = eSVD_obj, which_fit = "param", what_obj = "init_case_control_variable")
library_size_variable <- .get_object(eSVD_obj = eSVD_obj, which_fit = "param", what_obj = "init_library_size_variable")
bool_library_includes_interept <- .get_object(eSVD_obj = eSVD_obj, which_fit = "param", what_obj = "nuisance_bool_library_includes_interept")

bool_covariates_as_library = T
bool_library_includes_interept = T

covariates <- eSVD_obj$covariates
case_control_idx <- which(colnames(covariates) == case_control_variable)
library_size_variables <- library_size_variable
if(bool_covariates_as_library) library_size_variables <- unique(c(library_size_variables,
                                                                  setdiff(colnames(covariates),
                                                                          c("Intercept", case_control_variable))))
if(bool_library_includes_interept) library_size_variables <-  unique(c("Intercept", library_size_variables))

library_idx <- which(colnames(covariates) %in% library_size_variables)
idx_vec <- c(case_control_idx, library_idx)

esvd_res <- eSVD_obj$fit_Second
nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,-library_idx], esvd_res$z_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  covariates[,library_idx], esvd_res$z_mat[,library_idx]
))
if(!is.null(library_min)) library_mat <- pmax(library_mat, library_min)

nuisance_vec <- eSVD_obj$fit_Second$nuisance_vec
nuisance_vec[901:1000] <- 0
nuisance_lower_quantile = 0.01
nuisance_vec <- pmax(nuisance_vec,
                     stats::quantile(nuisance_vec, probs = nuisance_lower_quantile))
bool_stabilize_underdispersion = T
if(bool_stabilize_underdispersion & mean(log10(nuisance_vec)) > 0) {
  nuisance_vec <- 10^(scale(log10(nuisance_vec), center = T, scale = F))
}

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
quantile(Alpha[,950])
plot(Alpha[,950], ylim = c(0,10))
plot(mean_mat_nolib[,950], ylim = c(0,1e4))
Alpha <- pmin(Alpha, 200)
AplusAlpha <- input_obj + Alpha
quantile(AplusAlpha[,950])

if(!is.null(alpha_max)) AplusAlpha <- pmin(AplusAlpha, alpha_max)

SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
quantile(SplusBeta[,950])
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2

j <- 950
quantile(AplusAlpha[,j])
quantile(SplusBeta[,j])

plot(posterior_mean_mat[,j], eSVD_obj$dat[,j])
# plot(posterior_mean_mat[,j], denoised_mat[,j])

cc_vec <- eSVD_obj[["case_control"]]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])
mean(dat[control_idx,950])
mean(dat[case_idx,950])
sd(dat[control_idx,950])
sd(dat[case_idx,950])

individual_vec <- eSVD_obj[["individual"]]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
stopifnot(length(intersect(control_individuals, case_individuals)) == 0)

p <- ncol(posterior_mean_mat)

tmp <- .determine_individual_indices(case_individuals = case_individuals,
                                     control_individuals = control_individuals,
                                     individual_vec = individual_vec)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                       n = nrow(posterior_mean_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

j <- 950
plot(avg_posterior_mean_mat[,j], sqrt(avg_posterior_var_mat[,j]),
     pch = 16, asp = T, col = c(rep(2, length(case_individuals)),
                                rep(3, length(control_individuals))))

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
teststat_vec[j]

