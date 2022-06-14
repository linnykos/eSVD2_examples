plot(eSVD_obj$fit_First$z_mat[,"covariate_1"])

nat_mat1 <- tcrossprod(eSVD_obj$fit_First$x_mat, eSVD_obj$fit_First$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates, eSVD_obj$fit_First$z_mat)
pred_mat <- exp(nat_mat1 + nat_mat2)
image(pred_mat)
plot(pred_mat[,1])
plot(pred_mat[,10])
plot(pred_mat[,100])
plot(pred_mat[,150])

esvd_res <- eSVD_obj$fit_First
covariates <- eSVD_obj$covariates
library_size_variable <- "Log_UMI"
# library_idx <- which(colnames(covariates) %in% c("Intercept", library_size_variable))
# library_idx <- which(colnames(covariates) %in% library_size_variable)
relevant_idx <- which(colnames(covariates) %in% c("Intercept", "case_control"))
nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,relevant_idx], esvd_res$z_mat[,relevant_idx])
pred_mat <- exp(nat_mat1 + nat_mat2)
image(pred_mat)

plot(pred_mat[,1])
plot(pred_mat[,10])
plot(pred_mat[,100])
plot(pred_mat[,150])
image(log(pred_mat))

plot(esvd_res$z_mat[,"case_control"])
plot(esvd_res$z_mat[,"gender"])
plot(esvd_res$z_mat[,"Intercept"])
plot(esvd_res$z_mat[,"Log_UMI"])
plot(eSVD_obj$covariates[,"Log_UMI"])

esvd_res <- eSVD_obj$fit_Init
plot(eSVD_obj$initial_Reg$z_mat1[,"covariate_1"])
plot(eSVD_obj$initial_Reg$log_pval)
plot(esvd_res$z_mat[,"covariate_1"])
plot(esvd_res$z_mat[,"covariate_2"])
plot(esvd_res$z_mat[,"covariate_3"])
plot(esvd_res$z_mat[,"covariate_4"])
plot(esvd_res$z_mat[,"Intercept"])
plot(esvd_res$z_mat[,"Log_UMI"])

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(eSVD_obj$covariates, esvd_res$z_mat)
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_idx <- c(1,6)
library_mat <- exp(tcrossprod(
  eSVD_obj$covariates[,library_idx], esvd_res$z_mat[,library_idx]
))
# library_mat <- pmin(library_mat, 500)

Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = eSVD_obj$fit_First$nuisance_vec, FUN = "*")
AplusAlpha <- as.matrix(eSVD_obj$dat + Alpha)
image(AplusAlpha)
image(log(AplusAlpha))
quantile(AplusAlpha)
plot(AplusAlpha[,1])

SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = eSVD_obj$fit_First$nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2


###################

hist(eSVD_obj$fit_First$nuisance_vec)
plot(eSVD_obj$fit_First$nuisance_vec, jitter(nuisance_vec), asp = T)

zz <- eSVD_obj$fit_First$posterior_mean_mat
plot(zz[,1])
plot(zz[,100])
plot(zz[,150])
plot(zz[,90])

image(eSVD_obj$fit_First$posterior_mean_mat)
image(eSVD_obj$fit_First$posterior_var_mat)
image(log(eSVD_obj$fit_First$posterior_mean_mat))
image(log(eSVD_obj$fit_First$posterior_var_mat))
zz <- eSVD_obj$fit_First$posterior_mean_mat
zz <- pmin(zz, 1)
image(zz)

eSVD_obj2 <- compute_posterior(input_obj = eSVD_obj,
                               bool_adjust_covariates = F)
image(eSVD_obj2$fit_First$posterior_mean_mat)
zz <- eSVD_obj2$fit_First$posterior_mean_mat
zz <- pmin(zz, 1)
image(zz)

######################
case_control_variable <- "covariate_1"
case_control_idx <- which(colnames(covariates) == case_control_variable)
esvd_res <- eSVD_obj2$fit_First
covariates <- eSVD_obj2$covariates
library_size_variable <- "Log_UMI"
library_idx <- which(colnames(covariates) %in% library_size_variable)

hist(covariates[,library_size_variable])
hist(esvd_res$z_mat[,library_size_variable])
hist(esvd_res$z_mat[,library_size_variable])

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,-library_idx],
                       esvd_res$z_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_tmp <- tcrossprod(
  covariates[,library_idx], esvd_res$z_mat[,library_idx]
)
image(library_tmp)
library_mat <- exp(library_tmp)
library_mat <- pmin(library_mat, 500)

plot(as.numeric(library_mat), rep(Matrix::rowSums(dat), times = p), asp = T)
plot(as.numeric(library_mat), rep(Matrix::rowSums(dat), times = p))

#########

library_tmp <- tcrossprod(
  covariates[,"Log_UMI"], esvd_res$z_mat[,"Log_UMI"]
)
library_mat <- exp(library_tmp)
library_mat <- pmin(library_mat, 500)
plot(as.numeric(library_mat), rep(Matrix::rowSums(dat), times = p))

##########################

mean_mat_nolib <- exp(nat_mat1 + nat_mat2)
nuisance_vec <- eSVD_obj$fit_First$nuisance_vec
Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- as.matrix(eSVD_obj$dat + Alpha)
image(AplusAlpha)
plot(AplusAlpha[,1])
plot(AplusAlpha[,50])
plot(AplusAlpha[,100])
plot(AplusAlpha[,150])
plot(AplusAlpha[,90])

##########################
##########################
##########################
##########################
##########################

k <- 150
plot(eSVD_obj$fit_First$posterior_mean_mat[,k],
     eSVD_obj$fit_First$posterior_var_mat[,k], asp = T,
     col = covariates[,"case_control"]+1,
     main = k)

##########################
##########################

input_obj <- eSVD_obj
esvd_res <- eSVD_obj$fit_First
nuisance_vec <- esvd_res$nuisance_vec
case_control_variable <- .get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- .get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])

covariate_individual = "individual"
individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
stopifnot(length(intersect(control_individuals, case_individuals)) == 0)

latest_Fit <- .get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
relevant_idx <- which(colnames(covariates) %in% "case_control")
nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(covariates[,relevant_idx], esvd_res$z_mat[,relevant_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  covariates[,-relevant_idx], esvd_res$z_mat[,-relevant_idx]
))
nuisance_vec <- pmax(nuisance_vec,
                     stats::quantile(nuisance_vec, probs = 0.01))
Alpha <- sweep(mean_mat_nolib, MARGIN = 2,
               STATS = nuisance_vec, FUN = "*")
AplusAlpha <- as.matrix(eSVD_obj$dat) + Alpha
SplusBeta <- sweep(library_mat, MARGIN = 2,
                   STATS = nuisance_vec, FUN = "+")
posterior_mean_mat <- AplusAlpha/SplusBeta
posterior_var_mat <- AplusAlpha/SplusBeta^2

teststat_vec <- compute_test_statistic.default(
  input_obj = posterior_mean_mat,
  posterior_var_mat = posterior_var_mat,
  case_individuals = case_individuals,
  control_individuals = control_individuals,
  covariate_individual = covariate_individual,
  metadata = metadata,
  verbose = 1
)

plot(teststat_vec,
     pch = 16,
     col = true_cc_status)

k <- 100
plot(posterior_mean_mat[,k], posterior_var_mat[,k], asp = T,
     col = covariates[,"case_control"]+1)
