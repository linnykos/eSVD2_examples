initialize_esvd2 <- function(dat,
                             k,
                             covariates = NULL,
                             batches_names = c("Capbatch_", "Seqbatch_"),
                             bool_intercept = T,
                             lib_coef = NULL,
                             rescale_vars = "nFeature_RNA",
                             subject_names = "individual_",
                             tol = 1e-3,
                             verbose = 0){
  stopifnot(all(is.null(covariates)) || is.matrix(covariates),
            k <= ncol(dat), k > 0, k %% 1 == 0)
  stopifnot(colnames(covariates)[1] == "Intercept")

  # step: compute the proper covariate matrix
  if(is.null(lib_coef)){
    if("Log_UMI" %in% rescale_vars){
      n <- nrow(dat)
      lib_coef <- sqrt(sum(covariates[,"Log_UMI"]^2)/(n-1))
    } else {
      lib_coef <- 1
    }
  }
  if(length(rescale_vars) > 0){
    for(var_name in rescale_vars){
      covariates[,var_name] <- scale(covariates[,var_name],
                                     center = F,
                                     scale = T)
    }
  }
  offset_vec <- lib_coef*covariates[,which(colnames(covariates) == "Log_UMI")]
  covariates <- covariates[,-which(colnames(covariates) == "Log_UMI")]
  # move all the individual covariates to the end
  col_idx1 <- grep(subject_names, colnames(covariates))
  if(!all(is.null(batches_names))){
    col_idx2 <- unlist(lapply(batches_names, function(batch){
      grep(batch, colnames(covariates))
    }))
    covariates <- cbind(covariates[,-c(col_idx1, col_idx2)],
                        covariates[,col_idx2],
                        covariates[,col_idx1])
  } else {
    covariates <- cbind(covariates[,-c(col_idx1)],
                        covariates[,col_idx1])
  }

  covariates_tmp <- covariates[,-which(colnames(covariates) == "Intercept")]
  if(!bool_intercept) covariates <- covariates_tmp

  n <- nrow(dat); p <- ncol(dat)
  dat <- as.matrix(dat)

  ######
  # step: determine the coefficients via regression
  coef_mat <- t(sapply(1:p, function(j){
    if(verbose == 2) print(j)
    if(verbose == 1 && p > 10 && j %% floor(p/10) == 0) cat('*')

    glm_fit <- glmnet::glmnet(x = covariates_tmp,
                              y = dat[,j],
                              family = "poisson",
                              offset = offset_vec,
                              alpha = 0,
                              standardize = F,
                              intercept = bool_intercept,
                              lambda = exp(seq(log(10000), log(0.01), length.out = 100)))

    ## [[note to self: can be improved using the stat helper functions]]
    if(bool_intercept){
      c(glm_fit$a0[length(glm_fit$a0)], glm_fit$beta[,ncol(glm_fit$beta)])
    } else {
      glm_fit$beta[,ncol(glm_fit$beta)]
    }

  }))

  if(bool_intercept){
    colnames(coef_mat) <- c("Intercept", colnames(covariates_tmp))
  } else {
    colnames(coef_mat) <- colnames(covariates_tmp)
  }
  rownames(coef_mat) <- colnames(dat)

  ######
  # step: compute SVD
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(cbind(covariates, offset_vec), cbind(coef_mat, 1))
  residual_mat <- dat_transform - nat_mat

  svd_res <- eSVD2:::.svd_safe(mat = residual_mat,
                               check_stability = T, # boolean
                               K = k, # positive integer
                               mean_vec = NULL, # boolean, NULL or vector
                               rescale = F, # boolean
                               scale_max = NULL, # NULL or positive integer
                               sd_vec = NULL)
  x_init <- eSVD2:::.mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_init <- eSVD2:::.mult_mat_vec(svd_res$v, sqrt(svd_res$d))

  rownames(x_init) <- rownames(dat)
  rownames(y_init) <- colnames(dat)

  structure(list(x_mat = x_init, y_mat = y_init, z_mat = coef_mat,
                 covariates = covariates,
                 offset_vec = offset_vec,
                 offset_coef = lib_coef),
            class = "eSVD")
}
