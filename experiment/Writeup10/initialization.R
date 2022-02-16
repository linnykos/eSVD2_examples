initialize_esvd2 <- function(dat,
                             k,
                             family = "poisson",
                             covariates = NULL,
                             offset_vec = rep(0, nrow(dat)),
                             column_set_to_one = NULL,
                             tol = 1e-3,
                             verbose = 0){
  stopifnot(is.character(family),
            family == "poisson",
            all(is.null(covariates)) || is.matrix(covariates),
            length(offset_vec) == nrow(dat),
            k <= ncol(dat), k > 0, k %% 1 == 0)
  
  # step: compute the proper covariate matrix
  # regress all variables against intercept + Log_UMI + diagnosis_ASD
  covariates <- covariates[,-which(colnames(covariates) == "Intercept")]
  offset_vec <- covariates[,which(colnames(covariates) == "Log_UMI")]
  covariates <- covariates[,-which(colnames(covariates) == "Log_UMI")]
  covariates[,"nFeature_RNA"] <- scale(covariates[,"nFeature_RNA"], 
                                       center = F,
                                       scale = T)
  # move all the individual covariates to the end
  col_idx1 <- grep("individual_", colnames(covariates))
  col_idx2 <- c(grep("Capbatch_", colnames(covariates)),
               grep("Seqbatch_", colnames(covariates)))
  covariates <- cbind(covariates[,-c(col_idx1, col_idx2)], 
                      covariates[,col_idx2],
                      covariates[,col_idx1])
  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0
  
  # lm_res <- stats::lm(percent.mt ~ ., data = data.frame(covariates))
  # na_idx <- which(is.na(stats::coef(lm_res)))
  # if(length(na_idx) > 0){
  #   if(verbose >= 1) {
  #     print(paste0("Removed ", length(na_idx), " covariates due to singularity"))
  #     print(paste0("Removed covariates:"))
  #     print(colnames(covariates)[na_idx])
  #   }
  #   covariates <- covariates[,-na_idx]
  # }
  
  ######
  # step: determine the coefficients via regression
  coef_mat <- t(sapply(1:p, function(j){
    if(verbose == 1 && p > 10 && j %% floor(p/10) == 0) cat('*')
    
    # df <- as.data.frame(cbind(y = dat[,j], covariates))
    # colnames(df)[1] <- "tmp"
    # glm_fit <- stats::glm(tmp ~ . - 1, 
    #                       offset = offset_vec,
    #                       data = df, 
    #                       family = stats::poisson)
    glm_fit <- glmnet::glmnet(x = covariates,
                              y = dat[,j],
                              family = "poisson",
                              offset = offset_vec,
                              alpha = 0,
                              standardize = F,
                              intercept = F,
                              lambda = exp(seq(log(10000), log(0.01), length.out = 100)))
    if(verbose == 2) print(paste0("Iteration ", j, ": Deviance of ", 
                                  round(glm_fit$dev.ratio[length(glm_fit$dev.ratio)], 2)))
    
    ## [[note to self: can be improved using the stat helper functions]]
    c(glm_fit$beta[,ncol(glm_fit$beta)], glm_fit$dev.ratio[length(glm_fit$dev.ratio)])
  }))
  deviance_vec <- coef_mat[,ncol(coef_mat)]
  coef_mat <- coef_mat[,-ncol(coef_mat)]
  # coef_mat <- cbind(0, coef_mat)
  colnames(coef_mat) <- colnames(covariates)
  rownames(coef_mat) <- colnames(dat)
  
  ######
  # step: compute SVD 
  dat_transform <- log1p(dat)
  nat_mat <- tcrossprod(covariates, coef_mat)
  residual_mat <- dat_transform - nat_mat
  
  svd_res <- eSVD2:::.svd_truncated(residual_mat,
                                    K = k,
                                    symmetric = F,
                                    rescale = F,
                                    mean_vec = NULL,
                                    sd_vec = NULL,
                                    K_full_rank = F)
  x_init <- eSVD2:::.mult_mat_vec(svd_res$u, sqrt(svd_res$d))
  y_init <- eSVD2:::.mult_mat_vec(svd_res$v, sqrt(svd_res$d))
  
  rownames(x_init) <- rownames(dat)
  rownames(y_init) <- colnames(dat)
  
  structure(list(x_mat = x_init, y_mat = y_init, b_mat = coef_mat,
                 covariates = covariates,
                 nuisance_param_vec = rep(1, ncol(dat)),
                 offset_vec = offset_vec,
                 deviance_vec = deviance_vec),
            class = "eSVD")
}
