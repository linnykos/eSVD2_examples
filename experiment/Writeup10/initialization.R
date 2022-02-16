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
  keep_idx <- which(colnames(covariates) %in% c("Intercept", "Log_UMI", "diagnosis_ASD"))
  other_idx <- which(colnames(covariates) %in% c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt",
                                                 "nFeature_RNA", "region_PFC", "sex_F", "Seqbatch_SB2", "Seqbatch_SB1"))
  for(j in other_idx){
    df_tmp <- data.frame(covariates[,j], covariates[,keep_idx])
    colnames(df_tmp)[1] <- "tmp"
    lm_fit <- stats::lm("tmp ~ . ", data = df_tmp)
    vec_tmp <- stats::residuals(lm_fit)
    covariates[,j] <- vec_tmp
  }
  
  cols_regress_out <- grep("individual", colnames(covariates))
  covariates_new <- covariates[,-cols_regress_out,drop = F]
  for(i in 1:length(cols_regress_out)){
    df_tmp <- data.frame(covariates[,cols_regress_out[i]], covariates_new)
    colnames(df_tmp)[1] <- "tmp"
    lm_fit <- stats::lm("tmp ~ . - 1", data = df_tmp)
    vec_tmp <- stats::residuals(lm_fit)
    if(sum(abs(vec_tmp)) < 1e-6) break()
    covariates_new <- cbind(covariates_new, vec_tmp)
    colnames(covariates_new)[ncol(covariates_new)] <- paste0("individual_",i)
  }
  covariates <- covariates_new
  
  family <- eSVD2:::.string_to_distr_funcs(family)
  if(family$name != "gaussian") stopifnot(all(dat[!is.na(dat)] >= 0))
  
  n <- nrow(dat); p <- ncol(dat)
  dat[is.na(dat)] <- 0
  
  ######
  # step: determine the coefficients via regression
  
  if(!all(is.null(covariates))){
    b_init <- sapply(1:ncol(covariates), function(j){
      if(stats::sd(covariates[,j]) == 0) {
        log(matrixStats::colMeans2(dat)+tol)
      } else {
        if(colnames(covariates)[j] %in% column_set_to_one){
          rep(1, ncol(dat))
        } else {
          rep(0, ncol(dat))
        }
      }
    })
    
    colnames(b_init) <- colnames(covariates)
    nat_offset_mat <- tcrossprod(covariates, b_init)
  } else {
    b_init <- NULL
    nat_offset_mat <- 0
  }
  
  nat_mat <- family$dat_to_nat(dat, gamma = rep(1, ncol(dat)))
  residual_mat <- nat_mat - nat_offset_mat
  residual_mat <- sweep(residual_mat, 1, offset_vec, "-")
  
  ## [[NOTE TO SELF: Can we replace this with a poisson GLM?]]
  remaining_covarites <- which(!colnames(covariates) %in% column_set_to_one)
  if(length(remaining_covarites) > 0){
    tmp <- eSVD2:::.regress_out_matrix(residual_mat,
                                       covariates[,remaining_covarites,drop = F],
                                       verbose = verbose)
    residual_mat <- tmp$residual_mat
    b_init[,remaining_covarites] <- tmp$b_mat
    residual_mat[is.na(residual_mat)] <- 0
    b_init[is.na(b_init)] <- 0
    nat_offset_mat <- tcrossprod(covariates, b_init)
  }
  rownames(b_init) <- colnames(dat)
  
  ######
  # step: compute SVD 
  
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
  
  structure(list(x_mat = x_init, y_mat = y_init, b_mat = b_init,
                 covariates = covariates,
                 nuisance_param_vec = rep(1, ncol(dat)),
                 offset_vec = offset_vec),
            class = "eSVD")
}
