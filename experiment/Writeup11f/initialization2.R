initialize_esvd2 <- function(dat,
                             covariates,
                             subject_variables,
                             case_control_variable = NULL,
                             k = 30,
                             lambda = 0.01,
                             library_size_variable = "Log_UMI",
                             offset_variables = c("Intercept", "Log_UMI"),
                             verbose = 0){
  stopifnot(inherits(dat, c("dgCMatrix", "matrix")),
            nrow(dat) == nrow(covariates),
            is.matrix(covariates),
            k <= ncol(dat), k > 0, k %% 1 == 0,
            lambda <= 1e4, lambda >= 1e-4,
            library_size_variable %in% colnames(covariates),
            is.null(case_control_variable) || case_control_variable %in% colnames(covariates),
            all(subject_variables %in% colnames(covariates)),
            "Intercept" %in% colnames(covariates))
  stopifnot(all(is.null(offset_variables)) || all(offset_variables %in% colnames(covariates)))

  n <- nrow(dat); p <- ncol(dat)
  if(is.matrix(dat)) dat[is.na(dat)] <- 0
  param <- eSVD2:::.format_param_initialize(bool_intercept = NULL,
                                            case_control_variable = case_control_variable,
                                            k = k,
                                            lambda = lambda,
                                            library_size_variable = library_size_variable,
                                            offset_variables = offset_variables,
                                            subject_variables = subject_variables)

  if(verbose >= 1) print("Performing GLMs")
  z_mat <- .initialize_coefficient2(case_control_variable = case_control_variable,
                                    covariates = covariates,
                                    dat = dat,
                                    lambda = lambda,
                                    offset_variables = offset_variables,
                                    subject_variables = subject_variables,
                                    verbose = verbose)

  eSVD_obj <- list(dat = dat,
                   covariates = covariates,
                   param = param)
  class(eSVD_obj) <- "eSVD"

  if(verbose >= 1) print("Computing residuals")
  eSVD_obj[["fit_Init"]] <- eSVD2:::.initialize_residuals(
    covariates = covariates,
    dat = dat,
    k = k,
    z_mat = z_mat
  )

  eSVD_obj[["latest_Fit"]] <- "fit_Init"
  eSVD_obj
}


.initialize_coefficient2 <- function(case_control_variable,
                                     covariates,
                                     dat,
                                     lambda,
                                     offset_variables,
                                     subject_variables,
                                     verbose = 0){
  n <- nrow(dat); p <- ncol(dat)

  z_mat <- matrix(1, nrow = p, ncol = ncol(covariates))
  colnames(z_mat) <- colnames(covariates)
  rownames(z_mat) <- colnames(dat)

  stopifnot(!all(is.null(offset_variables)))
  covariates_tmp <- covariates[,which(!colnames(covariates) %in% offset_variables), drop=F]
  z_mat[,"Intercept"] <- log(Matrix::colMeans(dat))

  for(j in 1:p){
    if(verbose == 1 && p >= 10 && j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(paste0("Working on variable ", j , " of ", p, ": Intercept of ",
                                  round(z_mat[j,"Intercept"], 2)))

    offset_vec <- as.numeric(tcrossprod(covariates[,offset_variables], z_mat[j,offset_variables,drop = F]))

    glm_fit <- glmnet::glmnet(x = covariates_tmp,
                              y = as.numeric(dat[,j]),
                              family = "poisson",
                              offset = offset_vec,
                              alpha = 0,
                              standardize = F,
                              intercept = F,
                              lambda = exp(seq(log(1e4), log(lambda), length.out = 100)))

    z_mat[j, colnames(covariates_tmp)] <- glm_fit$beta[,ncol(glm_fit$beta)]
  }

  z_mat
}
