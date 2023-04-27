simulation_generate_data <- function(seed){
  set.seed(seed)
  n_each <- 100; s <- 20; p <- 1000; k <- 2
  n <- n_each*s

  # form latent variables
  x_mat <- cbind(runif(n, min = -0.1, max = 1),
                 runif(n, min = -0.1, max = 1),
                 runif(n, min = -0.1, max = 1))
  # have a block structure for the y's
  y_block_assignment <- rep(c(1:3), times = ceiling(p/3))[1:p]
  y_centers <- matrix(c(1.5, 0.1, 0.1,
                        0.1, 1.5, 0.1,
                        0.1, 0.1, 1.5), ncol = 3, nrow = 3, byrow = T)
  y_mat <- y_centers[y_block_assignment,] + matrix(rnorm(p*3, mean = 0, sd = 0.1), ncol = 3, nrow = p)
  gene_plot_idx <- c(which(y_block_assignment == 1), which(y_block_assignment == 2), which(y_block_assignment == 3))

  # form covariates
  # first form the table
  df <- cbind(1,
              0,
              rep(c(0,1), each = s/2),
              rep(c(0,1), times = s/2),
              scale(round(rnorm(s, mean = 30, sd = 5)), center = F, scale = T),
              1:s)
  colnames(df) <- c("Intercept", "Log_UMI", "CC", "Sex", "Age", "Individual")
  # expand to covariate matrix
  covariate <- do.call(rbind, lapply(1:nrow(df), function(i){
    matrix(rep(df[i,], each = n_each), nrow = n_each, ncol = ncol(df))
  }))
  colnames(covariate) <- colnames(df)[1:ncol(df)]
  z_mat <- cbind(rep(0, p), # intercept
                 rep(0.1, p), # library
                 c(rep(1,10),rep(0, p-10)), # cc
                 rnorm(p, mean = 0, sd = 0.2), # sex
                 rnorm(p, mean = 0, sd = 0.5)) # age
  colnames(z_mat) <- colnames(df)[1:(ncol(df)-1)]
  # form nuisance
  dispersion_vec <- sample(rep(c(0.5, 1, 10), each = ceiling(p/3))[1:p])
  dispersion_vec[1:10] <- 10

  # generate data
  nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariate[,"CC"], z_mat[,"CC"])

  # for each gene, shrink all the cells in an individual to its mean
  cell_individual_list <- lapply(1:s, function(i){which(covariate[,"Individual"] == i)})
  shrink_percentage <- 0.6 # higher means we shrink more
  for(j in 1:p){
    for(idx in cell_individual_list){
      mean_val <- mean(nat_mat[idx,j])
      nat_mat[idx,j] <- mean_val*shrink_percentage + nat_mat[idx,j]*(1-shrink_percentage)
    }
  }

  # manually force more correlation
  shrink_percentage <- 0.7 # higher means we shrink more. 0.4 "works" but then there is no difference visible
  for(j in 11:p){
    target_idx <- sample(intersect(1:10, which(y_block_assignment == y_block_assignment[j])),1)
    tmp_df <- data.frame(x = nat_mat[,target_idx], y = nat_mat[,j])
    lm_res <- stats::lm(y ~ x - 1, data = tmp_df)
    pred_y <- lm_res$fitted.values
    nat_mat[,j] <- pred_y*shrink_percentage + nat_mat[,j]*(1-shrink_percentage)
  }

  # manually force the off-genes to have no DE
  case_idx <- which(covariate[,"CC"] == 1)
  for(j in 11:p){
    mean_val_case <- mean(nat_mat[case_idx,j])
    mean_val_control <- mean(nat_mat[-case_idx,j])

    if(mean_val_control < mean_val_case){
      nat_mat[case_idx,j] - mean_val_case + mean_val_control
    } else {
      nat_mat[-case_idx,j] - mean_val_control + mean_val_case
    }
  }

  gamma_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    gamma_mat[,j] <- stats::rgamma(
      n = n,
      shape = exp(nat_mat[,j])*dispersion_vec[j],
      rate = dispersion_vec[j])
  }
  gamma_mat <- pmin(gamma_mat, 50)

  lib_mat <- tcrossprod(covariate[,c("Intercept", "Log_UMI", "Sex", "Age")], z_mat[,c("Intercept", "Log_UMI", "Sex", "Age")])
  lib_mat <- exp(lib_mat)
  obs_mat <- matrix(0, nrow = n, ncol = p)
  for(j in 1:p){
    obs_mat[,j] <- stats::rpois(n = n,
                                lambda = lib_mat[,j]*gamma_mat[,j])
  }
  covariate[,"Log_UMI"] <- log1p(Matrix::rowSums(obs_mat))

  rownames(obs_mat) <- paste0("c", 1:nrow(obs_mat))
  colnames(obs_mat) <- paste0("g", 1:ncol(obs_mat))
  rownames(covariate) <- rownames(obs_mat)

  seurat_obj <- Seurat::CreateSeuratObject(counts = t(obs_mat),
                                           meta.data = as.data.frame(covariate[,c("CC", "Sex", "Age", "Individual")]))

  list(covariate = covariate,
       seurat_obj = seurat_obj)
}

###########################

simulation_run_esvd <- function(covariate,
                                seurat_obj){
  set.seed(10)
  mat <- Matrix::t(seurat_obj[["RNA"]]@counts)
  covariate_dat <- seurat_obj@meta.data[,c("CC", "Sex", "Age", "Individual")]
  covariate_df <- data.frame(covariate_dat)
  covariates <- eSVD2:::format_covariates(dat = mat,
                                          covariate_df = covariate_df)
  covariate_df[,"Individual"] <- as.factor(covariate_df[,"Individual"])

  eSVD_obj <- eSVD2:::initialize_esvd(dat = mat,
                                      covariates = covariate[,c("Intercept", "CC", "Log_UMI", "Sex", "Age")],
                                      case_control_variable = "CC",
                                      bool_intercept = T,
                                      k = 2,
                                      lambda = 0.1,
                                      metadata_case_control = covariate[,"CC"],
                                      metadata_individual = covariate_df[,"Individual"],
                                      verbose = 0)

  eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Init",
    omitted_variables = "Log_UMI"
  )

  eSVD_obj <- eSVD2:::opt_esvd(
    input_obj = eSVD_obj,
    l2pen = 0.1,
    max_iter = 100,
    offset_variables = setdiff(colnames(eSVD_obj$covariates), "CC"),
    tol = 1e-6,
    verbose = 0,
    fit_name = "fit_First",
    fit_previous = "fit_Init")

  eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_First",
    omitted_variables = "Log_UMI"
  )

  eSVD_obj <- eSVD2:::opt_esvd(
    input_obj = eSVD_obj,
    l2pen = 0.1,
    max_iter = 100,
    offset_variables = NULL,
    tol = 1e-6,
    verbose = 0,
    fit_name = "fit_Second",
    fit_previous = "fit_First")

  eSVD_obj <- eSVD2:::.reparameterization_esvd_covariates(
    input_obj = eSVD_obj,
    fit_name = "fit_Second",
    omitted_variables = NULL
  )

  eSVD_obj <- eSVD2:::estimate_nuisance(input_obj = eSVD_obj,
                                        bool_covariates_as_library = T,
                                        bool_library_includes_interept = T,
                                        bool_use_log = F,
                                        verbose = 0)

  eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                        bool_adjust_covariates = F,
                                        alpha_max = NULL,
                                        bool_covariates_as_library = T,
                                        bool_stabilize_underdispersion = T,
                                        library_min = 1,
                                        pseudocount = 1)

  eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                             verbose = 0)

  df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
  teststat_vec <- eSVD_obj$teststat_vec
  p <- length(teststat_vec)
  gaussian_teststat <- sapply(1:p, function(j){
    qnorm(pt(teststat_vec[j], df = df_vec[j]))
  })

  multtest_res <- eSVD2:::multtest(gaussian_teststat)

  eSVD_obj$dat <- NULL
  list(eSVD_obj = eSVD_obj,
       gaussian_teststat = gaussian_teststat,
       multtest_res = multtest_res)
}
