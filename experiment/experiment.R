eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj)
metadata <- sns@meta.data
metadata[,"individual"] <- as.factor(metadata[,"individual"])
covariate_individual = "individual"

input_obj <- eSVD_obj
case_control_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])

individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
stopifnot(length(intersect(control_individuals, case_individuals)) == 0)

latest_Fit <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
posterior_mean_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
posterior_var_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)
input_obj = posterior_mean_mat

verbose = 2
p <- ncol(posterior_mean_mat)
individual_stats <- lapply(1:p, function(j){
  if(verbose == 1 && p > 10 && j %% floor(p/10) == 0) cat('*')
  if(verbose >= 2) print(j)

  # next find the cells, then compute one gaussian per individual
  case_gaussians <- sapply(case_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata[,covariate_individual] == indiv)]
    cell_idx <- which(rownames(posterior_mean_mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata[,covariate_individual] == indiv)]
    cell_idx <- which(rownames(posterior_mean_mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  list(case_gaussians = case_gaussians,
       control_gaussians = control_gaussians)
})
