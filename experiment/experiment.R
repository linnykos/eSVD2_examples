gene_name <- "AURKA"
covariate_individual <- "Subject_Identity"
input_obj <- eSVD_obj
metadata <- adams@meta.data

############

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

quantile(posterior_mean_mat[,gene_name])
quantile(posterior_var_mat[,gene_name])

#############

stopifnot(all(dim(posterior_mean_mat) == dim(posterior_var_mat)),
          covariate_individual %in% colnames(metadata),
          all(rownames(metadata) == rownames(posterior_mean_mat)))

p <- ncol(posterior_mean_mat)

if(verbose >= 1) print("Computing individual-level statistics")
tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                             control_individuals = control_individuals,
                                             covariate_individual = covariate_individual,
                                             metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- .construct_averaging_matrix(idx_list = all_indiv_idx,
                                       n = nrow(posterior_mean_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

avg_posterior_mean_mat[1:length(tmp$case_indiv_idx),gene_name]
avg_posterior_mean_mat[(length(tmp$case_indiv_idx)+1):nrow(avg_posterior_mean_mat),gene_name]

###

case_row_idx <- 1:length(case_individuals)
control_row_idx <- (length(case_individuals)+1):nrow(avg_posterior_mean_mat)
case_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[case_row_idx,,drop = F])
control_gaussian_mean <- Matrix::colMeans(avg_posterior_mean_mat[control_row_idx,,drop = F])
case_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[case_row_idx,,drop = F],
  avg_posterior_var_mat = avg_posterior_var_mat[case_row_idx,,drop = F]
)
control_gaussian_var <- eSVD2:::.compute_mixture_gaussian_variance(
  avg_posterior_mean_mat = avg_posterior_mean_mat[control_row_idx,,drop = F],
  avg_posterior_var_mat = avg_posterior_var_mat[control_row_idx,,drop = F]
)

case_gaussian_mean[gene_name]
control_gaussian_mean[gene_name]
sqrt(case_gaussian_var[gene_name])
sqrt(control_gaussian_var[gene_name])

n1 <- length(case_individuals)
n2 <- length(control_individuals)
teststat_vec <- (case_gaussian_mean - control_gaussian_mean) /
  (sqrt(case_gaussian_var/n1 + control_gaussian_var/n2))
names(teststat_vec) <- colnames(posterior_mean_mat)
teststat_vec[gene_name]

