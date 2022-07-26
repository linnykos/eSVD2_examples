intersect(which(x_vec_noninflamed <= -2), which(logpvalue_vec_noninflamed <= 1))
teststat_vec[37]

input_obj <- eSVD_obj_noninflamed
metadata <- regevEpi_noninflamed@meta.data
covariate_individual <- "Sample"

case_control_variable <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "init_case_control_variable", which_fit = "param")
covariates <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "covariates", which_fit = NULL)
cc_vec <- covariates[,case_control_variable]
cc_levels <- sort(unique(cc_vec), decreasing = F)
stopifnot(length(cc_levels) == 2)
control_idx <- which(cc_vec == cc_levels[1])
case_idx <- which(cc_vec == cc_levels[2])

latest_Fit <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "latest_Fit", which_fit = NULL)
posterior_mean_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_mean_mat", which_fit = latest_Fit)
posterior_var_mat <- eSVD2:::.get_object(eSVD_obj = input_obj, what_obj = "posterior_var_mat", which_fit = latest_Fit)

individual_vec <- metadata[,covariate_individual]
control_individuals <- unique(individual_vec[control_idx])
case_individuals <- unique(individual_vec[case_idx])
tmp <- eSVD2:::.determine_individual_indices(case_individuals = case_individuals,
                                             control_individuals = control_individuals,
                                             covariate_individual = covariate_individual,
                                             metadata = metadata)
all_indiv_idx <- c(tmp$case_indiv_idx, tmp$control_indiv_idx)
avg_mat <- eSVD2:::.construct_averaging_matrix(idx_list = all_indiv_idx,
                                               n = nrow(posterior_mean_mat))
avg_posterior_mean_mat <- as.matrix(avg_mat %*% posterior_mean_mat)
avg_posterior_var_mat <- as.matrix(avg_mat %*% posterior_var_mat)

#######################

y1 <- logpvalue_vec_noninflamed[unique(c(idx_inflamed, idx_noninflamed))]
y2 <- logpvalue_vec_inflamed[unique(c(idx_inflamed, idx_noninflamed))]

idx <- unique(c(which(is.infinite(y1)), which(is.infinite(y2))))
y1 <- y1[-idx]; y2 <- y2[-idx]

y1 <- eSVD_obj_inflamed$teststat_vec[unique(c(idx_inflamed, idx_noninflamed))]
y2 <- eSVD_obj_noninflamed$teststat_vec[unique(c(idx_inflamed, idx_noninflamed))]
stats::cor(y1, y2)

y1 <- gaussian_teststat_inflamed[unique(c(idx_inflamed, idx_noninflamed))]
y2 <- gaussian_teststat_noninflamed[unique(c(idx_inflamed, idx_noninflamed))]
stats::cor(y1, y2)

y1 <- gaussian_teststat_inflamed[unique(c(idx_inflamed, idx_noninflamed))]
y2 <- gaussian_teststat_noninflamed[unique(c(idx_inflamed, idx_noninflamed))]
y1 <- y1 - locfdr_res_inflamed$fp0["mlest", "delta"]
y2 <- y2 - locfdr_res_noninflamed$fp0["mlest", "delta"]
sum(y1*y2)/sqrt(sum(y1^2)*sum(y2^2))

y1 <- gaussian_teststat_inflamed[unique(c(hk_idx))]
y2 <- gaussian_teststat_noninflamed[unique(c(hk_idx))]
stats::cor(y1, y2)
y1 <- y1 - locfdr_res_inflamed$fp0["mlest", "delta"]
y2 <- y2 - locfdr_res_noninflamed$fp0["mlest", "delta"]
sum(y1*y2)/sqrt(sum(y1^2)*sum(y2^2))

hk_idx2 <- setdiff(hk_idx, c(idx_inflamed, idx_noninflamed))
y1 <- gaussian_teststat_inflamed[unique(c(hk_idx2))]
y2 <- gaussian_teststat_noninflamed[unique(c(hk_idx2))]
stats::cor(y1, y2)
y1 <- y1 - locfdr_res_inflamed$fp0["mlest", "delta"]
y2 <- y2 - locfdr_res_noninflamed$fp0["mlest", "delta"]
sum(y1*y2)/sqrt(sum(y1^2)*sum(y2^2))

y1 <- gaussian_teststat_inflamed[unique(c(cycling_idx))]
y2 <- gaussian_teststat_noninflamed[unique(c(cycling_idx))]
stats::cor(y1, y2)
y1 <- y1 - locfdr_res_inflamed$fp0["mlest", "delta"]
y2 <- y2 - locfdr_res_noninflamed$fp0["mlest", "delta"]
sum(y1*y2)/sqrt(sum(y1^2)*sum(y2^2))


###################

y1 <- gaussian_teststat_adams[unique(c(idx_adams, idx_habermann))]
y2 <- gaussian_teststat_habermann[unique(c(idx_adams, idx_habermann))]
stats::cor(y1, y2)
y1 <- y1 - locfdr_res_adams$fp0["mlest", "delta"]
y2 <- y2 - locfdr_res_habermann$fp0["mlest", "delta"]
sum(y1*y2)/sqrt(sum(y1^2)*sum(y2^2))

y1 <- gaussian_teststat_adams[hk_idx]
y2 <- gaussian_teststat_habermann[hk_idx]
stats::cor(y1, y2)
y1 <- y1 - locfdr_res_adams$fp0["mlest", "delta"]
y2 <- y2 - locfdr_res_habermann$fp0["mlest", "delta"]
sum(y1*y2)/sqrt(sum(y1^2)*sum(y2^2))

