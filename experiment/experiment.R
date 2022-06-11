set.seed(123)
n <- 100
p <- 150
k <- 5
x_mat <- matrix(abs(rnorm(n * k))/10, nrow = n, ncol = k)
y_mat <- matrix(abs(rnorm(p * k))/10, nrow = p, ncol = k)
covariates <- cbind(sample(c(0,1), size = n, replace = T),
                    matrix(rnorm(n * 3, mean = 1, sd = 0.1), nrow = n, ncol = 3))
colnames(covariates) <- paste0("covariate_", 1:4)
z_mat <- cbind(c(rep(0, p/2), rep(100, p/2)), rep(10,p), rep(10,p), rep(1,p))
colnames(z_mat) <- paste0("covariate_", 1:4)
nat_mat <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, z_mat)/50

# Simulate data
dat <- generate_data(nat_mat,
                     family = "poisson",
                     nuisance_param_vec = NA,
                     library_size_vec = 1)
dat <- Matrix::Matrix(dat, sparse = T)
colnames(dat) <- paste0("g", 1:ncol(dat))
rownames(dat) <- paste0("c", 1:nrow(dat))

# res <- .initialize_coefficient(bool_intercept = T,
#                                case_control_variable = "covariate_1",
#                                covariates = covariates,
#                                dat = dat,
#                                lambda = 0.1,
#                                mixed_effect_variables = c("covariate_2", "covariate_3"),
#                                offset_variables = NULL)

#############
bool_intercept = T
case_control_variable = "covariate_1"
lambda = 0.1
mixed_effect_variables = c("covariate_2", "covariate_3")
offset_variables = NULL
verbose = 2

n <- nrow(dat); p <- ncol(dat)
covariates_nooffset <- covariates[,which(!colnames(covariates) %in% offset_variables)]
offset_vec <- Matrix::rowSums(covariates[,offset_variables,drop = F])

log_pval <- rep(NA, p)
names(log_pval) <- colnames(dat)

if(bool_intercept){
  colnames_vec <- c("Intercept", colnames(covariates_nooffset))
} else {
  colnames_vec <- colnames(covariates_nooffset)
}

z_mat1 <- matrix(NA, nrow = p, ncol = length(colnames_vec))
z_mat2 <- matrix(NA, nrow = p, ncol = length(colnames_vec))
colnames(z_mat1) <- colnames_vec
colnames(z_mat2) <- colnames_vec
rownames(z_mat2) <- colnames(dat)
rownames(z_mat1) <- colnames(dat)

# for(j in 1:p){
#   print(j)
#   if(verbose == 1 && p >= 10 && j %% floor(p/10) == 0) cat('*')
#   if(verbose >= 2) print(paste0("Finished variable ", j , " of ", p))
#   tmp <- .lrt_coefficient(bool_intercept = bool_intercept,
#                           case_control_variable = case_control_variable,
#                           covariates = covariates_nooffset,
#                           lambda = lambda,
#                           mixed_effect_variables = mixed_effect_variables,
#                           offset_vec = offset_vec,
#                           vec = as.numeric(dat[,j]),
#                           verbose = verbose)
#   z_mat1[j,] <- tmp$coef_vec1
#   z_mat2[j,] <- tmp$coef_vec2
#   log_pval[j] <- tmp$log_pval
# }

##############

covariates = covariates_nooffset
vec = as.numeric(dat[,j])

penalty_factor1 <- rep(0, ncol(covariates))
penalty_factor1[colnames(covariates) %in% mixed_effect_variables] <- 1
glm_fit1 <- glmnet::glmnet(x = covariates,
                           y = vec,
                           alpha = 0,
                           family = "poisson",
                           intercept = bool_intercept,
                           lambda = exp(seq(log(1e4), log(lambda), length.out = 100)),
                           offset = offset_vec,
                           penalty.factor = penalty_factor1,
                           standardize = F)

if(bool_intercept){
  coef_vec1 <- c(glm_fit1$a0[length(glm_fit1$a0)], glm_fit1$beta[,ncol(glm_fit1$beta)])
  mean_vec1 <- exp(covariates %*% coef_vec1[-1] + offset_vec + coef_vec1[1])
} else {
  coef_vec1 <- glm_fit1$beta[,ncol(glm_fit1$beta)]
  mean_vec1 <- exp(covariates %*% coef_vec1 + offset_vec)
}

log_vec <- vec/mean_vec1; log_vec[vec != 0] <- log(log_vec[vec != 0])
deviance1 <- 2*sum(vec*log_vec - (vec - mean_vec1))

covariates2 <- covariates[,which(colnames(covariates) != case_control_variable), drop = F]
penalty_factor2 <- rep(0, ncol(covariates2))
penalty_factor2[colnames(covariates2) %in% mixed_effect_variables] <- 1
glm_fit2 <- glmnet::glmnet(x = covariates2,
                           y = vec,
                           alpha = 0,
                           family = "poisson",
                           intercept = bool_intercept,
                           lambda = exp(seq(log(1e4), log(lambda), length.out = 100)),
                           offset = offset_vec,
                           penalty.factor = penalty_factor2,
                           standardize = F)
if(bool_intercept){
  coef_vec2 <- c(glm_fit2$a0[length(glm_fit2$a0)], glm_fit2$beta[,ncol(glm_fit2$beta)])
  mean_vec2 <- exp(covariates2 %*% coef_vec2[-1] + offset_vec + coef_vec2[1])
} else {
  coef_vec2 <- glm_fit2$beta[,ncol(glm_fit2$beta)]
  mean_vec2 <- exp(covariates2 %*% coef_vec2 + offset_vec)
}
log_vec <- vec/mean_vec2; log_vec[vec != 0] <- log(log_vec[vec != 0])
deviance2 <- 2*sum(vec*log_vec - (vec - mean_vec2))

residual_deviance <- max(deviance2 - deviance1, 0)
log_pval <- stats::pchisq(residual_deviance, df = 1,
                          lower.tail = FALSE,
                          log.p = T)

