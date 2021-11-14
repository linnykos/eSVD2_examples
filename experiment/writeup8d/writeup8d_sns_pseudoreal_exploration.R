rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_pseudoreal_data.RData")
# load("../../out/writeup8d/writeup8d_sns_pseudoreal_data.RData") # for local mac

true_mean <- eSVD2:::.mult_vec_mat(s_vec, lambda_mat)

quantile_vec <- sapply(1:ncol(mat), function(j){
  residual <- true_mean[,j] - mat[,j]
  idx <- which(abs(residual) <= 3*sqrt(true_mean[,j]))
  length(idx)/nrow(mat)
})

quantile(quantile_vec)
quantile(true_esvd$nuisance_param_vec)

######################

j <- 8
true_esvd$nuisance_param_vec[j]
residual <- lambda_mat[,j] - exp(nat_mat[,j])
plot(exp(nat_mat[,j]), residual, asp = T)
plot(exp(nat_mat[,j]), residual)
plot(exp(nat_mat[,j]), lambda_mat[,j], asp = T)
plot(exp(nat_mat[,j]), lambda_mat[,j])

hist(lambda_mat[,j])

######################

j <- 1761
true_esvd$nuisance_param_vec[j]
plot(exp(nat_mat[,j]), lambda_mat[,j])

######################

j <- which.min(abs(true_esvd$nuisance_param_vec - 0.6))
true_esvd$nuisance_param_vec[j]
plot(exp(nat_mat[,j]), lambda_mat[,j])

######################

n <- nrow(lambda_mat)
set.seed(10)
vec <- sapply(1:n, function(i){
  stats::rgamma(1,
                shape = true_esvd$nuisance_param_vec[j],
                scale = exp(nat_mat[i,j])/true_esvd$nuisance_param_vec[j])

})
plot(exp(nat_mat[,j]), vec)

##################

x_vec <- seq(0, 500, length.out = 1000)
shape <- 10
y_vec <- stats::dgamma(x_vec, shape = shape,
                       scale = 500/shape)
plot(x_vec,y_vec)

shape <- 0.5
scale <- 10
n <- 10000
set.seed(10)
samples <- stats::rgamma(n, shape = shape, scale = scale/shape)
length(which(samples <= 1))/10000
hist(samples)
plot(samples)
mean(samples)
var(samples)
500^2/shape

##################################3

set.seed(10)
vec <- stats::rnbinom(10000, size = 1, mu = 10)
hist(vec)

vec <- stats::rnbinom(10000, size = 0.1, mu = 10)
hist(vec)



