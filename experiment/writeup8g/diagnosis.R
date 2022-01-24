rm(list=ls())
load("../../out/writeup8g/tmp.RData")

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

idx <- which(colnames(mat) == "POLA2")
plot(mean_mat[,idx], mat[,idx], asp = T)
plot(mean_mat_nolib[,idx]*library_mat[,idx], mat[,idx], asp = T)

mean_vec <- mean_mat_nolib[,idx]*library_mat[,idx]
var_vec <- (mat[,idx] - mean_vec)^2
plot(mean_vec, var_vec, asp = T)
tmp_df <- data.frame(mean = mean_vec, var = var_vec)
tmp_lm <- stats::lm(var ~ mean - 1, data = tmp_df)
stats::coef(tmp_lm)

res_list[1,idx]

hist(mat[,idx], breaks = 50)
n <- nrow(mat)
val <- res_list[1,idx]
set.seed(10)
zz <- stats::rgamma(n, shape = val*mean_mat_nolib[,idx], rate = val)
yy <- stats::rpois(n, lambda = zz*library_mat[,idx])
hist(yy, breaks = 50)
table(mat[,idx])
table(yy)

############################

mu <- mean_mat_nolib[,idx]
sf <- library_mat[,idx]
y <- mat[,idx]
calculate_fano_parameter(y, mu, sf)

val <- 61189
set.seed(10)
zz <- stats::rgamma(n, shape = val*mean_mat_nolib[,idx], rate = val)
yy <- stats::rpois(n, lambda = zz*library_mat[,idx])
hist(yy, breaks = 50)
table(mat[,idx])
table(yy)
hist(zz, breaks = 50)
hist(library_mat[,idx], breaks = 50)

####################################
####################################
####################################

idx <- which(colnames(mat) == "GALNTL6")
hist(mat[,idx], breaks = 50)
n <- nrow(mat)
val <- res_list[1,idx]
set.seed(10)
zz <- stats::rgamma(n, shape = val*mean_mat_nolib[,idx], rate = val)
yy <- stats::rpois(n, lambda = zz*library_mat[,idx])
hist(yy, breaks = 50)


mu <- mean_mat_nolib[,idx]
sf <- library_mat[,idx]
y <- mat[,idx]
calculate_fano_parameter(y, mu, sf)
