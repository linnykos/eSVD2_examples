rm(list=ls())
set.seed(20)
cov_x <- cov(MASS::mvrnorm(n = 10, rep(0, 5), diag(5)))
cov_y <- cov(MASS::mvrnorm(n = 10, rep(0, 5), toeplitz(5:1)))

res <- .identification(cov_x, cov_y, check = T)

##############
check = T; tol = 1e-6

eigen_x <- eigen(cov_x)
eigen_y <- eigen(cov_y)

Vx <- eigen_x$vectors
Vy <- eigen_y$vectors

if(any(eigen_x$values <= tol) | any(eigen_y$values <= tol)) warning("Detecting rank defficiency in reparameterization step")

Dx <- eigen_x$values
Dy <- eigen_y$values

# form R
tmp2 <- diag(sqrt(Dy)) %*% t(Vy) %*% Vx %*% diag(sqrt(Dx))
tmp <- crossprod(.mult_mat_vec(Vy, sqrt(Dy)), .mult_mat_vec(Vx, sqrt(Dx)))
stopifnot(sum(abs(tmp2 - tmp)) <= 1e-6)
svd_tmp <- svd(tmp)
Q <- tcrossprod(svd_tmp$u, svd_tmp$v)

# run a check
if(check){
  sym_mat <- crossprod(Q, tmp)

  stopifnot(sum(abs(sym_mat - t(sym_mat))) <= 1e-6)
}
