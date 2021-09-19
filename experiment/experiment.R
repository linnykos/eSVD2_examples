set.seed(123)
n <- 10
p <- 15
k <- 2
nuisance_param_vec <- runif(p, 0, 10)
x_mat <- matrix(rnorm(n * k), nrow = n, ncol = k)
y_mat <- matrix(rnorm(p * k), nrow = p, ncol = k)
nat_mat <- tcrossprod(x_mat, y_mat)

dat <- eSVD2::generate_data(
  nat_mat, family = "neg_binom2", nuisance_param_vec = nuisance_param_vec,
  library_size_vec = 1
)

gamma = nuisance_param_vec
A <- dat
theta <- x_mat%*%t(y_mat)
s <- NULL
zz1 <- .log_prob.neg_binom2(A, theta, s, gamma)
i <- 5
zz2 <- .log_prob_row.neg_binom2(A[i,], theta[i,], si = NULL, gamma)
cbind(zz1[i,], zz2)
j <- 5
zz3 <- .log_prob_col.neg_binom2(A[,j], thetaj = theta[,j], s = NULL, gamma[j])
cbind(zz1[,j], zz3)
