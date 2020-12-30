rm(list=ls())

# generate data
set.seed(123)
n <- 100
p <- 150
k <- 5
x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
y_mat <- -matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
nat_mat <- tcrossprod(x_mat, y_mat)
library_size_vec <- 1:n
nuisance_param_vec <- c(1:p)*10

dat <- eSVD2::generate_data(
  nat_mat, family = "neg_binom", nuisance_param_vec = nuisance_param_vec,
  library_size_vec = library_size_vec
)

# image(log(dat+1))

####################

res <- eSVD2::initialize_esvd(dat, k = 5, family = "poisson",
                              library_size_vec = library_size_vec)
res2 <- eSVD2::opt_esvd(res$x_mat, res$y_mat, dat, family = "poisson",
                 library_size_vec = res$library_size_vec)
