# Curved Gaussian example
library(eSVD2)
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

## Determine initialization
init <- eSVD2::initialize_esvd(dat, k = 2, family = "neg_binom2")


## Optimization
res <- eSVD2::opt_esvd(init$x_mat, init$y_mat, dat, family = "neg_binom2",
                       nuisance_param_vec = init$nuisance_param_vec,
                       library_size_vec = 1,
                       max_iter = 10,
                       verbose = 1)

