# Curved Gaussian example
library(eSVD2)
set.seed(123)
n <- 100
p <- 150
k <- 5
x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
nat_mat <- tcrossprod(x_mat, y_mat)

## Simulate data
dat <- eSVD2::generate_data(nat_mat, family = "curved_gaussian", nuisance_param_vec = 2,
                            library_size_vec = 1, tol = 1e-3)

## Determine initialization
init <- eSVD2::initialize_esvd(
    dat, k = k, family = "curved_gaussian", nuisance_param_vec = 2, library_size_vec = 1,
    config = eSVD2::initialization_options()
)

## Optimization
res <- opt_esvd(init$x_mat, init$y_mat, dat, family = "curved_gaussian",
                nuisance_param_vec = 2, library_size_vec = 1,
                max_iter = 10, verbose = 1)
nat_est <- tcrossprod(res$x, res$y)
plot(res$loss)



# Exponential example
library(eSVD2)
set.seed(123)
n <- 100
p <- 150
k <- 5
x_mat <- matrix(abs(rnorm(n * k)), nrow = n, ncol = k)
y_mat <- matrix(abs(rnorm(p * k)), nrow = p, ncol = k)
nat_mat <- -tcrossprod(x_mat, y_mat) / 10

## Simulate data
dat <- eSVD2::generate_data(nat_mat, family = "exponential", nuisance_param_vec = NA,
                            library_size_vec = 1)

## Determine initialization
init <- eSVD2::initialize_esvd(
    dat, k = k, family = "exponential", nuisance_param_vec = NA, library_size_vec = 1,
    config = eSVD2::initialization_options()
)

## Optimization
res <- opt_esvd(init$x_mat, init$y_mat, dat, family = "exponential",
                nuisance_param_vec = NA, library_size_vec = 1,
                max_iter = 15, verbose = 1)
nat_est <- tcrossprod(res$x, res$y)
plot(res$loss)
