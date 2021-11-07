rm(list=ls())
load("../../../../out/writeup8d/writeup8d_sns_pseudoreal.RData")

x_init = init_res$x_mat
y_init = init_res$y_mat
dat = mat
family = "neg_binom2"
nuisance_param_vec = init_res$nuisance_param_vec
library_size_vec = 1
method = "newton"
b_init = init_res$b_mat
covariates = init_res$covariates
offset_vec = init_res$offset_vec
reestimate_nuisance = T
global_estimate = T
reparameterize = T
max_iter = 50
l2pen = 0.1
verbose = 1

n <- nrow(dat)
p <- ncol(dat)
k <- ncol(x_init)

storage.mode(x_init) <- "double"
storage.mode(y_init) <- "double"
storage.mode(dat) <- "double"

family <- eSVD2:::.string_to_distr_funcs(family)
library_size_vec <- eSVD2:::.parse_library_size(dat, library_size_vec)
if(all(!is.na(nuisance_param_vec)) && length(nuisance_param_vec) == 1) {
  nuisance_param_vec <- rep(nuisance_param_vec[1], ncol(dat))
}
library_size_vec <- as.numeric(library_size_vec)
nuisance_param_vec <- as.numeric(nuisance_param_vec)

opt_fun <- if(method == "newton") eSVD2:::constr_newton else eSVD2:::constr_lbfgs

r <- ncol(covariates)
b_mat <- matrix(0, p, r)
x_mat <- x_init
y_mat <- y_init
losses <- c()
i =1

##############################

X0 = x_mat
Y = y_mat
B = b_mat
Z = covariates
A = dat
family = family
s = library_size_vec
gamma = nuisance_param_vec
i = 1
Zi <- if(is.null(Z)) NULL else Z[i, ]

opt <- opt_fun(
  x0 = X0[i, ],
  f = eSVD2:::objfn_Xi,
  gr = eSVD2:::grad_Xi,
  hn = eSVD2:::hessian_Xi,
  direc = eSVD2:::direction_Xi,
  feas = eSVD2:::feas_Xi,
  eps_rel = 1e-3,
  verbose = (verbose >= 3),
  Y = Y, B = B, Zi = Zi, Ai = A[i, ],
  family = family,
  si = s[i],
  gamma = gamma,
  offseti = offset_vec[i],
  l2pen = l2pen)

opt

