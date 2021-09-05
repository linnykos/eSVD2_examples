rm(list=ls())
load("2021-09-05-nb_hessian.RData")
yb_mat <- cbind(y_mat, b_mat)
YB0 <- yb_mat
X <- x_mat
Z <- covariates
A <- dat
family <- eSVD2:::.string_to_distr_funcs("neg_binom")
s <- library_size_vec
gamma <- nuisance_param_vec
opt_fun <- eSVD2:::constr_newton

j <- 73
p <- ncol(A)
YB <- YB0
XZ <- cbind(X, Z)
# the following line gives the error:
# Newton iter = 0, fx = 0.114529, ||grad|| = 85.474231
# Error in solve.default(Hess, grad) :
#   system is computationally singular: reciprocal condition number = 1.61942e-16
opt <- opt_fun(
  x0 = YB0[j, ],
  f = eSVD2:::objfn_Yj,
  gr = eSVD2:::grad_Yj,
  hn = eSVD2:::hessian_Yj,
  feas = eSVD2:::feas_Yj,
  eps_rel = 1e-3,
  verbose = 3,
  X = XZ,
  Bj = NULL,
  Z = NULL,
  Aj = A[, j],
  family = family,
  s = s,
  gammaj = gamma[j])

#######

# investigate the inputs
quantile(matrixStats::colSds(X))
table(A[,j])
gamma[j]

# to investigate the hessian:
hess_mat <- eSVD2:::hessian_Yj(Yj = YB0[j, ],
                               X = XZ,
                               Bj = NULL,
                               Z = NULL,
                               Aj = A[, j],
                               family = family,
                               s = s,
                               gammaj = gamma[j])
eigen(hess_mat)$values # we see a lot of very small values

