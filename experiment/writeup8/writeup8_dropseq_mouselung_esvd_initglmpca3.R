rm(list=ls())
load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")

library(eSVD2); library(glmGamPoi); library(scran); library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat2 <- as.matrix(Matrix::t(mat))
nuisance_vec <- rep(glmpca_res$glmpca_family$nb_theta[1], ncol(mat2))

###################

Y <- mat
L <- 30
fam <- "nb"
minibatch <- "stochastic"
optimizer <- "avagrad"
ctl = list()
sz=NULL
nb_theta=NULL
X=NULL
Z=NULL
init=list(factors=NULL, loadings=NULL)

if(length(fam)>1){ fam<-fam[1] }

N<-ncol(Y); J<-nrow(Y)
#sanity check inputs
if(fam %in% c("poi","nb","nb2","binom")){ stopifnot(min(Y) >= 0) }

ic <- glmpca:::init_ctl(N, fam, minibatch, optimizer, ctl)
fam<-ic$fam; minibatch<-ic$minibatch; optimizer<-ic$optimizer; ctl<-ic$ctl

#create glmpca_family object
gnt <- glmpca:::glmpca_init(Y, fam, sz=sz, nb_theta=nb_theta)
gf<-gnt$gf; rfunc<-gnt$rfunc; offsets<-gnt$offsets

#initialize factors and loadings matrices
uv <- glmpca:::uv_init(N, J, L, gnt$intercepts, X=X, Z=Z, init=init)
U<-uv$U; V<-uv$V; lid<-uv$lid; uid<-uv$uid; vid<-uv$vid

#[[note: the intercept and its coefficient got bundled into U and V]]
#[[see U[,1] and V[,1] for this]]

ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("mat", "mat2",
                                "nuisance_vec",
                                "date_of_run",
                                "session_info",
                                "U", "V",
                                "offsets",
                                "lung")]
rm(list=ls_vec)

########################################

covariates <- matrix(1, nrow = nrow(mat2), ncol = 1)
colnames(covariates) <- c("Intercept")
b_init <- matrix(V[,1], nrow = nrow(V), ncol = 1)
colnames(b_init) <- c("Intercept")

K <- 30
n <- nrow(mat2)
p <- ncol(mat2)
set.seed(10)
x_mat <- matrix(stats::rnorm(n*K), nrow = n, ncol = K)
y_mat <- matrix(stats::rnorm(p*K), nrow = p, ncol = K)
tmp <- eSVD2:::.reparameterize(x_mat, y_mat, equal_covariance = T)
x_mat <- tmp$x_mat; y_mat <- tmp$y_mat
tmp <- tcrossprod(x_mat, y_mat) + tcrossprod(covariates, b_init)
x_mat <- x_mat/max(tmp)*3
y_mat <- y_mat/max(tmp)*3

print("Estimating NB via eSVD via Newton's")
time_start1 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(x_mat, y_mat, mat2,
                            family = "neg_binom2",
                            nuisance_param_vec = nuisance_vec,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = b_init,
                            covariates = covariates,
                            reestimate_nuisance = T,
                            global_estimate = T,
                            reparameterize = T,
                            offset_vec = offsets,
                            max_iter = 100,
                            tol = 1e-8,
                            verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_initglmpca3.RData")
