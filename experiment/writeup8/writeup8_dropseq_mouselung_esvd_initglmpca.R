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
                                "offsets")]

########################################

covariates <- cbind(1, offsets)
colnames(covariates) <- c("Intercept", "GLMPCA-offsets")
b_init <- cbind(V[,1], 1)
colnames(b_init) <- c("Intercept", "GLMPCA-offsets")

print("Estimating NB via eSVD")
time_start3 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(U[,-1], V[,-1], mat2,
                            family = "neg_binom2",
                            nuisance_param_vec = nuisance_vec,
                            library_size_vec = 1,
                            b_init = b_init,
                            covariates = covariates,
                            max_iter = 100,
                            tol = 1e-8,
                            verbose = 1)
time_end3 <- Sys.time()
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_initglmpca.RData")

###################

print("Estimating NB via GLM-PCA")
set.seed(10)
K <- 30
time_start4 <- Sys.time()
glmpca_res <- glmpca::glmpca(mat, L = K, fam = "nb",
                             ctl = list(verbose = T),
                             minibatch = "stochastic")
time_end4 <- Sys.time()
print("Finished")
save.image("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_initglmpca.RData")

