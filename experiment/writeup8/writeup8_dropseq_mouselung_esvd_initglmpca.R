rm(list=ls())
load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")

library(eSVD2); library(glmGamPoi); library(scran); library(Seurat)
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

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
fam <- match.arg(fam)
minibatch <- match.arg(minibatch)
optimizer <- match.arg(optimizer)

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
