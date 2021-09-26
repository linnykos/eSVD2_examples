rm(list=ls())

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
load("../../../../data/10x_mouseretinal/10x_mouseretinal_formatted.RData")

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

retinal <- Seurat::CreateSeuratObject(counts = Matrix::t(dat),
                                      meta.data = metadata, min.cells = 5)
largest_batch <- names(which.max(table(retinal@meta.data$BatchID)))
cell_keep <- rownames(retinal@meta.data[retinal@meta.data$BatchID == largest_batch,])
retinal <- retinal[,cell_keep]
retinal <- Seurat::NormalizeData(retinal, normalization.method = "LogNormalize", scale.factor = 10000)
retinal <-  Seurat::FindVariableFeatures(retinal, selection.method = "vst", nfeatures = 2000)

mat <- retinal[["RNA"]]@counts[Seurat::VariableFeatures(retinal),]
mat <- Matrix::t(mat)
mat <- as.matrix(mat)

##############
# initialization
K <- 30
n <- nrow(mat)
p <- ncol(mat)

covariates <- cbind(1, log(matrixStats::rowMeans2(mat)))
colnames(covariates) <- c("Intercept", "Log-UMI")

init_res <- eSVD2::initialize_esvd(mat,
                                   k = K,
                                   family = "neg_binom2",
                                   covariates = covariates,
                                   offset_vec = rep(0, nrow(mat)))

###################3

print("Estimating NB via eSVD, round 1")
time_start1 <- Sys.time()
set.seed(10)
esvd_res <- eSVD2::opt_esvd(init_res$x_mat,
                            init_res$y_mat,
                            mat,
                            family = "neg_binom2",
                            nuisance_param_vec = init_res$nuisance_param_vec,
                            library_size_vec = 1,
                            method = "newton",
                            b_init = init_res$b_mat,
                            covariates = init_res$covariates,
                            offset_vec = init_res$offset_vec,
                            reestimate_nuisance = T,
                            global_estimate = T,
                            reparameterize = T,
                            max_iter = 100,
                            tol = 1e-8,
                            verbose = 1)
time_end1 <- Sys.time()
save.image("../../../../out/writeup8c/writeup8c_10x_mouseretinal_esvd.RData")

print("Estimating NB via eSVD, round 2")
time_start2 <- Sys.time()
set.seed(10)
esvd_res2 <- eSVD2::opt_esvd(esvd_res$x_mat,
                             esvd_res$y_mat, mat,
                             family = "neg_binom2",
                             nuisance_param_vec = esvd_res$nuisance_param_vec,
                             library_size_vec = 1,
                             method = "newton",
                             b_init = esvd_res$b_mat,
                             covariates = esvd_res$covariates,
                             offset_vec = esvd_res$offset_vec,
                             reestimate_nuisance = T,
                             global_estimate = F,
                             reparameterize = T,
                             max_iter = 100,
                             tol = 1e-8,
                             verbose = 1)
time_end2 <- Sys.time()
save.image("../../../../out/writeup8c/writeup8c_10x_mouseretinal_esvd.RData")
