library(NMF)
library(dimRed)
library(destiny)
library(Rtsne)
library(umap)
library(pCMF)
library(SummarizedExperiment)
library(zinbwave)

method_svd <- function(dat, k = 2){
  stopifnot(k >= 2)
  tmp <- svd(dat)
  fit <- tmp$u[,1:k] %*% diag(sqrt(tmp$d[1:k]))

  list(fit = fit)
}

method_tsne_oracle <- function(dat, cell_truth, paramMat, k = 2){
  fit_list <- lapply(1:nrow(paramMat), function(i){
    set.seed(10)
    tmp <- Rtsne::Rtsne(dat, dims = k, perplexity = paramMat[i, "perplexity"])
    tmp$Y
  })

  dist_mat_truth <- as.matrix(stats::dist(cell_truth))

  quality_vec <- sapply(fit_list, function(fit){
    dist_mat_est <- as.matrix(stats::dist(fit[,1:2]))

    mean(sapply(1:nrow(dist_mat_est), function(i){
      cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
    }))
  })

  idx <- which.max(quality_vec)

  fit <- fit_list[[idx]]

  list(fit = fit, perplexity = paramMat[idx, "perplexity"])
}

method_umap_oracle <- function(dat, cell_truth, paramMat, k = 2){
  fit_list <- lapply(1:nrow(paramMat), function(i){
    custom.settings <- umap::umap.defaults
    custom.settings$n_components <- k
    custom.settings$n_neighbors <- paramMat[i, "n_neighbors"]
    custom.settings$min_dist <- paramMat[i, "min_dist"]
    custom.settings$init <- "random"

    set.seed(10)
    tmp <- umap::umap(dat, config = custom.settings)
    tmp$layout
  })

  dist_mat_truth <- as.matrix(stats::dist(cell_truth))

  quality_vec <- sapply(fit_list, function(fit){
    dist_mat_est <- as.matrix(stats::dist(fit[,1:2]))

    mean(sapply(1:nrow(dist_mat_est), function(i){
      cor(dist_mat_truth[i,], dist_mat_est[i,], method = "kendall")
    }))
  })

  idx <- which.max(quality_vec)

  fit <- fit_list[[idx]]

  list(fit = fit, min_dist = paramMat[idx, "min_dist"], n_neighbors = paramMat[idx, "n_neighbors"])
}

method_zinbwave <- function(dat, k = 2){
  set.seed(10)
  dat_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = t(dat)))
  tmp <- zinbwave::zinbwave(dat_se, K = 2, maxiter.optimize = 100, normalizedValues = F,
                            commondispersion = F)
  fit <- SingleCellExperiment::reducedDims(tmp)$zinbwave

  list(fit = fit)
}

method_pcmf <- function(dat, k = 2){
  set.seed(10)
  tmp <- pCMF::pCMF(dat, K = 2, sparsity = F, verbose = F)
  fit <- tmp$factor$U

  list(fit = fit)
}

method_isomap <- function(dat, k = 2){
  set.seed(10)
  dimRed_obj <- dimRed::dimRedData(dat)
  isomap_obj <- dimRed::Isomap()
  isomap_obj@stdpars$ndim <- k
  isomap_obj@stdpars$knn <- round(nrow(dat)/10)
  suppressMessages(emb <- isomap_obj@fun(dimRed_obj, pars = isomap_obj@stdpars))

  fit <- emb@data@data

  list(fit = fit)
}

# requires the "fastICA" package
method_ica <- function(dat, k = 2){
  set.seed(10)
  dimRed_obj <- dimRed::dimRedData(dat)
  fastica_obj <- dimRed::FastICA()
  fastica_obj@stdpars$ndim <- k
  suppressMessages(emb <- fastica_obj@fun(dimRed_obj, pars = fastica_obj@stdpars))

  fit <- emb@data@data

  list(fit = fit)
}

# requires the "NMF" package -- be sure to explicitly load this package prior to usage
# https://support.bioconductor.org/p/110844/ -- be sure to install the version from the below line
## to avoid a conflict
# BiocManager::install("renozao/NMF", ref = "devel")
method_nmf <- function(dat, k = 2){
  set.seed(10)
  dimRed_obj <- dimRed::dimRedData(dat)
  nnmf_obj <- dimRed::NNMF()
  nnmf_obj@stdpars$ndim <- k
  suppressMessages(emb <- nnmf_obj@fun(dimRed_obj, pars = nnmf_obj@stdpars))

  fit <- emb@data@data

  list(fit = fit)
}

# requires the "diffusionMap" function
method_diffusion <- function(dat, k = 2){
  set.seed(10)
  tmp <- destiny::DiffusionMap(dat)

  fit <- cbind(tmp$DC1, tmp$DC2)

  list(fit = fit)
}



