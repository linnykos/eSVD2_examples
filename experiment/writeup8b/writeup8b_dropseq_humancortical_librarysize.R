rm(list=ls())

library(Seurat)
load("../../../../out/writeup6b/writeup6b_dropseq_humancortical_glmpca_nb.RData")
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% "glmpca_res"]
rm(list = ls_vec)
load("../../../../out/writeup8b/writeup8b_dropseq_humancortical_esvd.RData")
source("plot_library_size.R")

umi_vec <- matrixStats::rowSums2(mat)
mat2 <- log(eSVD2:::.mult_vec_mat(1/umi_vec, mat)*1e4+1)
avg_expression <- matrixStats::colMeans2(mat2)

set.seed(10)
kmean_res <- stats::kmeans(avg_expression, centers = 6)
cbind(kmean_res$centers, kmean_res$size)
gene_grouping <- as.factor(kmean_res$cluster)

print("Smoothing original data")
smoothing_original <- smooth_gene_vs_umi(mat,
                                         gene_grouping)
save.image("../../../../out/writeup8b/writeup8b_dropseq_humancortical_librarysize.RData")

print("Smoothing eSVD fit 1")
pred_mat <- exp(tcrossprod(esvd_res$x_mat, esvd_res$y_mat))
smoothing_esvd <- smooth_gene_vs_umi(pred_mat,
                                     gene_grouping,
                                     umi_vec = umi_vec)
save.image("../../../../out/writeup8b/writeup8b_dropseq_humancortical_librarysize.RData")

print("Smoothing eSVD fit 3")
pred_mat <- exp(tcrossprod(esvd_res3$x_mat, esvd_res3$y_mat))
smoothing_esvd3 <- smooth_gene_vs_umi(pred_mat,
                                     gene_grouping,
                                     umi_vec = umi_vec)
save.image("../../../../out/writeup8b/writeup8b_dropseq_humancortical_librarysize.RData")

print("Smoothing SCTransform fit")
residual_mat <- t(as.matrix(cortical[["SCT"]]@data[Seurat::VariableFeatures(cortical, assay = "SCT"),]))
smoothing_sctransform <- smooth_gene_vs_umi(residual_mat,
                                            gene_grouping,
                                            umi_vec = umi_vec)
save.image("../../../../out/writeup8b/writeup8b_dropseq_humancortical_librarysize.RData")

print("Smoothing GLM-PCA")
pred_mat <- exp(tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)))
smoothing_glmpca <- smooth_gene_vs_umi(pred_mat,
                                       gene_grouping,
                                       umi_vec = umi_vec)
save.image("../../../../out/writeup8b/writeup8b_dropseq_humancortical_librarysize.RData")

