rm(list=ls())

library(Seurat)
load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData") # because we forgot to save lung
load("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_initglmpca.RData")

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                         key = "esvdfactorumap_",
                                                         assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2, initialized w/ GLM-PCA), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca_factor_esvd.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

##########

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["glmpcafactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                           key = "glmpcafactorumap_",
                                                           assay = "RNA")

plot1 <- Seurat::DimPlot(lung,
                         reduction = "glmpcafactorumap",
                         group.by = "celltype",
                         label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\nGLM-PCA (initialized w/ eSVD), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca_factor_glmpca.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

#################################

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
mean(eSVD2:::.log_prob.neg_binom2(mat2,
                                  theta = nat_mat,
                                  s = NULL,
                                  gamma = nuisance_vec))
quantile(esvd_res$b_mat[,1])
quantile(esvd_res$b_mat[,2])
nuisance_vec[1]

nat_mat <- tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)) + tcrossprod(as.matrix(glmpca_res$X), as.matrix(glmpca_res$coefX))
nat_mat <- sweep(nat_mat, 1, glmpca_res$offsets, "+")
mean(eSVD2:::.log_prob.neg_binom2(mat2,
                                  theta = nat_mat,
                                  s = NULL,
                                  gamma = rep(glmpca_res$glmpca_family$nb_theta[1], ncol(mat2))))
quantile(glmpca_res$coefX[,1])
glmpca_res$glmpca_family$nb_theta

#################

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
mean_mat <- exp(nat_mat)
j <- 10
MASS::theta.ml(y = mat2[,j], mu = mean_mat[,j])

# the following line is too slow
# MASS::theta.ml(y = as.numeric(mat2), mu = as.numeric(mean_mat))


zz <- glmGamPoi::overdispersion_mle(y = t(mat2),
                                    mean = t(mean_mat),
                                    global_estimate = F,
                                    verbose = T)


