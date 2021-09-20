rm(list=ls())

library(Seurat)
load("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb2_glmpca.RData")
quantile(nuisance_vec)
quantile(esvd_res$b_mat[,1])
quantile(esvd_res$b_mat[,2])

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2, initialized w/ GLM-PCA), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_factor_nb2_glmpca_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


####

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
pred_mat <- exp(nat_mat)
rownames(pred_mat) <- rownames(mat)
colnames(pred_mat) <- colnames(mat)

glm_assay <- Seurat::CreateAssayObject(counts = t(pred_mat))
lung[["pred"]] <- glm_assay

Seurat::DefaultAssay(lung) <- "pred"
lung <- Seurat::NormalizeData(lung, normalization.method = "LogNormalize", scale.factor = 10000)
lung[["pred"]]@var.features <- rownames(lung)
lung <-  Seurat::ScaleData(lung)
lung <- Seurat::RunPCA(lung, features = Seurat::VariableFeatures(lung), verbose = F,
                       reduction.name = "esvdpca")

set.seed(10)
lung <- Seurat::RunUMAP(lung, reduction = "esvdpca", dims = 1:30, reduction.name = "esvdumap")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdpca", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2, initialized w/ GLM-PCA), Full")
ggplot2::ggsave(filename = "../../../../out/fig/writeup6b/dropseq_mouselung_esvd_full_nb2_glmpca_umap.png",
                plot1, device = "png", width = 5, height = 5, units = "in")

###

set.seed(10)
pred_mat2 <- exp(tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates[,1], esvd_res$b_mat[,1]))
svd_res <- irlba::irlba(pred_mat2, nv = 30)
u_mat <- eSVD2:::.mult_mat_vec(svd_res$u, svd_res$d)
set.seed(10)
tmp <- Seurat::RunUMAP(u_mat)@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap2"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "esvdfactorumap2_", assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap2", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2, initialized w/ GLM-PCA), Full2")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_full2_nb2_glmpca_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


#####################

png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmpca_factor_heatmap.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res$x_mat,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, initialized GLM-PCA) Log-scale")
graphics.off()

max_value <- quantile(mat[mat > 2], probs = 0.999)
nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
mean_mat <- exp(nat_mat)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmpca_scatter_onlymean.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scatterplot_poisson(mat, mean_mat = mean_mat,
                         xlab = "Predicted mean",
                         ylab = "Observed value",
                         xlim = c(0, max_value),
                         main = "Mouse Lung (Dropseq)\neSVD (NB2, initialized va GLM-PCA, only mean)")
graphics.off()

png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmpca_scatter.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                    mean_mat = mean_mat,
                    size_vec = nuisance_vec,
                    log_scale = F,
                    xlim = c(0, max_value),
                    xlab = "Predicted mean",
                    ylab = "Observed value",
                    main = "Mouse Lung (Dropseq)\neSVD (NB2, initialized via GLM-PCA):",
                    include_percentage_in_main = T)
graphics.off()

########################

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates[,1], esvd_res$b_mat[,1])
mean(eSVD2:::.log_prob.neg_binom2(mat,
                                 theta = nat_mat,
                                 s = NULL,
                                 gamma = nuisance_vec))

load("../../../../out/writeup6/writeup6_dropseq_mouselung_glmpca_nb.RData")
mat <- as.matrix(Matrix::t(mat))
nat_mat <- tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)) + tcrossprod(as.matrix(glmpca_res$X), as.matrix(glmpca_res$coefX))
nat_mat <- sweep(nat_mat, 1, glmpca_res$offsets, "+")
mean(eSVD2:::.log_prob.neg_binom2(mat,
                                  theta = nat_mat,
                                  s = NULL,
                                  gamma = rep(glmpca_res$glmpca_family$nb_theta[1], ncol(mat))))


