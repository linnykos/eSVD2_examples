rm(list=ls())

library(Seurat)
load("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_nb2_glmgampoi.RData")
quantile(nuisance_vec)

###################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2, initialized w/ GLMGamPoi), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_factor_nb2_glmgampoi_umap.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

#####################

png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmgampoi_factor_heatmap.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res$x_mat,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2, initialized glmGamPoi) Log-scale")
graphics.off()

max_value <- quantile(mat[mat > 2], probs = 0.999)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmgampoi_scatter_init.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scatterplot_poisson(mat, mean_mat = t(glmGamPoi_res$Mu),
                         xlab = "Predicted mean",
                         ylab = "Observed value",
                         xlim = c(0, max_value),
                         main = "Mouse Lung (Dropseq)\neSVD (Initial mean, via glmGamPoi)")
graphics.off()

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
mean_mat <- exp(nat_mat)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmgampoi_scatter_onlymean.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scatterplot_poisson(mat, mean_mat = mean_mat,
                         xlab = "Predicted mean",
                         ylab = "Observed value",
                         xlim = c(0, max_value),
                         main = "Mouse Lung (Dropseq)\neSVD (NB2, initialized va glmGamPoi, only mean)")
graphics.off()

png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_nb2_glmgampoi_scatter.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                    mean_mat = mean_mat,
                    size_vec = nuisance_vec,
                    log_scale = F,
                    xlim = c(0, max_value),
                    xlab = "Predicted mean",
                    ylab = "Observed value",
                    main = "Mouse Lung (Dropseq)\neSVD (NB2, initialized via glmGamPoi):",
                    include_percentage_in_main = T)
graphics.off()
