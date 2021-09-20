rm(list=ls())

library(Seurat)
load("../../../../out/writeup8b/writeup8b_dropseq_mouselung_esvd.RData")

quantile(esvd_res$nuisance_param_vec)
quantile(esvd_res$b_mat[,1])
quantile(esvd_res$b_mat[,2])

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                         key = "esvdfactorumap_",
                                                         assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8b/dropseq_mouselung_esvd_factor1.png",
                plot1, device = "png", width = 6, height = 5, units = "in")

#########

png("../../../../out/fig/writeup8b/dropseq_mouselung_esvd_heatmap1.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res$x_mat,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2)")
graphics.off()

max_value <- quantile(mat[mat > 2], probs = 0.999)
nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
nat_mat <- sweep(nat_mat, 1, esvd_res$offset_vec, "+")
mean_mat <- exp(nat_mat)
png("../../../../out/fig/writeup8b/dropseq_mouselung_esvd_scatter1.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat,
                            mean_mat = mean_mat,
                            size_vec = esvd_res$nuisance_param_vec,
                            log_scale = F,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Mouse Lung (Dropseq)\neSVD (NB2):",
                            include_percentage_in_main = T)
graphics.off()

#####################

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
mean_mat <- exp(nat_mat)
quantile(mean_mat[mat != 0])

