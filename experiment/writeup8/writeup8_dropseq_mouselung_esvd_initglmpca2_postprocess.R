rm(list=ls())

library(Seurat)
load("../../../../out/writeup8/writeup8_dropseq_mouselung_esvd_initglmpca2.RData")

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                         key = "esvdfactorumap_",
                                                         assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2 w/ Newton, initialized w/ GLM-PCA), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_factor_esvd_newton.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_heatmap_esvd_newton.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res$x_mat,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2 w/ Newton, initialized GLM-PCA)")
graphics.off()

max_value <- quantile(mat2[mat2 > 2], probs = 0.999)
nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
nat_mat <- sweep(nat_mat, 1, esvd_res$offset_vec, "+")
mean_mat <- exp(nat_mat)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_scatter_esvd_newton.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat2,
                            mean_mat = mean_mat,
                            size_vec = esvd_res$nuisance_param_vec,
                            log_scale = F,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Mouse Lung (Dropseq)\neSVD (NB2 w/ Newton, initialized via GLM-PCA):",
                            include_percentage_in_main = T)
graphics.off()

##########

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res2$x_mat))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                         key = "esvdfactorumap_",
                                                         assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\neSVD (NB2 w/ L-BFGS, initialized w/ GLM-PCA), Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_factor_esvd_lbfgs.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_heatmap_esvd_lbfgs.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res2$x_mat,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\neSVD (NB2 w/ L-BFGS, initialized GLM-PCA)")
graphics.off()

max_value <- quantile(mat2[mat2 > 2], probs = 0.999)
nat_mat <- tcrossprod(esvd_res2$x_mat, esvd_res2$y_mat) + tcrossprod(covariates, esvd_res2$b_mat)
nat_mat <- sweep(nat_mat, 1, esvd_res2$offset_vec, "+")
mean_mat <- exp(nat_mat)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_scatter_esvd_lbfgs.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat2,
                            mean_mat = mean_mat,
                            size_vec = esvd_res2$nuisance_param_vec,
                            log_scale = F,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Mouse Lung (Dropseq)\neSVD (NB2 w/ L-BFGS, initialized via GLM-PCA):",
                            include_percentage_in_main = T)
graphics.off()

#############################

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(glmpca_res$factors))@cell.embeddings
rownames(tmp) <- rownames(lung@meta.data)

lung[["glmpcafactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                         key = "glmpcafactorumap_",
                                                         assay = "RNA")

plot1 <- Seurat::DimPlot(lung, reduction = "glmpcafactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Mouse Lung (Dropseq)\nGLM-PCA, Factor")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_factor_glmpca.png",
                plot1, device = "png", width = 6, height = 5, units = "in")


tmp <- eSVD2:::.reparameterize(x_mat = as.matrix(glmpca_res$factors),
                               y_mat = as.matrix(glmpca_res$loadings),
                               equal_covariance = T)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_heatmap_glmpca.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(tmp$x_mat,
                            membership_vec = as.factor(lung@meta.data$celltype),
                            bool_log = T,
                            scaling_power = 1.5,
                            xlab = "Latent factor",
                            ylab = "Cells",
                            main = "Heatmap for Mouse Lung (Dropseq)\nGLM-PCA")
graphics.off()

max_value <- quantile(mat2[mat2 > 2], probs = 0.999)
nat_mat <- tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)) + tcrossprod(as.matrix(glmpca_res$X), as.matrix(glmpca_res$coefX))
nat_mat <- sweep(nat_mat, 1, glmpca_res$offsets, "+")
mean_mat <- exp(nat_mat)
png("../../../../out/fig/writeup8/dropseq_mouselung_esvd_initglmpca2_scatter_glmpca.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
eSVD2:::plot_scatterplot_nb(mat2,
                            mean_mat = mean_mat,
                            size_vec = rep(glmpca_res$glmpca_family$nb_theta, ncol(mat2)),
                            log_scale = F,
                            xlim = c(0, max_value),
                            xlab = "Predicted mean",
                            ylab = "Observed value",
                            main = "Mouse Lung (Dropseq)\nGLM-PCA:",
                            include_percentage_in_main = T)
graphics.off()

####################################

mat3 <- log(eSVD2:::.mult_vec_mat(1/matrixStats::rowSums2(mat2), mat2)*1e4+1)
avg_expression <- matrixStats::colMeans2(mat3)

set.seed(10)
kmean_res <- stats::kmeans(avg_expression, centers = 6)
cbind(kmean_res$centers, kmean_res$size)
gene_grouping <- as.factor(kmean_res$cluster)
umi_vec <- matrixStats::rowSums2(mat2)

# pred_mat <- exp(tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates[,1,drop=F], esvd_res$b_mat[,1,drop=F]))
pred_mat <- exp(tcrossprod(esvd_res$x_mat, esvd_res$y_mat))
zz <- which(apply(pred_mat, 2, function(x){all(abs(x - 1) <= 1e-4)}))

smoothing_original <- smooth_gene_vs_umi(pred_mat,
                                         gene_grouping,
                                         umi_vec = umi_vec)

individual_list <- smoothing_original$individual_list
median_list <- smoothing_original$median_list

xlim <- range(median_list[[1]][,"x"])
# ylim <- range(unlist(lapply(median_list, function(x){x[,-1]})))
png("../../../../out/fig/writeup8/writeup8_dropseq_mouselung_esvd_expression_vs_umi_original.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(NA,
       xlim = xlim,
       ylim = range(median_list[[i]][,-1]),
       xlab = "Total cell UMI count",
       ylab = "Gene UMI count",
       main = paste0("Group ", i, " (", kmean_res$size[i],
                     " genes,\nLog-exp of ",
                     round(kmean_res$centers[i,1], 2), ")"))
  polygon(x = c(median_list[[i]][,"x"], rev(median_list[[i]][,"x"])),
          y = c(median_list[[i]][,"10%"], rev(median_list[[i]][,"90%"])),
          col = grDevices::rgb(0,0,0,0.2),
          density = 30,
          angle = -45)
  lines(x = median_list[[i]][,"x"],
        y = median_list[[i]][,"50%"],
        col = "red",
        lwd = 2)
}
graphics.off()

########

pred_mat <- exp(tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)))

smoothing_original <- smooth_gene_vs_umi(pred_mat,
                                         gene_grouping,
                                         umi_vec = umi_vec)
individual_list <- smoothing_original$individual_list
median_list <- smoothing_original$median_list

xlim <- range(median_list[[1]][,"x"])
# ylim <- range(unlist(lapply(median_list, function(x){x[,-1]})))
png("../../../../out/fig/writeup8/writeup8_dropseq_mouselung_expression_vs_umi_glmpca.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(NA,
       xlim = xlim,
       ylim = range(median_list[[i]][,-1]),
       xlab = "Total cell UMI count",
       ylab = "Gene UMI count",
       main = paste0("Group ", i, " (", kmean_res$size[i],
                     " genes,\nLog-exp of ",
                     round(kmean_res$centers[i,1], 2), ")"))
  polygon(x = c(median_list[[i]][,"x"], rev(median_list[[i]][,"x"])),
          y = c(median_list[[i]][,"10%"], rev(median_list[[i]][,"90%"])),
          col = grDevices::rgb(0,0,0,0.2),
          density = 30,
          angle = -45)
  lines(x = median_list[[i]][,"x"],
        y = median_list[[i]][,"50%"],
        col = "red",
        lwd = 2)
}
graphics.off()
