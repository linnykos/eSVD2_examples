plotting_func <- function(seurat_obj,
                          mat,
                          esvd_res,
                          group.by,
                          max_value,
                          factor_title,
                          factor_filename,
                          heatmap_title,
                          heatmap_filename,
                          scatter_title,
                          scatter_filename){
  set.seed(10)
  tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
  rownames(tmp) <- rownames(seurat_obj@meta.data)

  seurat_obj[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                                 key = "esvdfactorumap_",
                                                                 assay = "RNA")

  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "esvdfactorumap", group.by = group.by, label = TRUE,
                           repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  plot1 <- plot1 + ggplot2::ggtitle(factor_title)
  ggplot2::ggsave(filename = factor_filename,
                  plot1, device = "png", width = 6, height = 5, units = "in")

  #########

  png(heatmap_filename,
      height = 2000, width = 2000, units = "px", res = 300)
  eSVD2:::plot_scores_heatmap(esvd_res$x_mat,
                              membership_vec = as.factor(seurat_obj@meta.data[,group.by]),
                              bool_log = T,
                              scaling_power = 1.5,
                              xlab = "Latent factor",
                              ylab = "Cells",
                              main = heatmap_title)
  graphics.off()

  nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(esvd_res$covariates, esvd_res$b_mat)
  nat_mat <- sweep(nat_mat, 1, esvd_res$offset_vec, "+")
  mean_mat <- exp(nat_mat)
  png(scatter_filename,
      height = 2000, width = 2000, units = "px", res = 300)
  set.seed(10)
  eSVD2:::plot_scatterplot_nb(mat,
                              mean_mat = mean_mat,
                              size_vec = esvd_res$nuisance_param_vec,
                              quantile_shoulder = 0.5,
                              xlim = c(0, max_value),
                              xlab = "Predicted mean",
                              ylab = "Observed value",
                              main = scatter_title,
                              include_percentage_in_main = T)
  graphics.off()

  invisible()
}

############################

plotting_func_glmpca <- function(seurat_obj,
                                 mat,
                                 glmpca_res,
                                 group.by,
                                 max_value,
                                 factor_title,
                                 factor_filename,
                                 heatmap_title,
                                 heatmap_filename,
                                 scatter_title,
                                 scatter_filename){

  reparam_res <- eSVD2:::.reparameterize(as.matrix(glmpca_res$factors),
                                         as.matrix(glmpca_res$loadings),
                                         equal_covariance = T)

  set.seed(10)
  tmp <- Seurat::RunUMAP(reparam_res$x_mat)@cell.embeddings
  rownames(tmp) <- rownames(seurat_obj@meta.data)

  seurat_obj[["glmpcafactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp,
                                                                   key = "glmpcafactorumap_",
                                                                   assay = "RNA")

  plot1 <- Seurat::DimPlot(seurat_obj, reduction = "glmpcafactorumap", group.by = group.by, label = TRUE,
                           repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  plot1 <- plot1 + ggplot2::ggtitle(factor_title)
  ggplot2::ggsave(filename = factor_filename,
                  plot1, device = "png", width = 6, height = 5, units = "in")

  #########

  png(heatmap_filename,
      height = 2000, width = 2000, units = "px", res = 300)
  eSVD2:::plot_scores_heatmap(reparam_res$x_mat,
                              membership_vec = as.factor(seurat_obj@meta.data[,group.by]),
                              bool_log = T,
                              scaling_power = 1.5,
                              xlab = "Latent factor",
                              ylab = "Cells",
                              main = heatmap_title)
  graphics.off()

  # nat_mat <- tcrossprod(as.matrix(glmpca_res$factors), as.matrix(glmpca_res$loadings)) +
  #   tcrossprod(as.matrix(glmpca_res$X), as.matrix(glmpca_res$coefX))
  # nat_mat <- sweep(nat_mat, 1, glmpca_res$offsets, "+")
  # mean_mat <- exp(nat_mat)
  # png(scatter_filename,
  #     height = 2000, width = 2000, units = "px", res = 300)
  # set.seed(10)
  # eSVD2:::plot_scatterplot_nb(mat,
  #                             mean_mat = mean_mat,
  #                             size_vec = rep(glmpca_res$glmpca_family$nb_theta, ncol(mean_mat)),
  #                             log_scale = F,
  #                             quantile_shoulder = 0.5,
  #                             xlim = c(0, max_value),
  #                             xlab = "Predicted mean",
  #                             ylab = "Observed value",
  #                             main = scatter_title,
  #                             include_percentage_in_main = T)
  # graphics.off()

  invisible()
}

