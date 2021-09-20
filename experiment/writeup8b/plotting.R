plotting_func <- function(seurat_obj,
                          esvd_res,
                          covariates,
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
                              membership_vec = as.factor(seurat_obj@meta.data$celltype),
                              bool_log = T,
                              scaling_power = 1.5,
                              xlab = "Latent factor",
                              ylab = "Cells",
                              main = heatmap_title)
  graphics.off()

  nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
  nat_mat <- sweep(nat_mat, 1, esvd_res$offset_vec, "+")
  mean_mat <- exp(nat_mat)
  png(scatter_filename,
      height = 2000, width = 2000, units = "px", res = 300)
  set.seed(10)
  eSVD2:::plot_scatterplot_nb(mat,
                              mean_mat = mean_mat,
                              size_vec = esvd_res$nuisance_param_vec,
                              log_scale = F,
                              xlim = c(0, max_value),
                              xlab = "Predicted mean",
                              ylab = "Observed value",
                              main = scatter_title,
                              include_percentage_in_main = T)
  graphics.off()

  invisible()
}
