rm(list=ls())
load("../../../../out/writeup8/writeup8_dropseq_humancortical_esvd_nb_glmgampoi.RData")

set.seed(10)
tmp <- Seurat::RunUMAP(as.matrix(esvd_res$x_mat))@cell.embeddings
rownames(tmp) <- rownames(cortical@meta.data)

cortical[["esvdfactorumap"]] <- Seurat::CreateDimReducObject(embedding = tmp, key = "UMAP", assay = "RNA")

options(ggrepel.max.overlaps = Inf)
plot1 <- Seurat::DimPlot(cortical, reduction = "esvdfactorumap", group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5) + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + ggplot2::ggtitle("Human cortical (Dropseq)\neSVD (via GLMGamPoi), Factor, NB")
ggplot2::ggsave(filename = "../../../../out/fig/writeup8/dropseq_humancortical_esvd_factor_nb_glmgampoi_umap.png",
                plot1, device = "png", width = 7, height = 5, units = "in")

#########################

png("../../../../out/fig/writeup8/dropseq_humancortical_esvd_factor_nb_glmgampoi_factor_heatmap.png",
    height = 2000, width = 2000, units = "px", res = 300)
eSVD2:::plot_scores_heatmap(esvd_res$x_mat,
                           membership_vec = as.factor(cortical@meta.data$celltype),
                           bool_log = T,
                           scaling_power = 1.5,
                           xlab = "Latent factor",
                           ylab = "Cells",
                           main = "Heatmap for Human cortical (Dropseq)\neSVD (NB, initialized glmGamPoi) Log-scale")
graphics.off()

max_value <- quantile(mat[mat > 2], probs = 0.999)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_humancortical_esvd_nb_glmgampoi_scatter_init.png",
    height = 2000, width = 2000, units = "px", res = 300)
plot_scatterplot_poisson(mat, mean_mat = t(glmGamPoi_res$Mu),
                 xlab = "Predicted mean",
                 ylab = "Observed value",
                 xlim = c(0, max_value),
                 main = "Human cortical (Dropseq)\neSVD (Initial mean, via glmGamPoi)")
graphics.off()

nat_mat <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat) + tcrossprod(covariates, esvd_res$b_mat)
prob_mat <- exp(nat_mat)
mean_mat <-  eSVD2:::.mult_mat_vec(prob_mat/(1-prob_mat), nuisance_vec)
set.seed(10)
png("../../../../out/fig/writeup8/dropseq_humancortical_esvd_nb_glmgampoi_scatter_onlymean.png",
    height = 2000, width = 2000, units = "px", res = 300)
plot_scatterplot_poisson(mat, mean_mat = mean_mat,
                         xlab = "Predicted mean",
                         ylab = "Observed value",
                         xlim = c(0, max_value),
                         main = "Human cortical (Dropseq)\neSVD (NB initialized va glmGamPoi, only mean)")
graphics.off()

png("../../../../out/fig/writeup8/dropseq_humancortical_esvd_nb_glmgampoi_scatter.png",
    height = 2000, width = 2000, units = "px", res = 300)
set.seed(10)
plot_scatterplot_nb(mat,
                    prob_mat = prob_mat,
                    size_vec = nuisance_vec,
                    log_scale = F,
                    xlim = c(0, max_value),
                    xlab = "Predicted mean",
                    ylab = "Observed value",
                    main = "Human cortical (Dropseq)\neSVD (NB, via glmGamPoi):",
                    include_percentage_in_main = T)
graphics.off()

########################################

x_vec <- log(matrixStats::rowSums2(mat)+1)
png("../../../../out/fig/writeup8/dropseq_humancortical_librarysize.png",
    height = 1500, width = 3000, units = "px", res = 300)
par(mfrow = c(1,2))
hist(x_vec, breaks = 50, xlab = "Log-library size", col = "gray",
     main = "Human cortical (Dropseq)")

plot(x = x_vec, y = covariates[,"Log-UMI"], pch = 16, asp = T,
     xlab = "Log-library size",
     ylab = "Scran estimate")
graphics.off()

mat2 <- log(eSVD2:::.mult_vec_mat(1/matrixStats::rowSums2(mat), mat)*1e4+1)
avg_expression <- matrixStats::colMeans2(mat2)
perc_expression <- apply(mat2, 2, function(x){length(which(x > 0))/length(x)})
png("../../../../out/fig/writeup8/dropseq_humancortical_expression.png",
    height = 1500, width = 3000, units = "px", res = 300)
par(mfrow = c(1,2))
hist(avg_expression, breaks = 50, xlab = "Mean expression (Log-normalized)", col = "gray",
     main = "Human cortical (Dropseq)")

plot(x = avg_expression, y = perc_expression,
     pch = 16,
     xlab = "Mean expression (Log-normalized)",
     ylab = "Percent expression")
graphics.off()

set.seed(10)
kmean_res <- stats::kmeans(avg_expression, centers = 6)
cbind(kmean_res$centers, kmean_res$size)
gene_grouping <- as.factor(kmean_res$cluster)

smoothing_original <- smooth_gene_vs_umi(mat,
                                         gene_grouping)
individual_list <- smoothing_original$individual_list
median_list <- smoothing_original$median_list

xlim <- range(median_list[[i]][,"x"])
ylim <- range(unlist(lapply(median_list, function(x){x[,-1]})))
png("../../../../out/fig/writeup8/dropseq_humancortical_expression_vs_umi_original.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(NA,
       xlim = xlim,
       ylim = ylim,
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
