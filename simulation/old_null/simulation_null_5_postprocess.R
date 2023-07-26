rm(list=ls())
load("../eSVD2_examples/simulation/simulation_null_5_esvd.RData")

gene_plot_idx <- c(which(y_block_assignment == 1), which(y_block_assignment == 2), which(y_block_assignment == 3))

# let's see what the heatmap looks like without posterior
denoised_mat <- tcrossprod(eSVD_obj$fit_Second$x_mat,
                           eSVD_obj$fit_Second$y_mat) +
  tcrossprod(eSVD_obj$covariates[,"CC"], eSVD_obj$fit_Second$z_mat[,"CC"])
image(denoised_mat[,gene_plot_idx])
image(cor(denoised_mat[,gene_plot_idx]))
cor_mat <- cor(denoised_mat)
image(cor_mat[which(y_block_assignment == 1), which(y_block_assignment == 1)])
image(cor_mat[which(y_block_assignment == 2), which(y_block_assignment == 2)])
image(cor_mat[which(y_block_assignment == 3), which(y_block_assignment == 3)])
