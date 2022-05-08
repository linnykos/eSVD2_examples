rm(list=ls())
load("../../../out/Writeup11/Writeup11_sns_layer4_esvd.RData")
library(Seurat)
source("../experiment/Writeup11b/multiple_testing.R")
source("../experiment/Writeup11b/plotting.R")
source("../experiment/Writeup11b/reparameterization.R")

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
case_control_variable <- "diagnosis_ASD"
nuisance_vec[is.na(nuisance_vec)] <- 0
res <- compute_posterior(alpha_max = 50,
                         case_control_variable = case_control_variable,
                         esvd_res = esvd_res_full,
                         mat = mat,
                         nuisance_lower_quantile = 0.01,
                         nuisance_vec = nuisance_vec)
######################

metadata <- sns@meta.data
case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
teststat_vec <- compute_test_statistics(case_individuals = case_individuals,
                                        control_individuals = control_individuals,
                                        covariate_individual = "individual",
                                        metadata = metadata,
                                        posterior_mean_mat = res$posterior_mean_mat,
                                        posterior_var_mat = res$posterior_var_mat,
                                        verbose = T)

#############

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "IN-VIP"),]
de_gene_specific <- tmp[,"Gene name"]
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list[-1], function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
de_genes <- de_genes[!de_genes %in% de_gene_specific]
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(colnames(mat) %in% c(hk_genes, cycling_genes))
de_idx <- which(colnames(mat) %in% de_gene_specific)
other_idx <- which(colnames(mat) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]


###########################

set.seed(10)
null_res <- logcondens::logConDens(teststat_vec[hk_idx],
                                   smoothed = T,
                                   print = FALSE,
                                   xs = seq(1.5*min(teststat_vec),
                                            1.5*max(teststat_vec),
                                            length.out = 1000))
dens_val <- null_res$f.smoothed
dens_val <- dens_val * 350/max(dens_val)
max_val <- max(abs(teststat_vec))
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05

teststat_vec <- pmax(pmin(teststat_vec, 30), -30)
max_val <- max(abs(teststat_vec))
png("../../../out/fig/main/sns_layer4_esvd_teststat_histogram_cleaned.png", height = 1200, width = 1200,
    units = "px", res = 300)
break_vec <- seq(-max_val-0.05, max_val+0.05, by = 0.1)
break_vec[1] <- -max_val-0.05; break_vec[length(break_vec)] <- max_val+0.05
hist(teststat_vec, breaks = break_vec,
     xlim = c(-max_val, max_val),
     main = "",
     xlab = "Test statistic", ylab = "Frequency", freq = T)
lines(null_res$xs, dens_val, lwd = 4, lty = 5, col = "white")
lines(null_res$xs, dens_val, lwd = 3, lty = 2, col = 3)
graphics.off()

######################################

set.seed(10)
idx <- which.min(abs(-4-teststat_vec))
case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
case_individuals <- sample(case_individuals, 3)
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
control_individuals <- sample(control_individuals, 3)
cell_list <- lapply(c(case_individuals, control_individuals), function(indiv){
  tmp <- which(sns$individual == indiv)
  sample(tmp, 100)
})
cell_idx <- unlist(cell_list)

esvd_res <- esvd_res_full
offset_var <- setdiff(colnames(esvd_res$covariates), case_control_variable)

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates[,case_control_variable,drop = F],
                       esvd_res$b_mat[,case_control_variable,drop = F])
nat_mat3 <- tcrossprod(esvd_res$covariates,
                       esvd_res$b_mat)
nat_mat_nolib <- nat_mat1 + nat_mat2
nat_mat_full <- nat_mat1 + nat_mat3
mean_mat_nolib <- exp(nat_mat_nolib)
library_mat <- exp(tcrossprod(
  esvd_res$covariates[,offset_var],
  esvd_res$b_mat[,offset_var]
))
y_vec <- mat[cell_idx,idx]/library_mat[cell_idx,idx]
x_vec_prior <- exp(nat_mat_nolib[cell_idx,idx])
x_vec_posterior <- res$posterior_mean_mat[cell_idx,idx]

quantile(y_vec)
quantile(x_vec_prior)
quantile(x_vec_posterior)

cor(y_vec, x_vec_prior)
cor(y_vec, x_vec_posterior)

png("../../../out/fig/main/illustration_sns_layer4_esvd_scatterplot.png",
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(1,1,0,0))
plot(NA,
     xlim = range(c(x_vec_prior, x_vec_posterior)),
     ylim = range(c(y_vec)),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
points(x = x_vec_prior, y = y_vec, lwd = 2, col = rgb(0.7, 0.7, 0.7), cex = 2)
points(x = x_vec_posterior, y = y_vec, pch = 16, col = "white", cex = 2.5)
points(x = x_vec_posterior, y = y_vec, pch = 16, col = rgb(0.2, 0.2, 0.2), cex = 2)
axis(side = 1, labels = F, lwd = 2)
axis(side = 2, labels = F, lwd = 2)
graphics.off()

x_vec_posterior_full <- res$posterior_mean_mat[,idx]
for(i in 1:3){
  png(paste0("../../../out/fig/main/illustration_sns_layer4_esvd_scatterplot_rug", i, ".png"),
      height = 2500, width = 2500,
      units = "px", res = 500)
  par(mar = c(1,1,0,0), bg = NA)
  plot(NA,
       xlim = range(c(y_vec, x_vec_prior, x_vec_posterior)),
       ylim = range(c(y_vec, x_vec_prior, x_vec_posterior)),
       xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "")
  rug(x_vec_posterior_full[cell_list[[i]]], lwd = 2, col = rgb(33, 108, 218, maxColorValue = 255))
  graphics.off()
}
for(i in 4:6){
  png(paste0("../../../out/fig/main/illustration_sns_layer4_esvd_scatterplot_rug", i, ".png"),
      height = 2500, width = 2500,
      units = "px", res = 500)
  par(mar = c(1,1,0,0), bg = NA)
  plot(NA,
       xlim = range(c(y_vec, x_vec_prior, x_vec_posterior)),
       ylim = range(c(y_vec, x_vec_prior, x_vec_posterior)),
       xaxt = "n", yaxt = "n", bty = "n",
       xlab = "", ylab = "")
  rug(x_vec_posterior_full[cell_list[[i]]], lwd = 2, col = rgb(199, 53, 53, maxColorValue = 255))
  graphics.off()
}

##################

xlim <- range(c(y_vec, x_vec_prior, x_vec_posterior))
gaussian_list <- lapply(cell_list, function(cell_idx){
  mean_val <- mean(res$posterior_mean_mat[cell_idx,idx])
  var_val <- mean(res$posterior_var_mat[cell_idx,idx])
  xseq <- seq(xlim[1], xlim[2], length.out = 1000)
  yseq <- stats::dnorm(xseq, mean = mean_val, sd = sqrt(var_val))
  yseq <- yseq - min(yseq)
  yseq <- yseq/3.5*2
  cbind(xseq, yseq)
})

png(paste0("../../../out/fig/main/illustration_sns_layer4_esvd_gaussians.png"),
    height = 2500, width = 2500,
    units = "px", res = 500)
par(mar = c(0,0,0,0), bg = NA)
plot(NA,
     xlim = c(0.25, 1.6),
     ylim = c(0, 7.5),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "", ylab = "")
for(i in 6:1){
  if(i < 4){
    col_val <- rgb(33, 108, 218, maxColorValue = 255)
  } else {
    col_val <- rgb(199, 53, 53, maxColorValue = 255)
  }

  graphics::polygon(x = c(gaussian_list[[i]][,1], rev(gaussian_list[[i]][,1])),
                    y = 0.5+(i-1)+c(gaussian_list[[i]][,2], rep(0, nrow(gaussian_list[[i]]))),
                    col = col_val)
}
graphics.off()
