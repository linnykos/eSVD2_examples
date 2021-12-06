rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson4.RData")
hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]
mat[which(mat == min(mat))] <- 0
#########

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates, esvd_res_full$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

# save(mat, mean_mat, file = "../../out/writeup8e/ns_layer23_example.RData")

library_idx <- which(colnames(esvd_res_full$covariates) == "Log_UMI")
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,-library_idx], esvd_res_full$b_mat[,-library_idx])
nat_mat_nolib <- nat_mat1 + nat_mat2
mean_mat_nolib <- exp(nat_mat_nolib)
mean_mat_nolib <- pmin(mean_mat_nolib, 1e4)

zero_prop <- apply(mat, 2, function(x){length(which(x == 0))/length(x)})
mode_bool <- apply(mat, 2, function(x){
  if(all(x != 0)) return(FALSE)
  tab <- table(x)
  zero_idx <- which(names(tab) == "0")
  tab[zero_idx] > max(tab[-zero_idx])
})

nuisance_param_vec <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')

  val1 <- MASS::theta.ml(y = mat[,j], mu = mean_mat[,j])
  # if(colnames(mat)[j] %in% de_genes){
  #   return(val1)
  # }
  val2 <- MASS::theta.mm(y = mat[,j], mu = mean_mat[,j], dfr = nrow(mat)-1)
  val3 <- glmGamPoi::overdispersion_mle(y = mat[,j], mean = mean_mat[,j])$estimate

  # min(c(val1, val2, val3), na.rm = T)

  vec <- c(val1, val2, val3)
  vec <- vec[!is.na(vec)]
  if(length(vec) == 1) return(vec[1])
  vec <- pmax(pmin(vec, 1e5), 0.1)

  # honestly -- this doesn't make too much sense. The mode can still be zero even if
  # the nuisance parameter isn't less than 1. (It's the converse that's true.)
  # if(mode_bool[j]){
    vec <- c(vec, c(0.1, 0.5, 1))

    obs_prob <- length(which(mat[,j] == 0))/nrow(mat)
    target_prob_vec <- sapply(vec, function(val){
      mean((1+mean_mat[,j]/val)^(-val))
    })
    return(vec[which.min(abs(target_prob_vec - obs_prob))])
  # } else {
  #   cat('*')
  #   target_quantile_vec <- sapply(vec, function(val){
  #     if(val < 1){
  #       mode_val <- rep(0, nrow(mat))
  #     } else {
  #       mode_val <- mean_mat[,j]*(val-1)/val
  #     }
  #
  #     mean(sapply(1:nrow(mat), function(i){
  #       if(mode_val[i] == 0){
  #         lower_val <- 0
  #         upper_val <- stats::qnbinom(0.75, size = val, mu = mean_mat[i,j])
  #       } else {
  #         quantile_val <- stats::pnbinom(mode_val[i], size = val, mu = mean_mat[i,j])
  #         lower_val <- stats::qnbinom(quantile_val*.25, size = val, mu = mean_mat[i,j])
  #         upper_val <- stats::qnbinom(quantile_val + (1-quantile_val)*.75, size = val, mu = mean_mat[i,j])
  #       }
  #
  #       lower_val <= mat[i,j] & mat[i,j] <= upper_val
  #     }))
  #   })
  #
  #   if(any(target_quantile_vec > 0.75)){
  #     vec <- vec[target_quantile_vec > 0.75]
  #     target_quantile_vec <- target_quantile_vec[target_quantile_vec > 0.75]
  #   }
  #
  #   vec[which.min(abs(target_quantile_vec - 0.75))]
  # }
})
quantile(nuisance_param_vec)
# length(intersect(which(zero_prop <= 0.2), which(nuisance_param_vec == 1e5)))

library_mat <- sapply(1:ncol(mat), function(j){
  exp(esvd_res_full$covariates[,"Log_UMI",drop = F]*esvd_res_full$b_mat[j,"Log_UMI"])
})
library_mat <- pmin(library_mat, 50000)
AplusR <- sweep(mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat_nolib, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + library_mat
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2
tmp <- posterior_mean_mat/sqrt(posterior_var_mat)
quantile(tmp)
tmp <- posterior_mean_mat/mean_mat_nolib
round(quantile(tmp),4)

###############

png("../../out/fig/writeup8e/tmp.png",
    height = 2000, width = 1500, units = "px", res = 300)
image(t(esvd_res_full$x_mat))
graphics.off()

png("../../out/fig/writeup8e/tmp2.png",
    height = 2000, width = 1500, units = "px", res = 300)
image(t(esvd_res_full$y_mat))
graphics.off()

######################

nat_mat1 <- tcrossprod(esvd_res_full$x_mat, esvd_res_full$y_mat)
nat_mat2 <- tcrossprod(esvd_res_full$covariates[,c("Intercept", "diagnosis_ASD")], esvd_res_full$b_mat[,c("Intercept", "diagnosis_ASD")])
nat_mat_clean <- nat_mat1 + nat_mat2
mean_mat_clean <- exp(nat_mat_clean)

ratio_mat <- mean_mat_clean/mean_mat_nolib
quantile(ratio_mat)
posterior_mean_mat2 <- posterior_mean_mat * ratio_mat

case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
case_idx <- which(esvd_res_full$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res_full$covariates[,"diagnosis_ASD"] == 0)

individual_stats <- lapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  r_val <- nuisance_param_vec[j]

  # next find the cells, then compute one gaussian per individual
  case_gaussians <- sapply(case_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
    var_val <- mean(posterior_var_mat[cell_idx,j])
    c(mean_val = mean_val, var_val = var_val)
  })

  list(case_gaussians = case_gaussians,
       control_gaussians = control_gaussians)
})

group_stats <- lapply(1:length(individual_stats), function(j){
  case_gaussians <- individual_stats[[j]]$case_gaussians
  control_gaussians <- individual_stats[[j]]$control_gaussians

  case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                        var_val = mean(case_gaussians[2,]),
                        n = ncol(case_gaussians))
  control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                           var_val = mean(control_gaussians[2,]),
                           n = ncol(control_gaussians))

  list(case_gaussian = case_gaussian,
       control_gaussian = control_gaussian)
})

p_val_vec <- sapply(1:length(group_stats), function(j){
  case_gaussian <- group_stats[[j]]$case_gaussian
  control_gaussian <- group_stats[[j]]$control_gaussian

  n1 <- control_gaussian$n; n2 <- case_gaussian$n
  mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
  cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val

  combined_cov <- cov1/n1 + cov2/n2
  test_stat <- abs(mean1 - mean2)/sqrt(combined_cov)
  df <- (cov1/n1 + cov2/n2)^2/((cov1/n1)^2/(n1-1) + (cov2/n2)^2/(n2-1))

  stats::pt(test_stat, df = df, lower.tail = F, log.p = T)/log(10) + log10(2)
})
hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)

x_vec <- sapply(1:ncol(mat), function(j){
  # log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))
  log2(mean(mat[case_idx,j])) - log2(mean(mat[control_idx,j]))
})

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(p_val_vec))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

### let's draw it nicer
png("../../out/fig/writeup8e/sns_layer23_volcano.png", height = 1200, width = 1200,
    units = "px", res = 300)
plot(NA, xlim = c(-2,2), ylim = range(0, 20), bty = "n",
     main = "Volcano plot for Layer 2/3",
     xlab = "Log2 fold change (i.e., Log2 mean difference)", ylab = "-Log10(P value)")
for(x in seq(-2,2,by=0.5)){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
lines(rep(0,2), c(-1e5,1e5), col = "gray")
for(y in seq(0,max(-p_val_vec),by=2)){
  lines(c(-1e5,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(x = x_vec[-unique(c(hk_idx,de_idx))],
       y = -p_val_vec[-unique(c(hk_idx,de_idx))],
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = "white", cex = 1.5)
points(x = x_vec[shuf_idx],
       y = -p_val_vec[shuf_idx],
       pch = 16, col = col_vec[shuf_idx])
legend("topright", c("Published DE gene", "Housekeeping gene", "Other"),
       fill = c(2,3,rgb(0.5,0.5,0.5)), cex = 0.6)
graphics.off()

val <- -2
length(which(p_val_vec[-de_idx] <= val))/(length(p_val_vec) - length(de_idx))
length(which(p_val_vec[de_idx] <= val))/length(de_idx)
length(which(p_val_vec <= val))

zz <- 10^p_val_vec/2
zz[x_vec > 0] <- .5 + (.5-zz[x_vec > 0])
zz <- pmax(pmin(stats::qnorm(zz), 10), -10)
png("../../out/fig/writeup8e/sns_layer23_zscore_histogram.png", height = 1200, width = 1200,
    units = "px", res = 300)
hist(zz, breaks = seq(-10.05, 10.05, by = 0.1), xlim = c(-5,5),
     main = "Histogram of two-sided Z-scores",
     xlab = "Z-score", ylab = "Frequency")
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(zz[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Housekeeping gene"),
       fill = c(2,3), cex = 0.6)
graphics.off()

png("../../out/fig/writeup8e/sns_layer23_pval_vs_nuisance.png",
    height = 900, width = 2500,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
plot(nuisance_param_vec, p_val_vec, col = rgb(0.5,0.5,0.5,0.5), pch = 16,
     xlab = "Nuisance parameter", ylab = "-Log10(P value)",
     main = "P value vs. nuisance,\nper gene")
plot(nuisance_param_vec, p_val_vec, xlim = c(0, 100), ylim = c(-10,0),
     col = rgb(0.5,0.5,0.5,0.3), pch = 16,
     xlab = "Nuisance parameter", ylab = "-Log10(P value)",
     main = "Zoom in 1")
plot(nuisance_param_vec, p_val_vec, xlim = c(0.1, 1), ylim = c(-2,0),
     col = rgb(0.5,0.5,0.5,0.1), pch = 16,
     xlab = "Nuisance parameter", ylab = "-Log10(P value)",
     main = "Zoom in 2")
graphics.off()


png("../../out/fig/writeup8e/sns_layer23_nuisance_vs_zeroprop.png",
    height = 900, width = 2500,
    units = "px", res = 300)
par(mfrow = c(1,3))
plot(zero_prop, nuisance_param_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.5),
     xlab = "0-percentage", ylab = "Nuisance parameter",
     main = "Nuisance vs.Percentage\nof 0's, per gene")
plot(zero_prop, nuisance_param_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.3),
     ylim = c(0,100),
     xlab = "0-percentage", ylab = "Nuisance parameter",
     main = "Zoom in 1")
plot(zero_prop, nuisance_param_vec, pch = 16, col = rgb(0.5,0.5,0.5,0.1),
     ylim = c(0,1),
     xlab = "0-percentage", ylab = "Nuisance parameter",
     main = "Zoom in 2")
graphics.off()

hist(10^p_val_vec)
rug(10^(p_val_vec[shuf_idx]), col = col_vec[shuf_idx], lwd = 2)

#############################

png("../../out/fig/writeup8e/sns_layer23_histograms.png",
    height = 1200, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
hist(esvd_res_full$b_mat[,"Intercept"],
     xlab = "Intercept", ylab = "Frequency", main = "Intercept")
hist(esvd_res_full$b_mat[,"Log_UMI"],
     xlab = "Log_UMI", ylab = "Frequency", main = "Log_UMI")
hist(esvd_res_full$b_mat[,"diagnosis_ASD"],
     xlab = "diagnosis_ASD", ylab = "Frequency", main = "diagnosis_ASD")
for(i in shuf_idx){
  rug(esvd_res_full$b_mat[i,"diagnosis_ASD"], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Housekeeping gene"),
       fill = c(2,3))
graphics.off()


#############################

# plot some example genes -- we'll plot the following things:
# 1) observed vs predicted values
# 2) observed vs posterior mean
# 3) predicted value (no lib) vs cleaned mean
# 4) standard deviation vs cleaned mean

case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
all_individuals <- c(case_individuals, control_individuals)
col_vec_individuals <- c(rep(2, length(case_individuals)), rep(3, length(case_individuals)))
shuf_idx <- sample(1:length(col_vec_individuals))
all_individuals <- all_individuals[shuf_idx]; col_vec_individuals <- col_vec_individuals[shuf_idx]
gene_indices <- c(which(colnames(mat) == "SAT2"),
                  which(colnames(mat) == "DEXI"),
                  which.min(p_val_vec),
                  which.max(abs(x_vec)),
                  which.min(zero_prop),
                  which.min(abs(zero_prop - 0.9)),
                  which.min(abs(zero_prop - 0.2)))
for(idx in gene_indices){
  png(paste0("../../out/fig/writeup8e/sns_layer23_gene-", colnames(mat)[idx], ".png"),
      height = 2500, width = 2500,
      units = "px", res = 300)
  par(mfrow = c(2,2), mar = c(4,4,4,0.5))

  tmp <- cbind(mean_mat[,idx], mat[,idx])
  xlim <- range(tmp)
  plot(tmp[,1], tmp[,2], asp = T,
       xlim = xlim, ylim = xlim,
       xlab = "Predicted mean (full)", ylab = "Observed value",
       main = paste0("eSVD fit for ", colnames(mat)[idx], "\n0-percentage: ",
                     round(zero_prop[idx], 2), ", Nuisance: ", round(nuisance_param_vec[idx], 2)),
       pch = 16, col = rgb(0.5,0.5,0.5,0.1))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = "white", pch = 16, cex = 3)
  }
  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = col_vec_individuals[i], pch = 16, cex = 2)
  }

  ###

  tmp <- cbind(posterior_mean_mat[,idx], mat[,idx]/library_mat[,idx])
  xlim <- range(tmp)
  plot(tmp[,1], tmp[,2], asp = T,
       xlim = xlim, ylim = xlim,
       xlab = "Posterior mean (full)", ylab = "Observed relative expression",
       main = paste0("Posterior mean fit\nSize coefficient: ", round(esvd_res_full$b_mat[idx,"Log_UMI"], 2)),
       pch = 16, col = rgb(0.5,0.5,0.5,0.1))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = "white", pch = 16, cex = 3)
  }
  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = col_vec_individuals[i], pch = 16, cex = 2)
  }

  ###

  tmp <- cbind(mean_mat_clean[,idx], mean_mat_nolib[,idx])
  xlim <- range(tmp)
  plot(tmp[,1], tmp[,2], asp = T,
       xlim = xlim, ylim = xlim,
       xlab = "Predicted mean (Clean)", ylab = "Predicted mean (No lib)",
       main = "Posterior mean with and without covariates",
       pch = 16, col = rgb(0.5,0.5,0.5,0.1))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = "white", pch = 16, cex = 3)
  }
  for(i in 1:length(col_vec_individuals)){
    indiv <- all_individuals[i]
    cell_idx <- which(metadata$individual == indiv)
    vec <- colMeans(tmp[cell_idx,])
    points(vec[1], vec[2], col = col_vec_individuals[i], pch = 16, cex = 2)
  }

  ###

  xlim <- range(posterior_mean_mat2[,idx], sqrt(posterior_var_mat[,idx]))
  plot(posterior_mean_mat2[,idx], sqrt(posterior_var_mat[,idx]),
       xlab = "Posterior mean (Clean)", ylab = "Posterior standard deviation",
       pch = 16, col = rgb(0.5,0.5,0.5,0.1),
       xlim = c(0, max(posterior_mean_mat2[,idx])),
       ylim = c(0, max(sqrt(posterior_var_mat[,idx]))),
       main = paste0("Posterior standard deviation vs. mean\n-Log10(P value): ",
                     round(-p_val_vec[idx], 2)))
  lines(c(0, 2*xlim[2]), c(0, 2*xlim[2]), col = 2, lty = 2, lwd = 2)

  tmp <- cbind(individual_stats[[idx]]$case_gaussians, individual_stats[[idx]]$control_gaussians)
  col_vec <- c(rep(2, ncol(individual_stats[[idx]]$case_gaussians)),
               rep(3, ncol(individual_stats[[idx]]$control_gaussians)))
  shuf_idx <- sample(1:length(col_vec))
  tmp <- tmp[,shuf_idx]; col_vec <- col_vec[shuf_idx]
  points(tmp[1,], sqrt(tmp[2,]), pch = 16, col = "white", cex = 3)
  points(tmp[1,], sqrt(tmp[2,]), pch = 16, col = col_vec, cex = 2)

  graphics.off()
}

##########################################

sns <- Seurat::CreateSeuratObject(counts = t(mat), meta.data = metadata)

# plot some umaps
set.seed(10)
umap_res <- Seurat::RunUMAP(esvd_res_full$x_mat)
sns[["tmp"]] <- Seurat::CreateDimReducObject(embeddings = umap_res@cell.embeddings)

covariates <- c("diagnosis", "sex", "individual", "region", "Capbatch", "Seqbatch")
for(covariate in covariates){
  plot1 <- Seurat::DimPlot(sns, reduction = "tmp",
                           group.by = covariate, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via eSVD: ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_esvd2_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

covariates <- c("nCount_RNA", "age", "post.mortem.hours")
for(covariate in covariates){
  plot1 <- Seurat::FeaturePlot(sns,
                               features = covariate,
                               reduction = "tmp")
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via eSVD: ", covariate))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_esvd2_umap_", covariate, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

####################

set.seed(10)
zz <- apply(esvd_res_full$x_mat, 2, eSVD2:::.l2norm)
tmp_vec1 <- esvd_res_full$covariates[,"diagnosis_ASD"]
tmp_vec2 <- esvd_res_full$b_mat[,"diagnosis_ASD"]
zz1 <- eSVD2:::.l2norm(tmp_vec1)
zz2 <- eSVD2:::.l2norm(tmp_vec2)
tmp_vec1 <- tmp_vec1*sqrt(zz2/zz1)
tmp_vec2 <- tmp_vec2*sqrt(zz1/zz2)
tmp_mat <- cbind(esvd_res_full$x_mat, tmp_vec1)
umap_res <- Seurat::RunUMAP(tmp_mat)
sns[["tmp"]] <- Seurat::CreateDimReducObject(embeddings = umap_res@cell.embeddings)
plot1 <- Seurat::DimPlot(sns, reduction = "tmp",
                         group.by = "diagnosis", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("SNS (Layer 2/3) via eSVD:\ndiagnosis_ASD, with covariate"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/fig/writeup8e/sns_layer23_esvd2_umap_diagnosis_ASD_withcovariate.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##############################
