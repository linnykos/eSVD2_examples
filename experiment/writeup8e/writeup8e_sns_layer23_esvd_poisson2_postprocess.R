rm(list=ls())
load("../../out/writeup8e/writeup8e_sns_layer23_esvd_poisson2.RData")

########

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates, esvd_res$b_mat)
nat_mat <- nat_mat1 + nat_mat2
mean_mat <- exp(nat_mat)

nat_mat_wlibrary <- sweep(nat_mat, MARGIN = 1, STATS = esvd_res$offset_vec, FUN = "+")
mean_mat_wlibrary <- exp(nat_mat_wlibrary)
p <- ncol(mat)
angle_vec <- sapply(1:p, function(j){
  tmp <- cbind(mat[,j], mean_mat[,j])
  eSVD2:::.compute_principal_angle(tmp)
})
quantile(angle_vec)

nuisance_param_vec <- sapply(1:p, function(j){
  MASS::theta.mm(y = mat[,j], mu = mean_mat_wlibrary[,j], dfr = nrow(mat)-1)
})
nuisance_param_vec <- pmin(nuisance_param_vec, 1e4)

plot(angle_vec, nuisance_param_vec)

library_vec <- matrixStats::rowSums2(mat)
A <- sweep(mat, MARGIN = 1, STATS = library_vec, FUN = "/")
AplusR <- sweep(A, MARGIN = 2, STATS = nuisance_param_vec, FUN = "+")
RoverMu <- 1/sweep(mean_mat, MARGIN = 2, STATS = nuisance_param_vec, FUN = "/")
RoverMuplusS <- RoverMu + 1
posterior_mean_mat <- AplusR/RoverMuplusS
posterior_var_mat <- AplusR/RoverMuplusS^2
tmp <- posterior_mean_mat/posterior_var_mat
quantile(tmp)

##########################

nat_mat1 <- tcrossprod(esvd_res$x_mat, esvd_res$y_mat)
nat_mat2 <- tcrossprod(esvd_res$covariates[,c("Intercept", "diagnosis_ASD")], esvd_res$b_mat[,c("Intercept", "diagnosis_ASD")])
nat_mat_clean <- nat_mat1 + nat_mat2
mean_mat_clean <- exp(nat_mat_clean)

ratio_mat <- mean_mat_clean/mean_mat
posterior_mean_mat2 <- posterior_mean_mat * ratio_mat

# gene_idx <- which(colnames(mat) == "SAT2")
gene_idx <- which(colnames(mat) == "MX2")
case_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 0)

vioplot::vioplot(posterior_mean_mat2[case_idx, gene_idx],
                 posterior_mean_mat2[control_idx, gene_idx])

vec <- jitter(mat[,gene_idx])
length(which(mat[,gene_idx] == 0))/length(vec)
posterior_vec <- posterior_mean_mat2[,gene_idx]
shuffle_idx <- sample(length(vec))
col_vec <- rep(2, length(vec))
col_vec[control_idx] <- 1
plot(posterior_vec[shuffle_idx], vec[shuffle_idx], col = col_vec[shuffle_idx])

#####################################

# let's try my janky idea now
# first find the individuals
case_individuals <- unique(metadata[which(metadata$diagnosis == "ASD"),"individual"])
control_individuals <- unique(metadata[which(metadata$diagnosis == "Control"),"individual"])
case_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 0)

summary_stats <- sapply(1:ncol(mat), function(j){
  if(j %% floor(ncol(mat)/10) == 0) cat('*')
  r_val <- nuisance_param_vec[j]

  # next find the cells, then compute one gaussian per individual
  case_gaussians <- sapply(case_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
    var_val <- sum(posterior_var_mat[cell_idx,j])/length(cell_idx)^2
    c(mean_val = mean_val, var_val = var_val)
  })

  control_gaussians <- sapply(control_individuals, function(indiv){
    cell_names <- rownames(metadata)[which(metadata$individual == indiv)]
    cell_idx <- which(rownames(mat) %in% cell_names)

    mean_val <- mean(posterior_mean_mat2[cell_idx,j])
    var_val <- sum(posterior_var_mat[cell_idx,j])/length(cell_idx)^2
    c(mean_val = mean_val, var_val = var_val)
  })

  case_gaussian <- list(mean_val = mean(case_gaussians[1,]),
                        var_val = sum(case_gaussians[2,])/ncol(case_gaussians)^2,
                        n = ncol(case_gaussians))
  control_gaussian <- list(mean_val = mean(control_gaussians[1,]),
                           var_val = sum(control_gaussians[2,])/ncol(control_gaussians)^2,
                           n = ncol(control_gaussians))

  n1 <- control_gaussian$n; n2 <- case_gaussian$n
  mean1 <- control_gaussian$mean_val; mean2 <- case_gaussian$mean_val
  cov1 <- control_gaussian$var_val; cov2 <- control_gaussian$var_val
  combined_cov <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
  test_stat <- (n1*n2)*(mean1 - mean2)^2/(combined_cov * (n1+n2))

  # using https://github.com/cran/Hotelling/blob/master/R/hotelling.test.r
  # and https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
  p <- 1; m <- n1+n2-2
  p_val <- 1-stats::pf((m-p+1)/(p*m)*test_stat, df1 = p, df2 = m-p+1)

  log_diff <- log(mean(mat[case_idx,j])) - log(mean(mat[control_idx,j]))

  c(p_val = p_val, log_diff = log_diff)
})

hk_genes <- read.csv("../../data/housekeeping.txt", header = F)[,1]
de_genes <- c("TTF2", "MX2", "ASCC1", "GLRA3", "CIRBP",
              "SAT2", "QTRT1", "CDH2", "LUC7L", "TCF25",
              "SSBP2", "WDR60", "CABP1", "FBLN7", "CDC14B",
              "GPM6A", "IGFBP5", "FAM153B", "GUCY1A2", "RAB3C",
              "SSX2IP", "HS6ST3", "TENM3", "DACH1", "PLA2G4C",
              "TOX3", "SPAG16", "FAM171B", "GALNTL6", "NUMB",
              "CAPZB", "DDRGK1", "RMST", "SUGP2", "FAM49A",
              "KCNH7", "BRINP3", "GABRB1", "GOLGA8B", "OR2L13",
              "IMMP2L", "ARPP19", "VWA8", "RPS15", "DPYSL2",
              "RFX3", "RSRP1", "NFIA", "SNRNP70", "SYN2",
              "SPIN1", "PLPPR4", "SYNPR", "SLC22A10", "LINC01378",
              "RP11-577H5.5", "GABRG2", "MIR99AHG", "PPP3CA",
              "MIR137HG", "TBRG1", "GGT7", "NLGN1",
              "GNG7", "FZD3", "LRRTM3", "CPE", "KCNJ3",
              "AQP4-AS1", "TRAF3", "PKIA", "MGAT4C", "HNRNPDL",
              "SLITRK4", "BMPR1B", "AHI1", "CDH9", "RAPGEFL1",
              "RPL34P18", "LINC00657", "COL26A1", "CNTN3", "FRMD6",
              "RP11-30J20.1", "SLITRK5", "SLC39A10", "STX1A", "RPLP2",
              "MAP2", "CES4A", "NEGR1", "SORBS1", "COL24A1",
              "VSTM2L", "ERBB4", "STARD4-AS1", "MAPK1", "HSP90AA1",
              "RPL34", "CNTNAP2", "EIF1", "OLFM3", "GRID2",
              "CHL1", "RAP1GAP", "CAMK2N1", "SERINC1", "RGS12",
              "ATP1B1")
de_genes <- de_genes[de_genes %in% colnames(mat)]
hk_idx <- which(colnames(mat) %in% hk_genes)
de_idx <- which(colnames(mat) %in% de_genes)
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), ncol(summary_stats))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
tol <- 1e-12
plot(NA, xlim = range(summary_stats[2,]), ylim = range(-log10(summary_stats[1,]+tol)))
points(x = summary_stats[2, -unique(c(hk_idx,de_idx))],
       y = -log10(summary_stats[1,-unique(c(hk_idx,de_idx))]+tol),
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = summary_stats[2,hk_idx],
       y = -log10(summary_stats[1,hk_idx]+tol),
       pch = 16, col = col_vec[hk_idx])
points(x = summary_stats[2,de_idx],
       y = -log10(summary_stats[1,de_idx]+tol),
       pch = 16, col = col_vec[de_idx])

###########################

hist(esvd_res$b_mat[,"diagnosis_ASD"])
rug(esvd_res$b_mat[which(colnames(mat) %in% hk_genes),"diagnosis_ASD"], col = 3)
rug(esvd_res$b_mat[which(colnames(mat) %in% de_genes),"diagnosis_ASD"], col = 2)

#########################

case_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 1)
control_idx <- which(esvd_res$covariates[,"diagnosis_ASD"] == 0)
p_val_vec <- sapply(1:ncol(posterior_mean_mat2), function(j){
  res <- stats::wilcox.test(x = posterior_mean_mat2[case_idx,j],
                            y = posterior_mean_mat2[control_idx,j])
  res$p.value
})
quantile(p_val_vec)
col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), ncol(summary_stats))
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
tol <- 1e-12
plot(NA, xlim = range(summary_stats[2,]), ylim = range(-log10(p_val_vec + tol)))
points(x = summary_stats[2, -unique(c(hk_idx,de_idx))],
       y = -log10(p_val_vec[-unique(c(hk_idx,de_idx))] + tol),
       pch = 16, col = col_vec[-unique(c(hk_idx,de_idx))])
points(x = summary_stats[2,hk_idx],
       y = -log10(p_val_vec[hk_idx] + tol),
       pch = 16, col = col_vec[hk_idx])
points(x = summary_stats[2,de_idx],
       y = -log10(p_val_vec[de_idx] + tol),
       pch = 16, col = col_vec[de_idx])


quantile(-log10(p_val_vec[hk_idx]+tol))
quantile(-log10(p_val_vec[de_idx]+tol))

png("../../out/fig/writeup8e/tmp.png", height = 1500,
    width = 1500, res = 300, units = "px")
plot(x = posterior_mean_mat2[,de_idx],
     y = posterior_var_mat[,de_idx],
     asp = T, pch = 16, col = rgb(0.5, 0.5, 0.5, 0.5))
graphics.off()


