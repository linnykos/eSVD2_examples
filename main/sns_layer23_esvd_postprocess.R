rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

eSVD_obj$fit_Second$posterior_mean_mat <- NULL
eSVD_obj$fit_Second$posterior_var_mat <- NULL
eSVD_obj$teststat_vec <- NULL
eSVD_obj <- eSVD2:::compute_posterior(input_obj = eSVD_obj,
                                      bool_adjust_covariates = F,
                                      alpha_max = NULL,
                                      bool_covariates_as_library = T,
                                      library_min = 1e-1)
metadata <- sns@meta.data
metadata[,"individual"] <- as.factor(metadata[,"individual"])
time_start5 <- Sys.time()
eSVD_obj <- eSVD2:::compute_test_statistic(input_obj = eSVD_obj,
                                           covariate_individual = "individual",
                                           metadata = metadata,
                                           verbose = 1)

#################

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
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

hk_idx <- which(names(eSVD_obj$teststat_vec) %in% c(hk_genes, cycling_genes))
de_idx <- which(names(eSVD_obj$teststat_vec) %in% de_gene_specific)
other_idx <- which(names(eSVD_obj$teststat_vec) %in% c(sfari_genes, de_genes))

col_vec <- rep(rgb(0.5, 0.5, 0.5, 0.5), length(eSVD_obj$teststat_vec))
col_vec[other_idx] <- 4
col_vec[hk_idx] <- 3
col_vec[de_idx] <- 2
shuf_idx <- c(hk_idx, de_idx, other_idx)
shuf_idx <- shuf_idx[sample(length(shuf_idx))]

teststat_vec <- pmax(pmin(eSVD_obj$teststat_vec, 30), -30)
max_val <- max(abs(eSVD_obj$teststat_vec))
png("../../../out/fig/main/sns_layer23_esvd_teststat_histogram.png",
    height = 1200, width = 1200,
    units = "px", res = 300)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
hist(teststat_vec, breaks = break_vec,
     xlim = c(-max_val, max_val),
     main = "Histogram of test statistic",
     xlab = "Z-score", ylab = "Frequency", freq = T)
lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
for(i in shuf_idx){
  rug(teststat_vec[i], col = col_vec[i], lwd = 2)
}
legend("topright", c("Published DE gene", "Other interest gene", "Housekeeping gene"),
       fill = c(2,4,3), cex = 0.6)
graphics.off()

png("../../../out/fig/main/sns_layer23_esvd_teststat_histogram_separate.png",
    height = 1000, width = 3000,
    units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
uniq_col_vec <- c(2,4,3)
break_vec <- seq(-max_val-0.15, max_val+0.15, by = 0.1)
main_vec <- c("Published DE", "Other interest\n(SFARI, DE other region)", "Housekeeping+Cell-cycle")
for(kk in 1:length(uniq_col_vec)){
  idx <- which(col_vec == uniq_col_vec[kk])
  hist(teststat_vec[idx], breaks = break_vec,
       xlim = c(-max_val, max_val),
       main = paste0("Z-scores: ", main_vec[kk]),
       xlab = "Z-score", ylab = "Frequency", freq = T)
  lines(rep(0,2), c(0, 1e5), lwd = 1, lty = 3)
  rug(teststat_vec[idx], col = col_vec[idx], lwd = 2)
}
graphics.off()

########################################

gene_names <- names(eSVD_obj$teststat_vec)
hk_idx <- setdiff(hk_idx, de_idx)
other_idx <- setdiff(other_idx, c(hk_idx, de_idx))

png(paste0("../../../out/fig/main/sns_layer23_diagnostic_gene.png"),
    height = 2500, width = 2500,
    units = "px", res = 300)
par(mfrow = c(2,2), mar = c(4,4,4,0.5))
eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "teststat",
                  gene_list = list(gene_names[de_idx],
                                   gene_names[hk_idx],
                                   gene_names[other_idx]),
                  color_palette = c(2,3,4))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "nuisance",
                  what_2 = "Log_UMI",
                  gene_list = list(gene_names[de_idx],
                                   gene_names[hk_idx],
                                   gene_names[other_idx]),
                  color_palette = c(2,3,4))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "diagnosis_ASD",
                  what_2 = "teststat",
                  gene_list = list(gene_names[de_idx],
                                   gene_names[hk_idx],
                                   gene_names[other_idx]),
                  color_palette = c(2,3,4))

eSVD2:::gene_plot(eSVD_obj,
                  what_1 = "sparsity",
                  what_2 = "teststat",
                  gene_list = list(gene_names[de_idx],
                                   gene_names[hk_idx],
                                   gene_names[other_idx]),
                  color_palette = c(2,3,4))

graphics.off()


# select 10 genes with smallest test-statistic in absolute values
set.seed(10)
gene_idx <- order(abs(eSVD_obj$teststat_vec), decreasing = F)[1:10]

tab <- table(sns$individual, sns$diagnosis)
indiv_cases <- rownames(tab)[which(tab[,"ASD"] != 0)]
indiv_controls <- rownames(tab)[which(tab[,"Control"] != 0)]
indiv_vec <- factor(as.character(sns$individual))

for(k in 1:length(gene_idx)){
  print(k)

  idx <- gene_idx[k]
  tmp <- as.numeric(eSVD_obj$dat[,idx])

  png(paste0("../../../out/fig/main/sns_layer23_esvd_diagnostic_nearzero-",
             k, ".png"),
      height = 2500, width = 2500,
      units = "px", res = 300)
  par(mfrow = c(2,2), mar = c(4,4,4,0.5))

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "dat",
                    bool_jitter_y = T,
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0(gene_names[idx], ": Zero% = ", round(100*length(which(tmp == 0))/length(tmp)))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "adjusted_relative_expression",
                    bool_jitter_y = T,
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Nuisance = ", round(eSVD_obj$fit_Second$nuisance_vec[idx], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "fit",
                    what_2 = "posterior_mean",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("ASD coefficient = ",
                                  round(eSVD_obj$fit_Second$z_mat[idx,"diagnosis_ASD"], 2))
  )

  eSVD2:::cell_plot(eSVD_obj,
                    variable = gene_names[idx],
                    what_1 = "posterior_mean",
                    what_2 = "posterior_sd",
                    indiv_cases = indiv_cases,
                    indiv_controls = indiv_controls,
                    indiv_vec = indiv_vec,
                    main = paste0("Test statistic = ", round(eSVD_obj$teststat_vec[idx], 2))
  )

  graphics.off()
}

