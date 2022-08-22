rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

##################

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
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

####################

plotting_func <- function(eSVD_obj,
                          cycling_genes,
                          de_genes,
                          de_gene_specific,
                          hk_genes,
                          sfari_genes){
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
}

find_de_genes <- function(eSVD_obj,
                          seurat_obj = sns,
                          covariate_individual = "individual"){
  set.seed(10)
  df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj,
                               metadata = seurat_obj@meta.data,
                               covariate_individual = covariate_individual)
  teststat_vec <- eSVD_obj$teststat_vec
  p <- length(teststat_vec)
  gaussian_teststat <- sapply(1:p, function(j){
    qnorm(pt(teststat_vec[j], df = df_vec[j]))
  })

  locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
  fdr_vec <- locfdr_res$fdr
  names(fdr_vec) <- names(gaussian_teststat)
  null_mean <- locfdr_res$fp0["mlest", "delta"]
  null_sd <- locfdr_res$fp0["mlest", "sigma"]
  print(paste0("Null mean: ", round(null_mean,2), ", null sd: ", round(null_sd, 2)))

  logpvalue_vec <- sapply(gaussian_teststat, function(x){
    if(x < null_mean) {
      Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
    } else {
      Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
    }
  })
  logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))

  # selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
  selected_genes <- names(fdr_vec)[order(fdr_vec, decreasing = F)[1:100]]
  selected_genes
}

original_selected_genes <- find_de_genes(eSVD_obj)
length(original_selected_genes)

downsample_values <- c(0.9, 0.8, 0.7, 0.6, 0.5)
downsampled_selected_genes <- lapply(downsample_values, function(downsample_value){
  print(downsample_value)
  load(paste0("../../../out/main/sns_layer23_esvd_downsampled-", downsample_value, ".RData"))
  print("Loaded")

  # png(paste0("../../../out/fig/main/sns_layer23_esvd_downsampled-", downsample_value, "_teststat_histogram.png"),
  #     height = 1200, width = 1200,
  #     units = "px", res = 300)
  # plotting_func(eSVD_obj,
  #               cycling_genes,
  #               de_genes,
  #               de_gene_specific,
  #               hk_genes,
  #               sfari_genes)
  # graphics.off()
  find_de_genes(eSVD_obj)
})


save(sns, original_selected_genes, downsampled_selected_genes,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_eSVD_downsampled_genes.RData")


######################

sapply(downsampled_selected_genes, length)
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, original_selected_genes))
})

length(intersect(original_selected_genes, sfari_genes))
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, sfari_genes))
})

length(intersect(original_selected_genes, bulk_de_genes))
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, bulk_de_genes))
})

length(intersect(original_selected_genes, hk_genes))
sapply(downsampled_selected_genes, function(x){
  length(intersect(x, hk_genes))
})

# ################################
#
# load("../../../out/main/sns_layer23_esvd.RData")
#
# find_de_genes_semisupervised <- function(eSVD_obj,
#                                          seurat_obj = sns,
#                                          covariate_individual = "individual"){
#   df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj,
#                                metadata = seurat_obj@meta.data,
#                                covariate_individual = covariate_individual)
#   teststat_vec <- eSVD_obj$teststat_vec
#   p <- length(teststat_vec)
#   gaussian_teststat <- sapply(1:p, function(j){
#     qnorm(pt(teststat_vec[j], df = df_vec[j]))
#   })
#
#   locfdr_res <- locfdr::locfdr(gaussian_teststat, plot = 0)
#   fdr_vec <- locfdr_res$fdr
#   names(fdr_vec) <- names(gaussian_teststat)
#   null_mean <- locfdr_res$fp0["mlest", "delta"]
#   null_sd <- locfdr_res$fp0["mlest", "sigma"]
#   gaussian_teststat <- (gaussian_teststat - null_mean)/null_sd
#
#   semisupervised_multtest(teststat_vec = gaussian_teststat,
#                           null_idx = which(names(teststat_vec) %in% hk_genes),
#                           alpha = 0.05)
# }
#
# original_selected_genes2 <- find_de_genes_semisupervised(eSVD_obj)
# length(original_selected_genes2)
#
# downsample_values <- c(0.9, 0.8, 0.7, 0.6, 0.5)
# downsampled_selected_genes2 <- lapply(downsample_values, function(downsample_value){
#   print(downsample_value)
#   load(paste0("../../../out/main/sns_layer23_esvd_downsampled-", downsample_value, ".RData"))
#   print("Loaded")
#
#   find_de_genes_semisupervised(eSVD_obj)
# })
#
# sapply(downsampled_selected_genes2, length)
# sapply(downsampled_selected_genes2, function(x){
#   length(intersect(x, original_selected_genes2))
# })
