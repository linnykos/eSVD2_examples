rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

##############

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

candidate_genes <- sort(unique(c(sfari_genes, bulk_de_genes)))

###############

load("../../../out/main/sns_layer23_deseq2.RData")
load("../../../out/main/sns_layer23_deseq2_downsampled.RData")

num_genes <- 50
gene_vec <- intersect(candidate_genes, rownames(deseq2_res))
selected_deseq_genes <- gene_vec[order(deseq2_res[gene_vec,"pvalue"], decreasing = F)[1:num_genes]]
deseq2_res[selected_deseq_genes,]
hk_genes2 <- intersect(hk_genes, rownames(deseq2_res))
base_mean <- mean(abs(deseq2_res[hk_genes2, "stat"]))
base_sd <- stats::sd(abs(deseq2_res[hk_genes2, "stat"]))
deseq_zmat <- sapply(selected_deseq_genes, function(x){
  (abs(deseq2_res[x, "stat"])-base_mean)/base_sd
})

for(i in 1:length(deseq_result_downsampled)){
  hk_genes2 <- intersect(hk_genes, rownames(deseq_result_downsampled[[i]]))

  base_mean <- mean(abs(deseq_result_downsampled[[i]][hk_genes2, "stat"]))
  base_sd <- stats::sd(abs(deseq_result_downsampled[[i]][hk_genes2, "stat"]))

  tmp <- sapply(selected_deseq_genes, function(x){
    (abs(deseq_result_downsampled[[i]][x,"stat"])-base_mean)/base_sd
  })

  deseq_zmat <- cbind(deseq_zmat, tmp)
}
colnames(deseq_zmat)[-1] <- names(deseq_result_downsampled)
round(deseq_zmat,2)
colMeans(deseq_zmat)
apply(deseq_zmat, 2, median)


###############################

load("../../../out/main/sns_layer23_sctransform.RData")
load("../../../out/main/sns_layer23_sctransform_downsampled.RData")

gene_vec <- intersect(candidate_genes, rownames(de_result))
selected_sct_genes <- gene_vec[order(abs(de_result[gene_vec,"z_score"]), decreasing = T)[1:num_genes]]
de_result[selected_sct_genes,]
hk_genes2 <- intersect(hk_genes,rownames(de_result))
base_mean <- mean(abs(de_result[hk_genes2, "z_score"]))
base_sd <- stats::sd(abs(de_result[hk_genes2, "z_score"]))
sct_zmat <- sapply(selected_sct_genes, function(x){
  (abs(de_result[x, "z_score"])-base_mean)/base_sd
})

for(i in 1:length(sctransform_result_downsampled)){
  hk_genes2 <- intersect(hk_genes, rownames(sctransform_result_downsampled[[i]]))

  base_mean <- mean(abs(sctransform_result_downsampled[[i]][hk_genes2, "z_score"]))
  base_sd <- stats::sd(abs(sctransform_result_downsampled[[i]][hk_genes2, "z_score"]))

  tmp <- sapply(selected_deseq_genes, function(x){
    (abs(sctransform_result_downsampled[[i]][x,"z_score"])-base_mean)/base_sd
  })

  sct_zmat <- cbind(sct_zmat, tmp)
}

colnames(sct_zmat)[-1] <- names(sctransform_result_downsampled)
round(sct_zmat,2)
colMeans(sct_zmat)
apply(sct_zmat, 2, median)

######################################

esvd_compute_testvec <- function(eSVD_obj){
  set.seed(10)
  df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
  teststat_vec <- eSVD_obj$teststat_vec
  p <- length(teststat_vec)
  gaussian_teststat <- sapply(1:p, function(j){
    qnorm(pt(teststat_vec[j], df = df_vec[j]))
  })
  gaussian_teststat
}

load("../../../out/main/sns_layer23_esvd.RData")
esvd_testvec_1 <- esvd_compute_testvec(eSVD_obj)

gene_vec <- intersect(candidate_genes, names(esvd_testvec_1))
selected_esvd_genes <- gene_vec[order(abs(esvd_testvec_1[gene_vec]), decreasing = T)[1:num_genes]]
esvd_testvec_1[selected_esvd_genes]
hk_genes2 <- intersect(hk_genes, names(esvd_testvec_1))
base_mean <- mean(abs(esvd_testvec_1[hk_genes2]))
base_sd <- stats::sd(abs(esvd_testvec_1[hk_genes2]))
esvd_zmat <- sapply(selected_esvd_genes, function(x){
  (abs(esvd_testvec_1[x])-base_mean)/base_sd
})

downsample_values <- seq(0.95, 0.6, by = -0.05)
esvd_testvec_list <- lapply(downsample_values, function(downsample_value){
  print(downsample_value)
  load(paste0("../../../out/main/sns_layer23_esvd_downsampled-", downsample_value, ".RData"))
  print("Loaded")

  esvd_compute_testvec(eSVD_obj)
})

for(i in 1:length(esvd_testvec_list)){
  hk_genes2 <- intersect(hk_genes, names(esvd_testvec_list[[i]]))
  base_mean <- mean(abs(esvd_testvec_list[[i]][hk_genes2]))
  base_sd <- stats::sd(abs(esvd_testvec_list[[i]][hk_genes2]))
  tmp <- sapply(selected_esvd_genes, function(x){
    (abs(esvd_testvec_list[[i]][x])-base_mean)/base_sd
  })

  esvd_zmat <- cbind(esvd_zmat, tmp)
}
colnames(esvd_zmat)[-1] <- downsample_values
round(esvd_zmat,2)

colMeans(esvd_zmat)
apply(esvd_zmat, 2, median)

########################

esvd_med <- apply(esvd_zmat, 2, mean)
deseq_med <- apply(deseq_zmat, 2, mean)
sct_med <- apply(sct_zmat, 2, mean)

esvd_lower <- apply(esvd_zmat, 2, stats::quantile, probs = 0.25)
deseq_lower <- apply(deseq_zmat, 2, stats::quantile, probs = 0.25)
sct_lower <- apply(sct_zmat, 2, stats::quantile, probs = 0.25)

esvd_upper <- apply(esvd_zmat, 2, stats::quantile, probs = 0.75)
deseq_upper  <- apply(deseq_zmat, 2, stats::quantile, probs = 0.75)
sct_upper  <- apply(sct_zmat, 2, stats::quantile, probs = 0.75)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
yellow_col <- rgb(255, 205, 114, maxColorValue = 255)
blue_col <- rgb(48, 174, 255, maxColorValue = 255)

orange_col_trans <- rgb(235, 134, 47, 0.4*255, maxColorValue = 255)
yellow_col_trans <- rgb(255, 205, 114, 0.4*255, maxColorValue = 255)
blue_col_trans <- rgb(48, 174, 255, 0.4*255, maxColorValue = 255)

n <- ncol(esvd_zmat)
ylim <- range(c(esvd_med, deseq_med, sct_med,
                esvd_lower, deseq_lower, sct_lower,
                esvd_upper, deseq_upper, sct_upper))
png(paste0("../../../out/fig/main/sns_layer23_downsample-signaldropff.png"),
    height = 1650, width = 3000,
    units = "px", res = 500)
par(mar = c(4,6,0,1))
plot(NA, xlim = c(1,n), ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n", xlab = "Downsample percentage",
     ylab = "Magnitude of DE genes\n(relative to HK genes)")

for(j in seq(0,10,by = .5)){
  lines(x = rep(j, 2), y = c(-1e4,1e4), col = "gray", lty = 2, lwd = 1)
}
for(j in seq(-10,10,by = .5)){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
lines(x = c(-1e4,1e4), y = rep(0, 2),  lwd = 2)


polygon(x = c(1:n, n:1),
        y = c(sct_upper, rev(sct_lower)),
        col = blue_col_trans,
        border = NA)
polygon(x = c(1:n, n:1),
        y = c(deseq_upper, rev(deseq_lower)),
        col = yellow_col_trans,
        border = NA)
polygon(x = c(1:n, n:1),
        y = c(esvd_upper, rev(esvd_lower)),
        col = orange_col_trans,
        border = NA)

lines(x = 1:n, y = sct_med, col = 1, lwd = 6)
lines(x = 1:n, y = sct_med, col = blue_col, lwd = 4, lty = 2)

lines(x = 1:n, y = deseq_med, col = 1, lwd = 6)
lines(x = 1:n, y = deseq_med, col = yellow_col, lwd = 4, lty = 2)

lines(x = 1:n, y = esvd_med, col = orange_col, lwd = 6)

axis(1, at = 1:n, labels = paste0(seq(0, 40, by = 5), "%"),
     cex.axis = 1, cex.lab = 1, lwd = 2)
axis(2, cex.axis = 1, cex.lab = 1, lwd = 2)
graphics.off()

