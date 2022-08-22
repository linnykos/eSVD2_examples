rm(list=ls())
library(Seurat)
library(eSVD2)

pie_custom <- function(x, offset = c(0,0), edges = 200, radius = 0.8,
                       clockwise = T,
                       init.angle = 90,
                       col = 1:length(x),
                       border = rep(1, length(x)), lwd = 1){
  if (!is.numeric(x) || any(is.na(x) | x < 0))
    stop("'x' values must be positive.")
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  if (length(border) == 1) border <- rep_len(border, nx)
  if (!is.null(lwd)) lwd <- rep_len(lwd, nx)
  twopi <- ifelse(clockwise, -2 * pi, 2 * pi)
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }

  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    graphics::polygon(c(P$x, 0) + offset[1], c(P$y, 0) + offset[2],
                      border = border[i], col = col[i], lwd = lwd[i])
    P <- t2xy(mean(x[i + 0:1]))
  }

  invisible()
}

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

de_gene_specific <- intersect(de_gene_specific, colnames(eSVD_obj$dat))
sfari_genes <- intersect(sfari_genes, colnames(eSVD_obj$dat))
bulk_de_genes <- intersect(bulk_de_genes, colnames(eSVD_obj$dat))
hk_genes <- intersect(hk_genes, colnames(eSVD_obj$dat))
lis <- list(de_gene_specific, sfari_genes, bulk_de_genes, hk_genes)
mat <- matrix(NA, 4, 4)
for(i in 1:4){
  for(j in 1:4){
    mat[i,j] <- length(intersect(lis[[i]], lis[[j]]))
  }
}
colnames(mat) <- c("Velmeshev", "SFARI", "Bulk", "HK")
rownames(mat) <- colnames(mat)
mat

sfari_genes <- setdiff(sfari_genes, de_gene_specific)
bulk_de_genes <- setdiff(bulk_de_genes, c(sfari_genes, de_gene_specific))
hk_genes <- setdiff(hk_genes, c(sfari_genes, de_gene_specific, bulk_de_genes))

#############

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
brown_col <- rgb(144, 100, 43,  maxColorValue = 255)

load("../../../out/main/sns_layer23_eSVD_downsampled_genes.RData")

downsampled_selected_genes <- c(list(original_selected_genes), downsampled_selected_genes)
downsample_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

for(i in 1:length(downsample_values)){
  downsample_value <- downsample_values[i]
  print(downsample_value)

  x_vec <- c(length(intersect(downsampled_selected_genes[[i]], de_gene_specific)),
             length(intersect(downsampled_selected_genes[[i]], sfari_genes)),
             length(intersect(downsampled_selected_genes[[i]], bulk_de_genes)),
             length(setdiff(downsampled_selected_genes[[i]], c(de_gene_specific, sfari_genes, bulk_de_genes, hk_genes))),
             length(intersect(downsampled_selected_genes[[i]], hk_genes))
  )

  png(paste0("../../../out/fig/main/sns_layer23_esvd_pieplot_downsample-", downsample_value, ".png"),
      height = 1500, width = 1500,
      units = "px", res = 500)
  par(mar = c(0,0,0,0))
  plot(NA, xlim = c(-5,5), ylim = c(-5,5),
       xaxt = "n", yaxt = "n", bty = "n", asp = T,
       xlab = "", ylab = "")
  pie_custom(x = x_vec, radius = 4.5,
             col = c(brown_col, purple_col, blue_col, rgb(0.6, 0.6, 0.6), green_col))
  graphics.off()
}

##################################

library(SummarizedExperiment)
library(DESeq2)
load("../../../out/main/sns_layer23_deseq2.RData")
load("../../../out/main/sns_layer23_deseq2_downsampled.RData")

original_selected_genes <- rownames(deseq2_res)[order(deseq2_res$padj, decreasing = F)[1:100]]
downsampled_selected_genes <- lapply(de_result_downsampled, function(x){
  rownames(x)[order(x$padj, decreasing = F)[1:100]]
})
downsampled_selected_genes <- c(list(original_selected_genes), downsampled_selected_genes)
downsample_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

for(i in 1:length(downsample_values)){
  downsample_value <- downsample_values[i]
  print(downsample_value)

  x_vec <- c(length(intersect(downsampled_selected_genes[[i]], de_gene_specific)),
             length(intersect(downsampled_selected_genes[[i]], sfari_genes)),
             length(intersect(downsampled_selected_genes[[i]], bulk_de_genes)),
             length(setdiff(downsampled_selected_genes[[i]], c(de_gene_specific, sfari_genes, bulk_de_genes, hk_genes))),
             length(intersect(downsampled_selected_genes[[i]], hk_genes))
  )

  png(paste0("../../../out/fig/main/sns_layer23_deseq2_pieplot_downsample-", downsample_value, ".png"),
      height = 1500, width = 1500,
      units = "px", res = 500)
  par(mar = c(0,0,0,0))
  plot(NA, xlim = c(-5,5), ylim = c(-5,5),
       xaxt = "n", yaxt = "n", bty = "n", asp = T,
       xlab = "", ylab = "")
  pie_custom(x = x_vec, radius = 4.5,
             col = c(brown_col, purple_col, blue_col, rgb(0.6, 0.6, 0.6), green_col))
  graphics.off()
}

##############################

load("../../../out/main/sns_layer23_sctransform_downsampled.RData")

original_selected_genes <- rownames(de_result)[order(de_result$p_val_adj, decreasing = F)[1:100]]
downsampled_selected_genes <- lapply(de_result_downsampled, function(x){
  rownames(x)[order(x$p_val_adj, decreasing = F)[1:100]]
})

downsampled_selected_genes <- c(list(original_selected_genes), downsampled_selected_genes)
downsample_values <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5)

for(i in 1:length(downsample_values)){
  downsample_value <- downsample_values[i]
  print(downsample_value)

  x_vec <- c(length(intersect(downsampled_selected_genes[[i]], de_gene_specific)),
             length(intersect(downsampled_selected_genes[[i]], sfari_genes)),
             length(intersect(downsampled_selected_genes[[i]], bulk_de_genes)),
             length(setdiff(downsampled_selected_genes[[i]], c(de_gene_specific, sfari_genes, bulk_de_genes, hk_genes))),
             length(intersect(downsampled_selected_genes[[i]], hk_genes))
  )

  png(paste0("../../../out/fig/main/sns_layer23_sctransform_pieplot_downsample-", downsample_value, ".png"),
      height = 1500, width = 1500,
      units = "px", res = 500)
  par(mar = c(0,0,0,0))
  plot(NA, xlim = c(-5,5), ylim = c(-5,5),
       xaxt = "n", yaxt = "n", bty = "n", asp = T,
       xlab = "", ylab = "")
  pie_custom(x = x_vec, radius = 4.5,
             col = c(brown_col, purple_col, blue_col, rgb(0.6, 0.6, 0.6), green_col))
  graphics.off()
}
