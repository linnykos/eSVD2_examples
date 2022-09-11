rm(list=ls())
library(Seurat)
library(eSVD2)

# load("../../../out/main/sns_layer23_esvd.RData")
load("../../../out/Writeup12/Writeup12_sns_layer23_esvd3.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

df_vec <- eSVD2:::compute_df(input_obj = eSVD_obj)
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
logpvalue_vec <- sapply(gaussian_teststat, function(x){
  if(x < null_mean) {
    Rmpfr::pnorm(x, mean = null_mean, sd = null_sd, log.p = T)
  } else {
    Rmpfr::pnorm(null_mean - (x-null_mean), mean = null_mean, sd = null_sd, log.p = T)
  }
})
logpvalue_vec <- -(logpvalue_vec/log(10) + log10(2))
logpvalue_vec <- pmin(logpvalue_vec, 15)

selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.0001)]

############

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
length(bulk_de_genes)
length(intersect(deg_df[,"external_gene_name"], colnames(eSVD_obj$dat)))
length(intersect(bulk_de_genes, sfari_genes))
length(intersect(bulk_de_genes, de_gene_specific))
length(intersect(bulk_de_genes, de_genes))
length(intersect(bulk_de_genes, hk_genes))

# ## see https://www.simplypsychology.org/brodmann-areas.html for brodmann areas
# col_idx <- grep("ASD_BA.*_FDR", colnames(deg_df))
# for(j in col_idx){
#   print(colnames(deg_df)[j])
#   tmp <- deg_df[which(deg_df[,j]<=0.05),"external_gene_name"]
#   print(length(tmp))
#   print(length(intersect(tmp, sfari_genes)))
# }
#
# bulk_de_genes2 <- deg_df[which(deg_df[,"ASD_BA41_42_22_FDR"]<=0.01),"external_gene_name"]
# length(bulk_de_genes2)

length(selected_genes)
length(intersect(selected_genes, de_gene_specific))
length(intersect(selected_genes, de_genes))
length(intersect(selected_genes, hk_genes))
length(intersect(selected_genes, sfari_genes))
length(intersect(selected_genes, cycling_genes))
length(intersect(selected_genes, bulk_de_genes))
length(intersect(selected_genes, bulk_de_genes2))


## https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
m <- length(sfari_genes)
n <- length(gaussian_teststat) - m
k <- length(selected_genes)
x <- length(intersect(selected_genes, c(sfari_genes)))
fisher <- stats::dhyper(x = x, m = m, n = n, k = k, log = F)
fisher

m <- length(bulk_de_genes)
n <- length(gaussian_teststat) - m
k <- length(selected_genes)
x <- length(intersect(selected_genes, c(bulk_de_genes)))
fisher <- stats::dhyper(x = x, m = m, n = n, k = k, log = F)
fisher

quantile(logpvalue_vec[intersect(sfari_genes, colnames(eSVD_obj$dat))])
quantile(logpvalue_vec[intersect(deg_df[,"external_gene_name"], colnames(eSVD_obj$dat))])

#######

idx <- which(names(eSVD_obj$teststat_vec) %in% selected_genes)
sfari_idx <- which(names(eSVD_obj$teststat_vec) %in% sfari_genes)
hk_idx <- which(names(eSVD_obj$teststat_vec) %in% hk_genes)

x_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)
xlim <- quantile(x_vec, probs = c(0.01, 1))
xlim <- c(-1,1)*max(abs(xlim))
y_vec <- logpvalue_vec
ylim <- range(y_vec)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
purple_col <- rgb(122, 49, 126, maxColorValue = 255)
blue_col <- rgb(129, 139, 191, maxColorValue = 255)
green_col <- rgb(70, 177, 70, maxColorValue = 255)
green_col_trans <- rgb(70, 177, 70, 255*.35, maxColorValue = 255)

# adjustment of selected genes for visual clarity
min_pthres <- min(y_vec[selected_genes])
selected_genes2 <- names(y_vec)[y_vec >= min_pthres]
idx <- which(names(y_vec) %in% selected_genes2)
sfari_idx <- which(names(y_vec) %in% sfari_genes)
bulk_idx <- which(names(y_vec) %in% bulk_de_genes)
hk_idx <- which(names(y_vec) %in% hk_genes)

png("../../../out/fig/main/sns_layer23_volcano.png",
    height = 3500, width = 2500,
    units = "px", res = 500)
par(mar = c(3,3,0.4,0.1))
plot(x = x_vec, y = y_vec,
     xaxt = "n", yaxt = "n", bty = "n",
     ylim = ylim, xlim = xlim,
     cex.lab = 1.25, type = "n")
for(j in seq(0,15,by = .5)){
  lines(x = c(-1e4,1e4), y = rep(j, 2), col = "gray", lty = 2, lwd = 1)
}
lines(x = c(-1e4,1e4), y = rep(min_pthres, 2), col = orange_col, lty = 2, lwd = 2)

points(x = x_vec, y = y_vec,
       col = rgb(0.6, 0.6, 0.6, 0.3), pch = 16)
points(x = x_vec[idx], y = y_vec[idx],
       col = orange_col, pch = 16, cex = 1.5)

# plot non-overlapping genes
points(x = x_vec[setdiff(bulk_idx, idx)], y = y_vec[setdiff(bulk_idx, idx)],
       col = blue_col, pch = 16, cex = 1, lwd = 2)
points(x = x_vec[setdiff(sfari_idx, idx)], y = y_vec[setdiff(sfari_idx, idx)],
       col = purple_col, pch = 16, cex = 1, lwd = 2)

# plot housekeeping
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = "white", pch = 16, cex = 1)
points(x = x_vec[hk_idx], y = y_vec[hk_idx],
       col = green_col_trans, pch = 16, cex = 1)

# circle overlapping gens
points(x = x_vec[intersect(idx, sfari_idx)], y = y_vec[intersect(idx, sfari_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx, sfari_idx)], y = y_vec[intersect(idx, sfari_idx)],
       col = purple_col, pch = 1, cex = 2, lwd = 2)
points(x = x_vec[intersect(idx, bulk_idx)], y = y_vec[intersect(idx, bulk_idx)],
       col = "white", pch = 1, cex = 2, lwd = 3)
points(x = x_vec[intersect(idx, bulk_idx)], y = y_vec[intersect(idx, bulk_idx)],
       col = blue_col, pch = 1, cex = 2, lwd = 2)

axis(1, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
lines(x = rep(0, 2), y = c(-10,100), lwd = 1.5, lty = 3, col = 1)
graphics.off()

###############

write_genes <- function(vec,
                        file){
  fileConn <- file(file)
  writeLines(vec, fileConn)
  close(fileConn)
}

write_genes(selected_genes, file = "../../../out/main/sns_layer23_esvd_de-genes.csv")
write_genes(de_gene_specific, file = "../../../out/main/sns_layer23_velmeshev_de-genes.csv")

#################3

library(org.Hs.eg.db)
library(clusterProfiler)

## see https://github.com/YuLab-SMU/clusterProfiler/issues/222
universe_vec <-  names(eSVD_obj$teststat_vec)
esvd_res <- clusterProfiler::enrichGO(
  gene = selected_genes,
  universe = universe_vec,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
velmeshev_res <- clusterProfiler::enrichGO(
  gene = de_gene_specific,
  universe = universe_vec,
  keyType = "SYMBOL",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

colnames(esvd_res@result)
head(esvd_res@result[,c("ID", "Description", "pvalue", "GeneRatio")], 50)
generatio <- sapply(esvd_res@result[,"GeneRatio"], function(x){
  zz <- strsplit(x, split = "/")[[1]]
  as.numeric(zz[1])/as.numeric(zz[2])
})
esvd_res@result[order(generatio, decreasing = T)[1:50],c("ID", "Description", "pvalue", "GeneRatio")]
esvd_res@result[order(esvd_res@result[,"pvalue"], decreasing = F)[1:50],c("ID", "Description", "pvalue", "GeneRatio")]

esvd_res@result["GO:0099536",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")] # synaptic signaling
esvd_res@result["GO:0061564",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")] # axon development
esvd_res@result["GO:0031175",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")] # neuron projection development
esvd_res@result["GO:0050877",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")] # nervous system process

esvd_res@result["GO:0006812",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:0060341",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:1901698",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:0034220",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:0051668",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]

# esvd_res@result["GO:0048812",] # neuron projection morphogenesis
# esvd_res@result["GO:0048667",] # cell morphogenesis involved in neuron differentiation

colnames(velmeshev_res@result)
head(velmeshev_res@result[,c("ID", "Description", "pvalue", "GeneRatio")])
velmeshev_res@result["GO:0099536",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
velmeshev_res@result["GO:0061564",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
velmeshev_res@result["GO:0031175",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
velmeshev_res@result["GO:0050877",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]

go_vec <- c("GO:0099536", "GO:0061564", "GO:0031175", "GO:0050877")
esvd_pvalue <- sapply(go_vec, function(x){
  esvd_res@result[x, "pvalue"]
})
velmeshev_pvalue <- sapply(go_vec, function(x){
  velmeshev_res@result[x, "pvalue"]
})

esvd_log10pvalue <- -log10(esvd_pvalue)
velmeshev_log10pvalue <- -log10(velmeshev_pvalue)

spacing <- .5
width <- 1
len <- length(esvd_pvalue)
total_len <- 2*width*len + spacing*(len-1)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
brown_col <- rgb(144, 100, 43, maxColorValue = 255)

png("../../../out/fig/main/sns_layer23_gsea.png",
    height = 880, width = 1500,
    units = "px", res = 500)
par(mar = c(0.1,2,0.1,0.1))
plot(NA, xlim = c(0, total_len),
     ylim = c(0, max(c(esvd_log10pvalue, velmeshev_log10pvalue))),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:len){
  x_left <- (2*width+spacing)*(i-1)
  polygon(x = c(x_left, x_left+width, x_left+width, x_left),
          y = c(0, 0, esvd_log10pvalue[i], esvd_log10pvalue[i]),
          col = orange_col, border = "black")
  polygon(x = c(x_left, x_left+width, x_left+width, x_left)+width,
          y = c(0, 0, velmeshev_log10pvalue[i], velmeshev_log10pvalue[i]),
          col = brown_col, border = "black")
}
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()

####################

go_vec <- c("GO:0099536", "GO:0061564", "GO:0031175", "GO:0050877")
gene_list <- sapply(go_vec, function(x){
  vec <- esvd_res@result[x,"geneID"]
  vec2 <- strsplit(vec, split = "/")[[1]]
  vec2
})
table(unlist(gene_list))
for(vec in gene_list){
  print(intersect(vec, sfari_genes))
  print(intersect(vec, bulk_de_genes))
  print("===")
}

x_vec <- log2(eSVD_obj$case_mean) - log2(eSVD_obj$control)

for(vec in gene_list){
  zz <- x_vec[vec]
  print(paste0("Neg: ", length(which(zz < 0)), ", Pos: ", length(which(zz > 0))))
}
