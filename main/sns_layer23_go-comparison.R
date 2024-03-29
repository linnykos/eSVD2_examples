rm(list=ls())
library(Seurat)
library(eSVD2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(DESeq2)


load("../../../out/main/sns_layer23_esvd.RData")
load("../../../out/main/sns_layer23_deseq2.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
esvd_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
gene_names <- names(fdr_vec)

#################

deseq_fdr_val <- stats::p.adjust(deseq2_res$pvalue, method = "BH")
names(deseq_fdr_val) <- rownames(deseq2_res)
deseq_genes <- names(deseq_fdr_val)[which(deseq_fdr_val <= 0.05)]

################

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
velmeshev_genes <- tmp[,"Gene name"]
velmeshev_genes <- velmeshev_genes[velmeshev_genes %in% gene_names]

deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% gene_names]

################

m <- length(bulk_de_genes)
n <- length(gene_names) - m
k <- length(esvd_genes)
x <- length(intersect(esvd_genes, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(gene_names)),1) )
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher

m <- length(bulk_de_genes)
n <- length(gene_names) - m
k <- length(deseq_genes)
x <- length(intersect(deseq_genes, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(gene_names)),1) )
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher

m <- length(bulk_de_genes)
n <- length(gene_names) - m
k <- length(velmeshev_genes)
x <- length(intersect(velmeshev_genes, bulk_de_genes))
paste0("#Author: ", m, ", #Bg: ", n, ", #Select: ", k, ", #Intersect: ", x,
       ", #Expected: ", round(m*(k/length(gene_names)),1) )
fisher <- sum(sapply(x:k, function(i){
  stats::dhyper(x = i, m = m, n = n, k = k, log = F)
}))
fisher

######################

esvd_ego <- clusterProfiler::enrichGO(gene          = esvd_genes,
                                      universe      = gene_names,
                                      OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                      keyType       = "SYMBOL",
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff  = 0.05,
                                      qvalueCutoff  = 0.05,
                                      readable      = TRUE)
zz <- esvd_ego@result
esvd_terms <- zz[which(zz$qvalue <= 0.05),"ID"]
# zz[which(zz$ID %in% esvd_terms),"Description"]

deseq_ego <- clusterProfiler::enrichGO(gene          = deseq_genes,
                                       universe      = gene_names,
                                       OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                       keyType       = "SYMBOL",
                                       ont           = "BP",
                                       pAdjustMethod = "BH",
                                       pvalueCutoff  = 0.05,
                                       qvalueCutoff  = 0.05,
                                       readable      = TRUE)
zz3 <- deseq_ego@result
deseq_terms <- zz3[which(zz3$qvalue <= 0.05),"ID"]

velmeshev_ego <- clusterProfiler::enrichGO(gene          = velmeshev_genes,
                                           universe      = gene_names,
                                           OrgDb         = org.Hs.eg.db::org.Hs.eg.db,
                                           keyType       = "SYMBOL",
                                           ont           = "BP",
                                           pAdjustMethod = "BH",
                                           pvalueCutoff  = 0.05,
                                           qvalueCutoff  = 0.05,
                                           readable      = TRUE)
zz2 <- velmeshev_ego@result
velmeshev_terms <- zz2[which(zz2$qvalue <= 0.05),"ID"]

###############

all_terms <- sort(unique(c(esvd_terms, velmeshev_terms, deseq_terms)))
all_terms <- intersect(intersect(intersect(all_terms, zz$ID), zz2$ID), zz3$ID)

esvd_better <- all_terms[zz[all_terms,"qvalue"] < zz2[all_terms,"qvalue"]]
velmeshev_better <- all_terms[zz[all_terms,"qvalue"] > zz2[all_terms,"qvalue"]]
deseq_better <- all_terms[zz[all_terms,"qvalue"] > zz3[all_terms,"qvalue"]]

###############

go_vec <- esvd_terms[c(1,5,6,25)]
esvd_pvalue <- sapply(go_vec, function(x){
  -log10(zz[x, "pvalue"])
})
velmeshev_pvalue <- sapply(go_vec, function(x){
  -log10(zz2[x, "pvalue"])
})
deseq_pvalue <- sapply(go_vec, function(x){
  -log10(zz3[x, "pvalue"])
})


spacing <- .5
width <- 1
len <- length(esvd_pvalue)
total_len <- 3*width*len + spacing*(len-1)

orange_col <- rgb(235, 134, 47, maxColorValue = 255)
brown_col <- rgb(144, 100, 43, maxColorValue = 255)
yellow_col <- rgb(255, 205, 114, maxColorValue = 255)

png("../../../out/fig/main/sns_layer23_gsea.png",
    height = 880, width = 1500,
    units = "px", res = 500)
par(mar = c(0.1,2,0.1,0.1))
plot(NA, xlim = c(0, total_len),
     ylim = c(0, max(c(esvd_pvalue, velmeshev_pvalue, deseq_pvalue))),
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:len){
  x_left <- (3*width+spacing)*(i-1)
  polygon(x = c(x_left, x_left+width, x_left+width, x_left),
          y = c(0, 0, esvd_pvalue[i], esvd_pvalue[i]),
          col = orange_col, border = "black")
  polygon(x = c(x_left, x_left+width, x_left+width, x_left)+width,
          y = c(0, 0, velmeshev_pvalue[i], velmeshev_pvalue[i]),
          col = brown_col, border = "black")
  polygon(x = c(x_left, x_left+width, x_left+width, x_left)+2*width,
          y = c(0, 0, deseq_pvalue[i], deseq_pvalue[i]),
          col = yellow_col, border = "black")
}
axis(2, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
graphics.off()
