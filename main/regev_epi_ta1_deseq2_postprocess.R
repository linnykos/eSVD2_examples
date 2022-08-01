rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/regevEpi_ta1-noninflamed_deseq2.RData")
noninflamed_deseq2 <- deseq2_res

load("../../../out/main/regevEpi_ta1-inflamed_deseq2.RData")
inflamed_deseq2 <- deseq2_res

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Non-Inflamed vs. He"))
sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Health"))
sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Non-In"))
noninf_de_genes <- sheet1[sheet1$ident == "TA 1","gene"]
inf_de_genes <- sheet2[sheet2$ident == "TA 1","gene"]
other_de_genes <- sheet3[sheet3$ident == "TA 1","gene"]
other_de_genes <- setdiff(other_de_genes, c(noninf_de_genes, inf_de_genes))
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)
cycling_genes <- setdiff(cycling_genes, c(other_de_genes, noninf_de_genes, inf_de_genes))
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]

###################################

inflamed_idx <- sapply(c(inf_de_genes, noninf_de_genes), function(x){
  zz <- which(rownames(inflamed_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
inflamed_de <- inflamed_deseq2[inflamed_idx,"log2FoldChange"]
noninflamed_idx <- sapply(c(inf_de_genes, noninf_de_genes), function(x){
  zz <- which(rownames(noninflamed_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
noninflamed_de <- noninflamed_deseq2[noninflamed_idx,"log2FoldChange"]
rm_idx <- unique(c(which(is.na(inflamed_de)), which(is.na(noninflamed_de))))
if(length(rm_idx) > 0){
  inflamed_de <- inflamed_de[-rm_idx]; noninflamed_de <- noninflamed_de[-rm_idx]
}
stats::cor(inflamed_de, noninflamed_de, method = "spearman")

inflamed_idx <- sapply(hk_genes, function(x){
  zz <- which(rownames(inflamed_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
inflamed_hk <- inflamed_deseq2[inflamed_idx,"log2FoldChange"]
noninflamed_idx <- sapply(hk_genes, function(x){
  zz <- which(rownames(noninflamed_deseq2)==x)
  if(length(zz) == 1) return(zz) else return(NA)
})
noninflamed_hk <- noninflamed_deseq2[noninflamed_idx,"log2FoldChange"]
rm_idx <- unique(c(which(is.na(inflamed_hk)), which(is.na(noninflamed_idx))))
if(length(rm_idx) > 0){
  inflamed_hk <- inflamed_hk[-rm_idx]; noninflamed_hk <- noninflamed_hk[-rm_idx]
}
stats::cor(inflamed_hk, noninflamed_hk, method = "spearman")

png("../../../out/fig/main/regevEpi_ta1-agreement_deseq2_de-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
xbnds <- range(c(inflamed_de, inflamed_hk))
ybnds <- range(c(noninflamed_de, noninflamed_hk))
bin <- hexbin::hexbin(inflamed_de, noninflamed_de, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(inflamed_de, noninflamed_de, method = "spearman"), 2)))
graphics.off()

png("../../../out/fig/main/regevEpi_ta1-agreement_deseq2_hk-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
xbnds <- range(c(inflamed_de, inflamed_hk))
ybnds <- range(c(noninflamed_de, noninflamed_hk))
bin <- hexbin::hexbin(inflamed_hk, noninflamed_hk, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(inflamed_hk, noninflamed_hk, method = "spearman"), 2)))
graphics.off()


