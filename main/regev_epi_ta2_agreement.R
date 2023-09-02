rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/regevEpi_ta2-inflamed_esvd.RData")
eSVD_obj_inflamed <- eSVD_obj
load("../../../out/main/regevEpi_ta2-noninflamed_esvd.RData")
eSVD_obj_noninflamed <- eSVD_obj

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Non-Inflamed vs. He"))
sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Health"))
sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                          sheet = "Epithelial (Inflamed vs. Non-In"))
noninf_de_genes <- sheet1[sheet1$ident == "TA 2","gene"]
inf_de_genes <- sheet2[sheet2$ident == "TA 2","gene"]
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]

gene_names <- names(eSVD_obj$teststat_vec)
inf_de_idx <- which(gene_names %in% inf_de_genes)
noninf_de_idx <- which(gene_names %in% noninf_de_genes)
hk_idx <- which(gene_names %in% hk_genes)
hk_idx <- setdiff(hk_idx, c(inf_de_idx, noninf_de_idx))

######################

gaussian_teststat_noninflamed <- eSVD_obj_noninflamed$teststat_vec
gaussian_teststat_inflamed <- eSVD_obj_inflamed$teststat_vec

#####################

png("../../../out/fig/main/regevEpi_ta2-agreement_de-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
y1 <- gaussian_teststat_noninflamed[unique(c(inf_de_idx, noninf_de_idx))]
y2 <- gaussian_teststat_inflamed[unique(c(inf_de_idx, noninf_de_idx))]
xbnds <- range(gaussian_teststat_noninflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
ybnds <- range(gaussian_teststat_inflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
bin <- hexbin::hexbin(y1, y2, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()
stats::cor(y1, y2, method = "spearman")

png("../../../out/fig/main/regevEpi_ta2-agreement_hk-genes.png",
    height = 1750, width = 1750,
    units = "px", res = 500)
par(mar = c(3,3,0.1,0.1))
y1 <- gaussian_teststat_noninflamed[hk_idx]
y2 <- gaussian_teststat_inflamed[hk_idx]
xbnds <- range(gaussian_teststat_noninflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
ybnds <- range(gaussian_teststat_inflamed[c(inf_de_idx, noninf_de_idx, hk_idx)])
bin <- hexbin::hexbin(y1, y2, xbins = 15, xbnds = xbnds, ybnds = ybnds)
my_colors <- colorRampPalette(viridis::viridis(11))
hexbin::plot(bin, colramp=my_colors , legend=F,
             xlab = "", ylab = "",
             main = paste0("Cor: ", round(stats::cor(y1, y2, method = "spearman"), 2)))
graphics.off()
stats::cor(y1, y2, method = "spearman")



