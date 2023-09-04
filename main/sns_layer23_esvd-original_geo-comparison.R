rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_esvd.RData")
# load("../../../out/Writeup12/Writeup12_sns_layer23_esvd3.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
selected_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
gene_names <- names(fdr_vec)

############

hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]

hk_genes <- hk_genes[hk_genes %in% gene_names]
sfari_genes <- sfari_genes[sfari_genes %in% gene_names]
bulk_de_genes <- bulk_de_genes[bulk_de_genes %in% gene_names]

load("../../../data/sns_autism/velmeshev_genes.RData")
tmp <- velmeshev_de_gene_df_list[[1]]
tmp <- tmp[which(tmp[,"Cell type"] == "L2/3"),]
de_gene_specific <- tmp[,"Gene name"]

####################################

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

esvd_res@result["GO:0099536",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:0007417",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:0034220",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
esvd_res@result["GO:0048878",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]

# esvd_res@result["GO:0048812",] # neuron projection morphogenesis
# esvd_res@result["GO:0048667",] # cell morphogenesis involved in neuron differentiation

colnames(velmeshev_res@result)
head(velmeshev_res@result[,c("ID", "Description", "pvalue", "GeneRatio")])
velmeshev_res@result["GO:0099536",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
velmeshev_res@result["GO:0061564",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
velmeshev_res@result["GO:0031175",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]
velmeshev_res@result["GO:0050877",c("ID", "Description", "GeneRatio", "BgRatio", "pvalue")]

go_vec <- c("GO:0099536", "GO:0007417", "GO:0034220", "GO:0048878")
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
