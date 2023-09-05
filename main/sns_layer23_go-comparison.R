rm(list=ls())
library(Seurat)
library(eSVD2)
library(clusterProfiler)
library(org.Hs.eg.db)

load("../../../out/main/sns_layer23_esvd.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

#################

fdr_vec <- eSVD_obj$pvalue_list$fdr_vec
esvd_genes <- names(fdr_vec)[which(fdr_vec <= 0.05)]
gene_names <- names(fdr_vec)

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

all_terms <- sort(unique(c(esvd_terms, velmeshev_terms)))
all_terms <- intersect(intersect(all_terms, zz$ID), zz2$ID)

esvd_better <- all_terms[zz[all_terms,"qvalue"] < zz2[all_terms,"qvalue"]]
velmeshev_better <- all_terms[zz[all_terms,"qvalue"] > zz2[all_terms,"qvalue"]]


