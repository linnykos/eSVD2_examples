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
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.005),"external_gene_name"]

candidate_genes <- sort(unique(c(sfari_genes, bulk_de_genes)))

###############

load("../../../out/Writeup12/Writeup12_sns_layer23_esvd3.RData")

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

num_genes <- 300
esvd_selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = F)[1:num_genes]]

load("../../../out/main/sns_layer23_deseq2.RData")
load("../../../out/main/sns_layer23_sctransform.RData")
load("../../../out/main/sns_layer23_mast.RData")

deseq_fdr_val <- stats::p.adjust(deseq2_res$pvalue, method = "BH")
names(deseq_fdr_val) <- rownames(deseq2_res)
deseq_selected_genes <- names(deseq_fdr_val)[order(deseq_fdr_val, decreasing = F)[1:num_genes]]

sct_fdr_val <- stats::p.adjust(de_result$p_val, method = "BH")
names(sct_fdr_val) <- rownames(de_result)
sct_selected_genes <- names(sct_fdr_val)[order(sct_fdr_val, decreasing = F)[1:num_genes]]

mast_fdr_val <- stats::p.adjust(mast_pval_glmer, method = "BH")
names(mast_fdr_val) <- names(mast_pval_glmer)
mast_selected_genes <- names(mast_fdr_val)[order(mast_fdr_val, decreasing = F)[1:num_genes]]

#################

length(intersect(esvd_selected_genes, sfari_genes))
length(intersect(esvd_selected_genes, bulk_de_genes))

length(intersect(deseq_selected_genes, sfari_genes))
length(intersect(deseq_selected_genes, bulk_de_genes))

length(intersect(sct_selected_genes, sfari_genes))
length(intersect(sct_selected_genes, bulk_de_genes))

length(intersect(mast_selected_genes, sfari_genes))
length(intersect(mast_selected_genes, bulk_de_genes))

#######

length(intersect(esvd_selected_genes, deseq_selected_genes))
length(intersect(esvd_selected_genes, sct_selected_genes))
length(intersect(esvd_selected_genes, mast_selected_genes))

####################

input <- c(
  sfari_genes = 2*10^3,
  bulk_de_genes = 5*10^3,
  esvd_selected_genes = 2*10^4,
  deseq_selected_genes = 5*10^4,
  mast_selected_genes = 2*10^5,
  sct_selected_genes = 5*10^5,
  "sfari_genes&esvd_selected_genes" = length(intersect(sfari_genes, esvd_selected_genes)),
  "bulk_de_genes&esvd_selected_genes" = length(intersect(bulk_de_genes, esvd_selected_genes)),
  "sfari_genes&deseq_selected_genes" = length(intersect(sfari_genes, deseq_selected_genes)),
  "bulk_de_genes&deseq_selected_genes" = length(intersect(bulk_de_genes, deseq_selected_genes)),
  "sfari_genes&mast_selected_genes" = length(intersect(sfari_genes, mast_selected_genes)),
  "bulk_de_genes&mast_selected_genes" =  length(intersect(bulk_de_genes, mast_selected_genes)),
  "sfari_genes&sct_selected_genes" = length(intersect(sfari_genes, sct_selected_genes)),
  "bulk_de_genes&sct_selected_genes" =  length(intersect(bulk_de_genes, sct_selected_genes))
)
input_mat <- UpSetR::fromExpression(input)

png("../../../out/fig/main/sns_layer23_upset_againstTruth.png",
    height = 2000, width = 2500,
    units = "px", res = 500)
UpSetR::upset(input_mat,
              nsets = 8,
              intersections = list(list("sfari_genes", "esvd_selected_genes"),
                                   list("bulk_de_genes", "esvd_selected_genes"),
                                   list("sfari_genes", "deseq_selected_genes"),
                                   list("bulk_de_genes", "deseq_selected_genes"),
                                   list("sfari_genes", "mast_selected_genes"),
                                   list("bulk_de_genes", "mast_selected_genes"),
                                   list("sfari_genes", "sct_selected_genes"),
                                   list("bulk_de_genes", "sct_selected_genes")),
              number.angles = 0,
              mb.ratio = c(0.5, 0.5),
              text.scale = 1.5,
              point.size = 2.8,
              line.size = 1)
graphics.off()

###################

input <- c(
  sfari_genes = 2*10^3,
  bulk_de_genes = 5*10^3,
  esvd_selected_genes = 2*10^4,
  deseq_selected_genes = 5*10^4,
  mast_selected_genes = 2*10^5,
  sct_selected_genes = 5*10^5,
  "esvd_selected_genes&deseq_selected_genes" = length(intersect(esvd_selected_genes, deseq_selected_genes)),
  "esvd_selected_genes&mast_selected_genes" = length(intersect(esvd_selected_genes, mast_selected_genes)),
  "esvd_selected_genes&sct_selected_genes" = length(intersect(esvd_selected_genes, sct_selected_genes)),
  "deseq_selected_genes&mast_selected_genes" = length(intersect(deseq_selected_genes, mast_selected_genes)),
  "deseq_selected_genes&sct_selected_genes" = length(intersect(deseq_selected_genes, sct_selected_genes)),
  "esvd_selected_genes&sfari_genes" = 1,
  "esvd_selected_genes&bulk_de_genes" = 1
)
input_mat <- UpSetR::fromExpression(input)

png("../../../out/fig/main/sns_layer23_upset_againstEachOther.png",
    height = 5200, width = 2500,
    units = "px", res = 500)
UpSetR::upset(input_mat,
              nsets = 7,
              intersections = list(list("esvd_selected_genes", "deseq_selected_genes"),
                                   list("esvd_selected_genes", "mast_selected_genes"),
                                   list("esvd_selected_genes", "sct_selected_genes"),
                                   list("deseq_selected_genes", "mast_selected_genes"),
                                   list("deseq_selected_genes", "sct_selected_genes"),
                                   list("deseq_selected_genes", "sfari_genes"),
                                   list("deseq_selected_genes", "bulk_de_genes")),
              number.angles = 0,
              mb.ratio = c(0.5, 0.5),
              text.scale = 1.5,
              point.size = 2.8,
              line.size = 1)
graphics.off()
