rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

###############

load("../../../out/main/sns_layer23_esvd.RData")
eSVD_obj <- eSVD2:::compute_pvalue(eSVD_obj)
logpvalue_vec <- eSVD_obj$pvalue_list$log10pvalue
num_genes <- 100
esvd_selected_genes <- names(logpvalue_vec)[order(logpvalue_vec, decreasing = T)[1:num_genes]]

##############

sfari_genes <- read.csv("../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
deg_df <- readxl::read_xlsx(
  path = "../../../data/bulkRNA-DEG-autism/SupplementaryTable3.xlsx",
  sheet = "DEGene_Statistics"
)
deg_df <- as.data.frame(deg_df)
bulk_de_genes <- deg_df[which(deg_df[,"WholeCortex_ASD_FDR"]<=0.05),"external_gene_name"]

candidate_genes <- sort(unique(c(sfari_genes, bulk_de_genes)))
gene_names <- names(logpvalue_vec)
candidate_genes <- candidate_genes[candidate_genes %in% gene_names]

################

load("../../../out/main/sns_layer23_deseq2.RData")
load("../../../out/main/sns_layer23_sctransform.RData")
load("../../../out/main/sns_layer23_mast.RData")

deseq_log10 <- -log10(deseq2_res$pvalue)
names(deseq_log10) <- rownames(deseq2_res)
deseq_selected_genes <- names(deseq_log10)[order(deseq_log10, decreasing = T)[1:num_genes]]

sct_log10 <- -log10(de_result$p_val)
names(sct_log10) <- rownames(de_result)
sct_selected_genes <- names(sct_log10)[order(sct_log10, decreasing = T)[1:num_genes]]

mast_log10 <- -log10(mast_pval_glmer)
names(mast_log10) <- names(mast_pval_glmer)
mast_selected_genes <- names(mast_log10)[order(mast_log10, decreasing = T)[1:num_genes]]

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
    height = 2000, width = 2250,
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
