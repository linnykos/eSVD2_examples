rm(list=ls())
library(Seurat)
library(eSVD2)
library(SummarizedExperiment)
library(DESeq2)

load("../../../out/main/adams_T_sctransform.RData")
sctransform_adams <- de_result
load("../../../out/main/habermann_T_sctransform.RData")
sctransform_habermann <- de_result

load("../../../out/main/adams_T_deseq2.RData")
adams_deseq2 <- deseq2_res
load("../../../out/main/habermann_T_deseq2.RData")
habermann_deseq2 <- deseq2_res

load("../../../out/main/habermann_T_esvd.RData")
eSVD_obj_habermann <- eSVD_obj
eSVD_obj_habermann <- eSVD2:::compute_pvalue(eSVD_obj_habermann)

load("../../../out/main/adams_T_esvd.RData")
eSVD_obj_adams <- eSVD_obj
eSVD_obj_adams <- eSVD2:::compute_pvalue(eSVD_obj_adams)

gene_names <- names(eSVD_obj_adams$teststat_vec)
all(gene_names %in% names(eSVD_obj_habermann$teststat_vec))

############

df_mat <- read.csv("~/project/eSVD/data/GSE136831_adams_lung/aba1983_Data_S8.txt",
                   sep = "\t")
adams_df_genes <- df_mat$gene[which(df_mat$cellType == "T")]
df_mat <- read.csv("~/project/eSVD/data/GSE135893_habermann_lung/Table_S4_DEG_analysis/Disease_vs_Control/T_Cells_disease_vs_control_.csv",
                   sep = ",")
habermann_df_genes <- df_mat$X

de_genes <- unique(c(adams_df_genes, habermann_df_genes))
hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
hk_genes <- setdiff(hk_genes, de_genes)

adams_df_genes <- intersect(adams_df_genes, gene_names)
habermann_df_genes <- intersect(habermann_df_genes, gene_names)

target_length <- length(unique(c(adams_df_genes, habermann_df_genes)))

#####################

logpvalue_vec <- eSVD_obj_adams$pvalue_list$log10pvalue
idx_adams <- order(logpvalue_vec, decreasing = T)[1:target_length]
eSVD_adams_de <- sort(names(logpvalue_vec)[idx_adams])

logpvalue_vec <- eSVD_obj_habermann$pvalue_list$log10pvalue
idx_habermann <- order(logpvalue_vec, decreasing = T)[1:target_length]
eSVD_habermann_de <- sort(names(logpvalue_vec)[idx_habermann])

###############################

deseq_adams_degenes <- rownames(adams_deseq2)[order(adams_deseq2[,"pvalue"], decreasing = F)[1:target_length]]
deseq_habermann_degenes <- rownames(habermann_deseq2)[order(habermann_deseq2[,"pvalue"], decreasing = F)[1:target_length]]

sctransform_adams_degenes <- rownames(sctransform_adams)[order(sctransform_adams[,"p_val"], decreasing = F)[1:target_length]]
sctransform_habermann_degenes <- rownames(sctransform_habermann)[order(sctransform_habermann[,"p_val"], decreasing = F)[1:target_length]]

de_vec <- unique(c(adams_df_genes, habermann_df_genes))

length(de_vec)
length(adams_df_genes)
length(habermann_df_genes)
length(intersect(adams_df_genes, habermann_df_genes))

##

length(eSVD_adams_de)
length(intersect(eSVD_adams_de, habermann_df_genes))
length(intersect(eSVD_adams_de, adams_df_genes))
length(intersect(eSVD_adams_de, de_vec))
length(intersect(eSVD_adams_de, hk_genes))

length(eSVD_habermann_de)
length(intersect(eSVD_habermann_de, adams_df_genes))
length(intersect(eSVD_habermann_de, habermann_df_genes))
length(intersect(eSVD_habermann_de, de_vec))
length(intersect(eSVD_habermann_de, hk_genes))

length(intersect(eSVD_adams_de, eSVD_habermann_de))

##
length(sctransform_adams_degenes)
length(intersect(sctransform_adams_degenes, habermann_df_genes))
length(intersect(sctransform_adams_degenes, adams_df_genes))
length(intersect(sctransform_adams_degenes, de_vec))
length(intersect(sctransform_adams_degenes, hk_genes))

length(sctransform_habermann_degenes)
length(intersect(sctransform_habermann_degenes, adams_df_genes))
length(intersect(sctransform_habermann_degenes, habermann_df_genes))
length(intersect(sctransform_habermann_degenes, de_vec))
length(intersect(sctransform_habermann_degenes, hk_genes))

length(intersect(sctransform_adams_degenes, sctransform_habermann_degenes))

##
length(deseq_adams_degenes)
length(intersect(deseq_adams_degenes, habermann_df_genes))
length(intersect(deseq_adams_degenes, adams_df_genes))
length(intersect(deseq_adams_degenes, de_vec))
length(intersect(deseq_adams_degenes, hk_genes))

length(deseq_habermann_degenes)
length(intersect(deseq_habermann_degenes, adams_df_genes))
length(intersect(deseq_habermann_degenes, habermann_df_genes))
length(intersect(deseq_habermann_degenes, de_vec))
length(intersect(deseq_habermann_degenes, hk_genes))

length(intersect(deseq_adams_degenes, deseq_habermann_degenes))

##################

input <- c(
  h.reported = 200,
  a.reported = 500,
  h.esvd = 2000,
  a.esvd = 5000,
  h.deseq = 20000,
  a.deseq = 50000,
  "h.reported&a.esvd" = length(intersect(habermann_df_genes, eSVD_adams_de)),
  "a.reported&h.esvd" = length(intersect(adams_df_genes, eSVD_habermann_de)),
  "a.esvd&h.esvd" = length(intersect(eSVD_adams_de, eSVD_habermann_de)),
  "h.reported&a.deseq" = length(intersect(habermann_df_genes, deseq_adams_degenes)),
  "a.reported&h.deseq" = length(intersect(adams_df_genes, deseq_habermann_degenes)),
  "h.deseq&a.deseq" =  length(intersect(deseq_adams_degenes, deseq_habermann_degenes))
)
input_mat <- UpSetR::fromExpression(input)
# input_mat <- input_mat[-c(1:(length(adams_df_genes)+length(habermann_df_genes)+length(eSVD_adams_de)+length(eSVD_habermann_de)+length(sctransform_adams_degenes)+length(sctransform_habermann_degenes))),]

png("../../../out/fig/main/adams_habermann_T_upset.png",
    height = 2000, width = 2000,
    units = "px", res = 500)
UpSetR::upset(input_mat,
              nsets = 6,
              intersections = list(list("h.reported", "a.esvd"),
                                   list("a.reported", "h.esvd"),
                                   list("a.esvd", "h.esvd"),
                                   list("h.reported", "a.deseq"),
                                   list("a.reported", "h.deseq"),
                                   list("h.deseq", "a.deseq")),
              number.angles = 0,
              mb.ratio = c(0.5, 0.5),
              text.scale = 1.5,
              point.size = 2.8,
              line.size = 1)
graphics.off()

#######################

input <- c(
  h.reported = length(habermann_df_genes),
  a.reported = length(adams_df_genes),
  h.esvd = length(eSVD_habermann_de),
  a.esvd = length(eSVD_adams_de),
  h.deseq = length(deseq_habermann_degenes),
  a.deseq = length(deseq_adams_degenes),
  h.sct = length(sctransform_habermann_degenes),
  a.sct = length(sctransform_adams_degenes),
  "h.reported&a.reported" = length(intersect(habermann_df_genes, adams_df_genes)),
  "h.reported&a.esvd" = length(intersect(habermann_df_genes, eSVD_adams_de)),
  "h.reported&h.esvd" = length(intersect(habermann_df_genes, eSVD_habermann_de)),
  "a.reported&h.esvd" = length(intersect(adams_df_genes, eSVD_habermann_de)),
  "a.reported&a.esvd" = length(intersect(adams_df_genes, eSVD_adams_de)),
  "a.esvd&h.esvd" = length(intersect(eSVD_adams_de, eSVD_habermann_de)),
  "h.reported&a.deseq" = length(intersect(habermann_df_genes, deseq_adams_degenes)),
  "h.reported&h.deseq" = length(intersect(habermann_df_genes, deseq_habermann_degenes)),
  "a.reported&h.deseq" = length(intersect(adams_df_genes, deseq_habermann_degenes)),
  "a.reported&a.deseq" = length(intersect(adams_df_genes, deseq_adams_degenes)),
  "h.deseq&a.deseq" =  length(intersect(deseq_habermann_degenes, deseq_adams_degenes)),
  "h.reported&a.sct" = length(intersect(habermann_df_genes, sctransform_adams_degenes)),
  "h.reported&h.sct" = length(intersect(habermann_df_genes, sctransform_habermann_degenes)),
  "a.reported&h.sct" = length(intersect(adams_df_genes, sctransform_habermann_degenes)),
  "a.reported&a.sct" = length(intersect(adams_df_genes, sctransform_adams_degenes)),
  "h.sct&a.sct" =  length(intersect(sctransform_habermann_degenes, sctransform_adams_degenes))
)
input_mat <- UpSetR::fromExpression(input)
# input_mat <- input_mat[-c(1:(length(adams_df_genes)+length(habermann_df_genes)+length(eSVD_adams_de)+length(eSVD_habermann_de)+length(sctransform_adams_degenes)+length(sctransform_habermann_degenes))),]

png("../../../out/fig/main/adams_habermann_T_upset-full.png",
    height = 2000, width = 5000,
    units = "px", res = 500)
UpSetR::upset(input_mat,
              nsets = 8,
              nintersects = 24,
              empty.intersections = T,
              number.angles = 0,
              mb.ratio = c(0.5, 0.5),
              text.scale = 1.5,
              point.size = 2.8,
              line.size = 1)
graphics.off()

