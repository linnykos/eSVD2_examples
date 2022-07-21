rm(list=ls())
library(Seurat)

load("../../../out/main/regevEpi_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

n <- ncol(regevEpi)
keep_vec <- rep(0, n)
keep_vec[regevEpi$Sample_Location == "Epi"] <- 1
regevEpi$keep <- keep_vec
regevEpi <- subset(regevEpi, keep == 1)

regevEpi_original <- regevEpi
celltype_labels <- c("Cycling TA", "Enterocyte Progenitors", "TA 1", "TA 2")
celltype_names <- c("cyclingta", "entprog", "ta1", "ta2")

for(ii in 1:length(celltype_labels)){
  regevEpi <- regevEpi_original
  celltype_label <- celltype_labels[ii]
  celltype_name <- celltype_names[ii]
  print(paste0("Working on ", celltype_name))

  keep_vec <- rep(0, ncol(regevEpi))
  keep_vec[which(regevEpi$Celltype == celltype_label)] <- 1
  regevEpi[["keep"]] <- keep_vec
  regevEpi <- subset(regevEpi, keep == 1)

  # keep only healthy subjects or subjects with enough inflamed/non-inflamed cells
  tab_mat <- table(regevEpi$Sample_Health, regevEpi$Subject)
  subj_healthy <- colnames(tab_mat)[which(tab_mat[1,]!=0)]
  subj_disease <- colnames(tab_mat)[intersect(which(tab_mat[2,]>=50), which(tab_mat[3,]>=50))]
  keep_vec <- rep(0, ncol(regevEpi))
  keep_vec[which(regevEpi$Subject %in% c(subj_healthy, subj_disease))] <- 1
  regevEpi[["keep"]] <- keep_vec
  regevEpi <- subset(regevEpi, keep == 1)
  mat <- as.matrix(Matrix::t(regevEpi[["RNA"]]@counts))
  print(dim(mat))

  # remove genes with all 0's
  gene_sum <- Matrix::colSums(mat)
  if(any(gene_sum == 0)) mat <- mat[,-which(gene_sum == 0)]
  print(dim(mat))

  # remove cells with too high counts
  n <- nrow(mat)
  rank_cutoff <- 5
  set.seed(10)
  high_expression_mat <- sapply(1:ncol(mat), function(j){
    order(mat[,j] + stats::runif(n, min = 0, max = .1), decreasing = T)[1:rank_cutoff]
  })
  upper_limit <- 4*(rank_cutoff/nrow(mat)*ncol(mat)) # inflation factor of 4
  cell_high_expression <- table(as.numeric(high_expression_mat))
  remove_cell_idx <- names(cell_high_expression)[cell_high_expression >= upper_limit]
  keep_vec <- rep(1, nrow(mat))
  keep_vec[as.numeric(remove_cell_idx)] <- 0
  regevEpi[["keep"]] <- keep_vec
  regevEpi <- subset(regevEpi, keep == 1)
  mat <- as.matrix(Matrix::t(regevEpi[["RNA"]]@counts))
  print(dim(mat))

  # remove genes that are too sparsely observed
  case_idx <- which(regevEpi$Subject_Disease == "Colitis")
  control_idx <- which(regevEpi$Subject_Disease == "HC")
  threshold <- 0.01
  gene_bool <- sapply(1:ncol(mat), function(j){
    case_percentage <- length(which(mat[case_idx,j] > 0))/length(case_idx)
    control_percentage <- length(which(mat[control_idx,j] > 0))/length(control_idx)

    max(case_percentage, control_percentage) > threshold
  })
  mat <- mat[,which(gene_bool)]
  print(dim(mat))

  regevEpi <- subset(regevEpi, features = colnames(mat))
  Seurat::DefaultAssay(regevEpi) <- "RNA"
  regevEpi <- Seurat::NormalizeData(regevEpi)
  regevEpi <- Seurat::FindVariableFeatures(regevEpi,
                                           selection.method = "vst",
                                           nfeatures = 5000)

  sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                            sheet = "Epithelial (Non-Inflamed vs. He"))
  sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                            sheet = "Epithelial (Inflamed vs. Health"))
  sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                            sheet = "Epithelial (Inflamed vs. Non-In"))
  de_genes <- sort(unique(c(sheet1$gene, sheet2$gene, sheet3$gene)))
  hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
  cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

  gene_keep <- unique(c(de_genes,  hk_genes, cycling_genes, Seurat::VariableFeatures(regevEpi)))
  gene_keep <- rownames(regevEpi[["RNA"]]@counts)[which(rownames(regevEpi[["RNA"]]@counts) %in% gene_keep)]
  regevEpi[["RNA"]]@var.features <- gene_keep

  ################

  # keep cells with non-trivial variance
  mat <- as.matrix(Matrix::t(regevEpi[["RNA"]]@counts))
  cell_variance <- matrixStats::rowSds(mat)
  keep_vec <- rep(1, nrow(mat))
  keep_vec[cell_variance < 0.1] <- 0
  regevEpi[["keep"]] <- keep_vec
  regevEpi <- subset(regevEpi, keep == 1)
  print(dim(regevEpi))

  save(regevEpi, date_of_run, session_info,
       file = paste0("../../../out/main/regevEpi_", celltype_name,"_preprocessed.RData"))
}


