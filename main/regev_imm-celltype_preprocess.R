rm(list=ls())
library(Seurat)

load("../../../out/main/regevImm_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

n <- ncol(regevImm)
keep_vec <- rep(0, n)
keep_vec[regevImm$Sample_Location == "LP"] <- 1
regevImm$keep <- keep_vec
regevImm <- subset(regevImm, keep == 1)

regevImm_original <- regevImm
celltype_labels <- c("Plasma", "Macrophages")
celltype_names <- c("plasma", "macro")

for(ii in 1:length(celltype_labels)){
  regevImm <- regevImm_original
  celltype_label <- celltype_labels[ii]
  celltype_name <- celltype_names[ii]
  print(paste0("Working on ", celltype_name))

  keep_vec <- rep(0, ncol(regevImm))
  keep_vec[which(regevImm$Celltype == celltype_label)] <- 1
  regevImm[["keep"]] <- keep_vec
  regevImm <- subset(regevImm, keep == 1)

  # keep only healthy subjects or subjects with enough inflamed/non-inflamed cells
  tab_mat <- table(regevImm$Sample_Health, regevImm$Subject)
  subj_healthy <- colnames(tab_mat)[which(tab_mat[1,]!=0)]
  subj_disease <- colnames(tab_mat)[intersect(which(tab_mat[2,]>=50), which(tab_mat[3,]>=50))]
  keep_vec <- rep(0, ncol(regevImm))
  keep_vec[which(regevImm$Subject %in% c(subj_healthy, subj_disease))] <- 1
  regevImm[["keep"]] <- keep_vec
  regevImm <- subset(regevImm, keep == 1)
  mat <- as.matrix(Matrix::t(regevImm[["RNA"]]@counts))
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
  regevImm[["keep"]] <- keep_vec
  regevImm <- subset(regevImm, keep == 1)
  mat <- as.matrix(Matrix::t(regevImm[["RNA"]]@counts))
  print(dim(mat))

  # remove genes that are too sparsely observed
  case_idx <- which(regevImm$Subject_Disease == "Colitis")
  control_idx <- which(regevImm$Subject_Disease == "HC")
  threshold <- 0.01
  gene_bool <- sapply(1:ncol(mat), function(j){
    case_percentage <- length(which(mat[case_idx,j] > 0))/length(case_idx)
    control_percentage <- length(which(mat[control_idx,j] > 0))/length(control_idx)

    max(case_percentage, control_percentage) > threshold
  })
  mat <- mat[,which(gene_bool)]
  print(dim(mat))

  regevImm <- subset(regevImm, features = colnames(mat))
  Seurat::DefaultAssay(regevImm) <- "RNA"
  regevImm <- Seurat::NormalizeData(regevImm)
  regevImm <- Seurat::FindVariableFeatures(regevImm,
                                           selection.method = "vst",
                                           nfeatures = 5000)

  if(celltype_label == "Plasma"){
    sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                              sheet = "Adaptive (Non-Inflamed vs. Heal"))
    sheet1 <- sheet1[sheet1$ident == "Plasma",]
    sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                              sheet = "Adaptive (Inflamed vs. Healthy)"))
    sheet2 <- sheet2[sheet2$ident == "Plasma",]
    sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                              sheet = "Adaptive (Inflamed vs. Non-Infl"))
    sheet3 <- sheet3[sheet3$ident == "Plasma",]

  } else if(celltype_label == "Macrophages"){
    sheet1 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                              sheet = "Innate (Non-Inflamed vs. Health"))
    sheet1 <- sheet1[sheet1$ident == "Macrophages",]
    sheet2 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                              sheet = "Innate (Inflamed vs. Healthy)"))
    sheet2 <- sheet2[sheet2$ident == "Macrophages",]
    sheet3 <- as.data.frame(readxl::read_xlsx("~/nzhanglab/data/SCP259_regev_colitis/NIHMS1532849-supplement-11.xlsx",
                                              sheet = "Innate (Inflamed vs. Non-Inflam"))
    sheet3 <- sheet3[sheet3$ident == "Macrophages",]
  }

  de_genes <- sort(unique(c(sheet1$gene, sheet2$gene, sheet3$gene)))
  hk_genes <- read.csv("../../../data/housekeeping/housekeeping.txt", header = F)[,1]
  cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

  gene_keep <- unique(c(de_genes,  hk_genes, cycling_genes, Seurat::VariableFeatures(regevImm)))
  gene_keep <- rownames(regevImm[["RNA"]]@counts)[which(rownames(regevImm[["RNA"]]@counts) %in% gene_keep)]
  regevImm[["RNA"]]@var.features <- gene_keep

  ################

  # keep cells with non-trivial variance
  mat <- as.matrix(Matrix::t(regevImm[["RNA"]]@counts))
  cell_variance <- matrixStats::rowSds(mat)
  keep_vec <- rep(1, nrow(mat))
  keep_vec[cell_variance < 0.1] <- 0
  regevImm[["keep"]] <- keep_vec
  regevImm <- subset(regevImm, keep == 1)
  print(dim(regevImm))

  save(regevImm, date_of_run, session_info,
       file = paste0("../../../out/main/regevImm_", celltype_name,"_preprocessed.RData"))
}


