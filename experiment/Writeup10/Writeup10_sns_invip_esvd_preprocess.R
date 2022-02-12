rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
print(table(sns@meta.data$celltype))
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "IN-VIP")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)
mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts))
print(dim(mat))

# remove genes with all 0's
gene_sum <- Matrix::colSums(mat)
if(any(gene_sum == 0)) mat <- mat[,-which(gene_sum == 0)]

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
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)
mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts))
print(dim(mat))

# remove genes that are too sparsely observed
case_idx <- which(sns$diagnosis == "ASD")
control_idx <- which(sns$diagnosis == "Control")
threshold <- 0.05
gene_bool <- sapply(1:ncol(mat), function(j){
  case_percentage <- length(which(mat[case_idx,j] > 0))/length(case_idx)
  control_percentage <- length(which(mat[control_idx,j] > 0))/length(control_idx)

  max(case_percentage, control_percentage) > threshold
})
mat <- mat[,which(gene_bool)]
print(dim(mat))

# remove genes with too high counts
gene_total <- matrixStats::colSums2(log1p(mat))
gene_threshold <- log1p(10)*nrow(mat)
gene_keep <- colnames(mat)[which(gene_total <= gene_threshold)]
mat <- mat[,gene_keep]
print(dim(mat))

sns <- subset(sns, features = colnames(mat))
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::FindVariableFeatures(sns,
                                    selection.method = "vst",
                                    nfeatures = 5000)

load("../../../../data/sns_autism/velmeshev_genes.RData")
de_genes1 <- velmeshev_marker_gene_df[,"Gene name"]
de_genes2 <- unlist(lapply(velmeshev_de_gene_df_list, function(de_mat){
  idx <- ifelse("Gene name" %in% colnames(de_mat), "Gene name", "HGNC Symbol")
  de_mat[,idx]
}))
de_genes <- sort(unique(c(de_genes1, de_genes2)))
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

gene_keep <- unique(c(de_genes,  hk_genes, sfari_genes, cycling_genes, Seurat::VariableFeatures(sns)))
gene_keep <- rownames(sns[["RNA"]]@counts)[which(rownames(sns[["RNA"]]@counts) %in% gene_keep)]
sns[["RNA"]]@var.features <- gene_keep
print(length(gene_keep))

################

# keep cells with non-trivial variance
cell_variance <- matrixStats::rowSds(mat)
keep_vec <- rep(1, nrow(mat))
keep_vec[cell_variance < 0.1] <- 0
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)

##############

categorical_var <- c("individual", "diagnosis", "region", "sex", "Seqbatch") #, "individual")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
n <- nrow(mat)
metadata <- sns@meta.data
covariates <- as.matrix(metadata[,numerical_var])
covariates <- cbind(1, log(matrixStats::rowSums2(mat)), covariates)
colnames(covariates)[1:2] <- c("Intercept", "Log_UMI")

for(variable in categorical_var){
  vec <- metadata[,variable]
  uniq_level <- unique(vec)
  for(i in uniq_level[-1]){
    tmp <- rep(0, n)
    tmp[which(vec == i)] <- 1

    var_name <- paste0(variable, "_", i)
    covariates <- cbind(covariates, tmp)
    colnames(covariates)[ncol(covariates)] <- var_name
  }
}

# regress all variables against intercept + Log_UMI + diagnosis_ASD
keep_idx <- which(colnames(covariates) %in% c("Intercept", "Log_UMI", "diagnosis_ASD"))
other_idx <- which(colnames(covariates) %in% c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt",
                                               "nFeature_RNA", "region_PFC", "sex_F", "Seqbatch_SB2", "Seqbatch_SB1"))
for(j in other_idx){
  df_tmp <- data.frame(covariates[,j], covariates[,keep_idx])
  colnames(df_tmp)[1] <- "tmp"
  lm_fit <- stats::lm("tmp ~ . ", data = df_tmp)
  vec_tmp <- stats::residuals(lm_fit)
  covariates[,j] <- vec_tmp
}

cols_regress_out <- grep("individual", colnames(covariates))
covariates_new <- covariates[,-cols_regress_out,drop = F]
for(i in 1:length(cols_regress_out)){
  df_tmp <- data.frame(covariates[,cols_regress_out[i]], covariates_new)
  colnames(df_tmp)[1] <- "tmp"
  lm_fit <- stats::lm("tmp ~ . - 1", data = df_tmp)
  vec_tmp <- stats::residuals(lm_fit)
  if(sum(abs(vec_tmp)) < 1e-6) break()
  covariates_new <- cbind(covariates_new, vec_tmp)
  colnames(covariates_new)[ncol(covariates_new)] <- paste0("individual_",i)
}
covariates <- covariates_new

ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("covariates", "sns", "date_of_run", "session_info")]
rm(list = ls_vec)
save(sns, covariates, date_of_run, session_info,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_processed.RData")
