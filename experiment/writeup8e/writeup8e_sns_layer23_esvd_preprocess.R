rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "L2/3")] <- 1
sns[["keep"]] <- keep_vec
sns <- subset(sns, keep == 1)

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts))
print(dim(mat))

# remove cells with too high counts
bool_mat <- sapply(1:ncol(mat), function(j){
  vec <- rep(0, nrow(mat))
  vec[order(mat[,j], decreasing = T)[1:5]] <- 1
  vec
})
cell_violation <- matrixStats::rowSums2(bool_mat)
quantile(cell_violation)
quantile(cell_violation, probs = seq(0,0.1,length.out=11))
quantile(cell_violation, probs = seq(0.9,1,length.out=11))
keep_vec <- rep(0, nrow(mat))
keep_vec[which(cell_violation <= 30)] <- 1
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
table(gene_bool)
mat <- mat[,which(gene_bool)]
print(dim(mat))

# remove genes with too high counts
gene_total <- matrixStats::colSums2(log1p(mat))
quantile(gene_total, probs = seq(0,0.1,length.out=11))
quantile(gene_total, probs = seq(0.9,1,length.out=11))
gene_keep <- colnames(mat)[which(gene_total <= 30000)]
mat <- mat[,gene_keep]
print(dim(mat))

sns <- subset(sns, features = colnames(mat))
Seurat::DefaultAssay(sns) <- "RNA"
sns <- Seurat::NormalizeData(sns)
sns <- Seurat::FindVariableFeatures(sns,
                                    selection.method = "vst",
                                    nfeatures = 5000)

de_genes <- c("TTF2", "MX2", "ASCC1", "GLRA3", "CIRBP",
              "SAT2", "QTRT1", "CDH2", "LUC7L", "TCF25",
              "SSBP2", "WDR60", "CABP1", "FBLN7", "CDC14B",
              "GPM6A", "IGFBP5", "FAM153B", "GUCY1A2", "RAB3C",
              "SSX2IP", "HS6ST3", "TENM3", "DACH1", "PLA2G4C",
              "TOX3", "SPAG16", "FAM171B", "GALNTL6", "NUMB",
              "CAPZB", "DDRGK1", "RMST", "SUGP2", "FAM49A",
              "KCNH7", "BRINP3", "GABRB1", "GOLGA8B", "OR2L13",
              "IMMP2L", "ARPP19", "VWA8", "RPS15", "DPYSL2",
              "RFX3", "RSRP1", "NFIA", "SNRNP70", "SYN2",
              "SPIN1", "PLPPR4", "SYNPR", "SLC22A10", "LINC01378",
              "RP11-577H5.5", "GABRG2", "MIR99AHG", "PPP3CA",
              "MIR137HG", "TBRG1", "GGT7", "NLGN1",
              "GNG7", "FZD3", "LRRTM3", "CPE", "KCNJ3",
              "AQP4-AS1", "TRAF3", "PKIA", "MGAT4C", "HNRNPDL",
              "SLITRK4", "BMPR1B", "AHI1", "CDH9", "RAPGEFL1",
              "RPL34P18", "LINC00657", "COL26A1", "CNTN3", "FRMD6",
              "RP11-30J20.1", "SLITRK5", "SLC39A10", "STX1A", "RPLP2",
              "MAP2", "CES4A", "NEGR1", "SORBS1", "COL24A1",
              "VSTM2L", "ERBB4", "STARD4-AS1", "MAPK1", "HSP90AA1",
              "RPL34", "CNTNAP2", "EIF1", "OLFM3", "GRID2",
              "CHL1", "RAP1GAP", "CAMK2N1", "SERINC1", "RGS12",
              "ATP1B1")
gene_keep <- unique(c(de_genes, Seurat::VariableFeatures(sns)))
gene_keep <- rownames(sns[["RNA"]]@counts)[which(rownames(sns[["RNA"]]@counts) %in% gene_keep)]
length(intersect(gene_keep, de_genes))
length(de_genes)
mat <- mat[,gene_keep]

################

# keep cells with non-trivial variance
cell_variance <- matrixStats::rowSds(mat)
quantile(cell_variance, probs = seq(0,0.1,length.out=11))

##############

categorical_var <- c("individual", "diagnosis", "region", "sex", "Seqbatch") #, "individual")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
n <- nrow(mat)
metadata <- sns@meta.data
covariates <- as.matrix(metadata[,numerical_var])
covariates <- cbind(1, log(matrixStats::rowMeans2(mat)), covariates)
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
de_genes <- de_genes[de_genes %in% colnames(mat)]

ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("covariates", "mat", "de_genes", "sns", "date_of_run", "session_info", "metadata")]
rm(list = ls_vec)
save.image("../../../../out/writeup8e/writeup8e_sns_layer23_processed.RData")
save(mat, covariates, de_genes, date_of_run, session_info, metadata,
     file = "../../../../out/writeup8e/writeup8e_sns_layer23_processed_onlymat.RData")
