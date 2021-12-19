rm(list=ls())
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/sns_autism/sns_formatted.RData")
head(sns@meta.data)
keep_vec <- rep(0, ncol(sns))
keep_vec[which(sns@meta.data$celltype == "L4")] <- 1
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
keep_vec[which(cell_violation <= 100)] <- 1
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

# these are from Supplement 4 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7678724/
# in the "ASD_DEGs" tab
de_genes <- c("CIRBP", "FSTL5", "TCF25", "NEGR1", "SAT2",
              "FOXP2", "NCAM2", "MLLT3", "FAM153B", "ATPIF1",
              "SATB2", "PTMAP5", "TCEAL2", "OLFM3", "CDH2",
              "CDC37", "NRN1", "COX5B", "TBRG1", "PTMAP2",
              "PPFIA2", "NLGN4X", "WDR60", "AKAP9", "BEX1",
              "GABRB1", "C16orf45", "CHCHD2", "ANO3", "KIAA1456",
              "GOLGA8B", "RP11-405A12.2", "RPS15", "SYNPO2", "HINT1",
              "TCF4", "AC007969.5", "MEF2C-AS1", "CABP1", "SRRM4",
              "AC105402.4", "ZNF208", "NKAIN3", "USP34", "CXXC4",
              "EPHA7", "MARCH1", "ZEB2", "RP11-166D19.1", "SH3D19",
              "SCN3B", "LRRC4C", "SORBS2", "EIF1", "BDP1",
              "RAPGEFL1", "SYN2", "HDAC4", "PTPN13", "GABRG1",
              "CSRNP3", "SYNPR", "THRB", "ADGRB3", "RAB3C",
              "SOX5", "ST6GAL2", "ADGRL2", "RP11-707M1.1", "PLPPR4",
              "PUM1", "RASGRF1", "MPP6", "OR2L13", "RP11-444D3.1",
              "SYT17", "NRXN1", "GUCY1A2", "FAM153A", "DGKZ",
              "MIR137HG", "LRRN3", "EML6", "ZDHHC21", "RORB",
              "CH17-472G23.1", "UNC5C", "PLCH1", "KCNH5", "SPIN1",
              "TMEFF2", "GABRA5", "GRIA4", "RP11-586K2.1", "RMST",
              "AK5", "PDE4DIP", "XIST", "GALNTL6", "KCND2",
              "PDZRN4", "ZFYVE28", "MGAT4C", "BRINP3", "NGFRAP1",
              "LMO3", "FUT9", "RBFOX3")
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
de_genes <- de_genes[de_genes %in% colnames(mat)]

ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("covariates", "mat", "de_genes", "sns", "date_of_run", "session_info", "metadata")]
rm(list = ls_vec)
save.image("../../../../out/writeup8e/writeup8e_sns_layer4_processed.RData")
save(mat, covariates, de_genes, date_of_run, session_info, metadata,
     file = "../../../../out/writeup8e/writeup8e_sns_layer4_processed_onlymat.RData")

