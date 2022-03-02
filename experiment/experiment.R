rm(list=ls())
rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed.RData")

library(Seurat)
set.seed(10)

mat <- as.matrix(Matrix::t(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,]))
covariate_dat <- sns@meta.data[,c("percent.mt", "individual", "region", "age", "sex",
                                  "RNA.Integrity.Number", "post.mortem.hours",
                                  "diagnosis", "Seqbatch")]
covariate_df <- data.frame(covariate_dat)
covariate_df[,"individual"] <- as.factor(covariate_df[,"individual"])
covariate_df[,"region"] <- as.factor(covariate_df[,"region"])
covariate_df[,"diagnosis"] <- factor(covariate_df[,"diagnosis"], levels = c("Control", "ASD"))
covariate_df[,"sex"] <- as.factor(covariate_df[,"sex"])
covariate_df[,"Seqbatch"] <- as.factor(covariate_df[,"Seqbatch"])
covariates <- eSVD2:::format_covariates(dat = mat,
                                        covariate_df = covariate_df,
                                        mixed_effect_variables = c("individual", "Seqbatch"))

load("../../../../out/Writeup11/Writeup11_sns_invip_esvd_coef_tmp.RData")
case_control_variable = "diagnosis_ASD"
offset_variables = "Log_UMI"
k = 10
verbose = 1
dat = mat

n <- nrow(mat)
if(verbose >= 1) print("Step 1b: Cleaning up coefficients")
covariates <- cbind(rep(1, n), covariates)
colnames(covariates)[1] <- "Intercept"
print(dim(covariates))
print(dim(b_mat))
col_idx <- sapply(colnames(covariates), function(i){which(colnames(b_mat) == i)})
b_mat <- b_mat[,as.numeric(col_idx)]

if(verbose >= 1) print("Step 2: Computing residuals")
tmp <- eSVD2:::.initialize_residuals(b_mat = b_mat,
                                     covariates = covariates,
                                     dat = dat,
                                     k = k)

esvd_init <- structure(list(x_mat = tmp$x_mat, y_mat = tmp$y_mat,
                            b_mat = b_mat,
                            covariates = covariates,
                            nuisance_param_vec = rep(0, ncol(dat))),
                       class = "eSVD")

save(mat, esvd_init,
     file = "../../../../out/Writeup11/Writeup11_sns_invip_esvd_coef.RData")

load("../../../../data/sns_autism/velmeshev_genes.RData")
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
hk_genes <- read.csv("../../../../data/housekeeping/housekeeping.txt", header = F)[,1]
sfari_genes <- read.csv("../../../../data/SFARI/SFARI-Gene_genes_09-02-2021release_01-06-2022export.csv", header = T)[,2]
cycling_genes <- c(cc.genes$s.genes, cc.genes$g2m.genes)

hk_idx <- which(colnames(mat) %in% c(hk_genes, cycling_genes))
de_idx <- which(colnames(mat) %in% de_gene_specific)
other_idx <- which(colnames(mat) %in% c(sfari_genes, de_genes))

zz <- esvd_init$b_mat[,"diagnosis_ASD"]
# zz <- esvd_res$b_mat[,"diagnosis_ASD"]
# zz <- esvd_res_full$b_mat[,"diagnosis_ASD"]
length(which(zz[hk_idx] == 0))/length(hk_idx)
tmp <- zz[hk_idx]; round(quantile(abs(tmp[tmp!=0])),3)
length(which(zz[de_idx] == 0))/length(de_idx)
tmp <- zz[de_idx]; round(quantile(abs(tmp[tmp!=0])),3)
length(which(zz[other_idx] == 0))/length(other_idx)
tmp <- zz[other_idx]; round(quantile(abs(tmp[tmp!=0])),3)
length(which(zz == 0))/length(zz)
round(quantile(abs(zz[zz!=0])),3)

##################################3




