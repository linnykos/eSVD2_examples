rm(list=ls())
library(Seurat)
library(SummarizedExperiment)
library(MAST)
library(lme4)

load("../../../out/simulation/simulation_1.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(seurat_obj[["RNA"]]@counts)
rds <- colSums(mat)
med_rds <- median(rds)
mat <- t(t(mat)/rds)*med_rds
tpms <- log1p(mat)

categorical_var <- c("gender", "tobacco", "cc", "individual")
numerical_var <- c("age")
metadata <- seurat_obj@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor(metadata[,var])
}
metadata[,"cc"] <- stats::relevel(metadata[,"cc"], "0")

# create the SingleCellAssay (sca) object. See https://www.rdocumentation.org/packages/MAST/versions/0.931/topics/SingleCellAssay
sca <- MAST::FromMatrix(exprsArray = tpms,
                        cData = data.frame(wellKey = rownames(metadata),
                                           grp = metadata),
                        fData = data.frame(primerid = rownames(mat)))

# put ngeneson into the sca object
ngeneson <- apply(mat, 2, function(x) mean(x>0))
CD <- SummarizedExperiment::colData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
SummarizedExperiment::colData(sca) <- CD

set.seed(10)
mast_res <- MAST::zlm(formula = ~ grp.cc + (1 | grp.individual) + cngeneson + grp.age + grp.tobacco + grp.gender,
                      sca = sca,
                      method = "glmer",
                      ebayes = FALSE,
                      silent = T)
save(seurat_obj, sca, mast_res,
     date_of_run, session_info,
     file = "../../../out/simulation/simulation_1_mast.RData")

set.seed(10)
mast_lrt <- MAST::lrTest(mast_res, "grp.cc")
mast_pval_glmer <- apply(mast_lrt, 1, function(x){x[3,3]})
save(seurat_obj, sca, mast_res, mast_lrt, mast_pval_glmer,
     date_of_run, session_info,
     file =  "../../../out/simulation/simulation_1_mast.RData")

