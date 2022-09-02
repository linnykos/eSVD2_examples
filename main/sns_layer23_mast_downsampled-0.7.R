rm(list=ls())
library(Seurat)
library(eSVD2)

load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

downsample_value <- 0.7
print(paste0("Working on downsample: ", downsample_value))
if("mat" %in% ls()) rm(list = "mat")

load(paste0("../../../out/main/sns_layer23_processed_downsampled-", downsample_value, ".RData"))
set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(mat)

rds <- colSums(mat)
med_rds <- median(rds)
mat <- t(t(mat)/rds)*med_rds
tpms <- log1p(mat)

categorical_var <- c("diagnosis", "individual", "region", "sex", "Seqbatch", "Capbatch")
numerical_var <- c("age", "percent.mt", "RNA.Integrity.Number", "post.mortem.hours")
metadata <- sns@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor(metadata[,var])
}
metadata[,"diagnosis"] <- stats::relevel(metadata[,"diagnosis"], "Control")

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
mast_res <- MAST::zlm(formula = ~ grp.diagnosis + (1 | grp.individual) + cngeneson + grp.region + grp.Seqbatch + grp.Capbatch + grp.sex + grp.age + grp.percent.mt + grp.RNA.Integrity.Number + grp.post.mortem.hours,
                      sca = sca,
                      method = "glmer",
                      ebayes = FALSE,
                      silent = T,
                      fitArgsC = list(control = lme4::glmerControl(calc.derivs = FALSE)),
                      fitArgsD = list(control = lme4::glmerControl(calc.derivs = FALSE), nAGQ = 0))
save(sns, sca, mast_res,
     date_of_run, session_info, downsample_value,
     file = paste0("../../../out/main/sns_layer23_mast_downsampled-", downsample_value, ".RData"))

set.seed(10)
mast_lrt <- MAST::lrTest(mast_res, "grp.diagnosis")
mast_pval_glmer <- apply(mast_lrt, 1, function(x){x[3,3]})
save(sns, sca, mast_res, mast_lrt, mast_pval_glmer,
     date_of_run, session_info, downsample_value,
     file = paste0("../../../out/main/sns_layer23_mast_downsampled-", downsample_value, ".RData"))

