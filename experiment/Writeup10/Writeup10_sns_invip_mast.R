rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
library(MAST)
library(lme4)
library(Seurat)

# see https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R
# following the analysis in https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/run_MAST.R
# and https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])
rds <- colSums(mat)
med_rds <- median(rds)
mat <- t(t(mat)/rds)*med_rds
tpms <- log1p(mat)

categorical_var <- c("diagnosis", "individual", "region", "sex", "Capbatch", "Seqbatch")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
metadata <- sns@meta.data[,c(categorical_var, numerical_var)]
for(var in categorical_var){
  metadata[,var] <- as.factor( metadata[,var])
}
metadata[,"diagnosis"] <- stats::relevel(metadata[,"diagnosis"], "Control")

# create the SingleCellAssay (sca) object. See https://www.rdocumentation.org/packages/MAST/versions/0.931/topics/SingleCellAssay
sca <- MAST::FromMatrix(exprsArray = tpms,
                        cData = data.frame(wellKey = rownames(metadata),
                                           grp = metadata),
                        fData = data.frame(primerid = rownames(mat)))

# put ngeneson into the sca object
ngeneson <- apply(mat, 2, function(x) mean(x>0))
CD <- colData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
colData(sca) <- CD

set.seed(10)
mast_res <- MAST::zlm(formula = ~ grp.diagnosis + (1 | grp.individual) + cngeneson + grp.RNA.Integrity.Number + grp.post.mortem.hours + grp.region + grp.Capbatch + grp.Seqbatch,
                      sca = sca,
                      method = "glmer",
                      ebayes = FALSE,
                      silent = F)
save(sns, sca, mast_res,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_mast.RData")

set.seed(10)
mast_lrt <- MAST::lrTest(mast_res, "grp.diagnosis")
save(sns, sca, mast_res, mast_lrt,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_mast.RData")

set.seed(10)
fit <- MAST::lrTest(zlmdata, CoefficientHypothesis("grp.diagnosisASD"))
save(sns, sca, mast_res, mast_lrt, fit,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_mast.RData")
