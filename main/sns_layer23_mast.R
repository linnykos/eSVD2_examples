rm(list=ls())
library(SummarizedExperiment)
library(MAST)
library(lme4)

load("../../../out/main/sns_layer23_mast-prepared.RData") # loads mat, metadata, categorical_var, numerical_var,

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# see https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R
# following the analysis in https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/run_MAST.R
# and https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html

mat <- as.matrix(mat)
rds <- colSums(mat)
med_rds <- median(rds)
mat <- t(t(mat)/rds)*med_rds
tpms <- log1p(mat)

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
                      silent = T)
save(mast_res,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_mast.RData")

set.seed(10)
mast_lrt <- MAST::lrTest(mast_res, "grp.diagnosis")
mast_pval_glmer <- apply(mast_lrt, 1, function(x){x[3,3]})
save(mast_res, mast_lrt, mast_pval_glmer,
     date_of_run, session_info,
     file = "../../../out/main/sns_layer23_mast.RData")
