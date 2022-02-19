rm(list=ls())
load("../../../../out/Writeup10/Writeup10_sns_invip_processed2.RData")
library(MAST)
library(ideas)

# see https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R
# following the analysis in https://github.com/himelmallick/BenchmarkSingleCell/blob/master/Library/run_MAST.R
# and https://www.bioconductor.org/packages/release/bioc/vignettes/MAST/inst/doc/MAITAnalysis.html

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat <- as.matrix(sns[["RNA"]]@counts[sns[["RNA"]]@var.features,])

cell_covariates <- data.frame(cell_id = rownames(sns@meta.data),
                              individual = sns@meta.data[,"individual"],
                              rd = covariates[,"Log_UMI"])
rownames(cell_covariates) <- colnames(mat)
indiv_vars <- c("diagnosis", "age", "sex", "Seqbatch", "RIN")
indiv_vec <- unique(sns@meta.data$individual)
indiv_covariates <- sapply(indiv_vec, function(indiv){
  idx <- which(sns@meta.data[,"individual"] == indiv)
  c(individual = indiv,
    diagnosis = sns@meta.data[idx[1],"diagnosis"],
    age = sns@meta.data[idx[1],"age"],
    Seqbatch = sns@meta.data[idx[1],"Seqbatch"],
    RIN = sns@meta.data[idx[1],"RNA.Integrity.Number"])
})
indiv_covariates <- as.data.frame(t(indiv_covariates))
indiv_covariates[,"individual"] <- as.factor(indiv_covariates[,"individual"])
indiv_covariates[,"diagnosis"] <- as.factor(indiv_covariates[,"diagnosis"])
indiv_covariates[,"Seqbatch"] <- as.factor(indiv_covariates[,"Seqbatch"])
indiv_covariates[,"age"] <- scale(as.numeric(indiv_covariates[,"age"]))
indiv_covariates[,"RIN"] <- scale(as.numeric(indiv_covariates[,"RIN"]))

# from https://github.com/Sun-lab/ideas_pipeline/blob/main/Autism/step1c_ideas.R
# https://github.com/Sun-lab/ideas_pipeline/blob/main/simulation/step2_evaluate_methods.R
dist1 <- ideas::ideas_dist(count_input = mat,
                           meta_cell = cell_covariates,
                           meta_ind = indiv_covariates,
                           var_per_cell = "rd",
                           var2test = "diagnosis",
                           var2test_type = "binary",
                           d_metric = "Was",
                           fit_method = "nb")

save(dist1, date_of_run, session_info,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_ideas.RData")

pval_res <- ideas::permanova(dist1,
                             meta_ind = indiv_covariates,
                             var2test = "diagnosis",
                             var2adjust = c("age", "sex", "Seqbatch", "RIN"),
                             var2test_type = "binary",
                             n_perm = 4999,
                             r.seed = 904)

save(dist1, pval_res, date_of_run, session_info,
     file = "../../../../out/Writeup10/Writeup10_sns_invip_ideas.RData")


