rm(list=ls())
library(Seurat)
library(nebula)
library(foreach)
library(future)
library(rngtools)

load("../../../out/main/sns_layer23_processed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# https://github.com/lhe17/nebula
# see example usage in https://github.com/KellisLab/AD_regulome_analysis/blob/f641d471a79eacb2f379afd2fdb9c3b449f8867e/integration_linking_modules/deg.R

gene_vec <- sns[["RNA"]]@var.features
rm_idx <- grep("^MT", gene_vec)
if(length(rm_idx) > 0) gene_vec <- gene_vec[-rm_idx]

sns <- subset(sns, features = gene_vec)

neb_data <- nebula::scToNeb(obj = sns,
                            assay = "RNA",
                            id = "individual",
                            pred = c("region","sex","Seqbatch","age","diagnosis"),
                            offset = "nCount_RNA")
df <- model.matrix( ~region + sex + Seqbatch + age + diagnosis,
                    data = neb_data$pred)
nebula_res <- nebula::nebula(count = neb_data$count,
                             id = neb_data$id,
                             pred = df,
                             offset = neb_data$offset,
                             model = "NBGMM",
                             verbose = TRUE)

save(date_of_run, session_info, sns,
     nebula_res,
     time_start1, time_end1,
     file = "../../../out/main/sns_layer23_nebula.RData")





