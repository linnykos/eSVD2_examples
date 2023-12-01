rm(list=ls())
library(Seurat)
library(nebula)
library(foreach)
library(future)
library(rngtools)

file_prefix <- "../../../out/main/sns_"
file_suffix <- "_processed.RData"
celltypes <- c("astpp", "endothelial", "insst", "invip", "layer4", "layer23",
               "layer56", "layer56cc", "microglia", "oligo", "opc")

# https://github.com/lhe17/nebula
# see example usage in https://github.com/KellisLab/AD_regulome_analysis/blob/f641d471a79eacb2f379afd2fdb9c3b449f8867e/integration_linking_modules/deg.R

for(kk in 1:length(celltypes)){

  celltype <- celltypes[kk]
  load(paste0(file_prefix, celltype, file_suffix))
  print(celltype)

  set.seed(10)
  date_of_run <- Sys.time()
  session_info <- devtools::session_info()

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
  save(nebula_res,
       date_of_run, session_info,
       file = paste0("../../../out/main/sns_", celltype, "_nebula.RData"))
}
