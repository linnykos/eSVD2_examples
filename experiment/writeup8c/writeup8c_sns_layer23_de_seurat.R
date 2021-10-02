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

## according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7678724/bin/NIHMS1053005-supplement-supplement.pdf,
# we want to regress out the following:
# [[age,  sex,  cortical  region, RIN  and  post-mortem  interval,
# as  well  as  10X  capture  and  sequencing batch and per-cell ribosomal RNA fraction]]
# a lot of these are categorical, so let's make indicators.
# We'll do a full set of indicators

categorical_var <- c("region", "sex", "Capbatch", "Seqbatch")
numerical_var <- c("age", "RNA.Integrity.Number", "post.mortem.hours", "percent.mt", "nFeature_RNA")
new_indicator_var <- c()
n <- ncol(sns)

for(variable in categorical_var){
  covariate <- sns@meta.data[,variable]
  uniq_level <- unique(covariate)
  for(i in uniq_level[-1]){
    tmp <- rep(0, n)
    tmp[which(covariate == i)] <- 1

    var_name <- paste0(variable, "_", i)
    sns[[var_name]] <- tmp
    new_indicator_var <- c(new_indicator_var, var_name)
  }
}

vars_to_regress <- c(numerical_var, new_indicator_var)
Seurat::DefaultAssay(sns) <- "RNA"
gene_names <- rownames(sns)
gene_names <- gene_names[-grep("^MT-", gene_names)]
set.seed(10)
sns <- Seurat::SCTransform(sns,
                           residual.features = gene_names,
                           vars.to.regress = vars_to_regress,
                           verbose = T)

set.seed(10)
sns_de <- Seurat::FindMarkers(sns,
                              assay = "SCT",
                              slot = "scale.data",
                              ident.1 = "Control",
                              ident.2 = "ASD",
                              group.by = "diagnosis",
                              logfc.threshold = 0,
                              min.pct = 0,
                              min.cells.feature = 0,
                              verbose = T)

save(sns, sns_de,
     file = "../../../../out/writeup8c/writeup8c_sns_layer23_de_seurat.RData")
