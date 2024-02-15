rm(list=ls())
library(Seurat)

global_suffix <- "~/kzlinlab/projects/eSVD2/out/main/"
file_names_list <- list(
  adams = list(
    prefix = "adams_",
    file_vec = c("T"),
    suffix = "_preprocessed.RData",
    object_name = "adams",
    cc_var = "Disease_Identity",
    individual_var = "Subject_Identity"
  ),
  habermann = list(
    prefix = "habermann_",
    file_vec = c("T"),
    suffix = "_preprocessed.RData",
    object_name = "habermann",
    cc_var = "Diagnosis",
    individual_var = "Sample_Name"
  ),
  regevEpi = list(
    prefix = "regevEpi_",
    file_vec = c("cyclingta-inflamed", "cyclingta-noninflamed",
                 "entprog-inflamed", "entprog-noninflamed",
                 "ta1-inflamed", "ta1-noninflamed",
                 "ta2-inflamed", "ta2-noninflamed"),
    suffix = "_esvd.RData",
    object_name = "regevEpi",
    cc_var = "Subject_Disease",
    individual_var = "Subject"
  ),
  sns = list(
    prefix = "sns_",
    file_vec = c("astpp", "endothelial", "insst", "invip", "layer23", "layer4",
                 "layer56cc", "layer56", "microglia", "oligo", "opc"),
    suffix = "_processed.RData",
    object_name = "sns",
    cc_var = "diagnosis",
    individual_var = "individual"
  )
)

for(lis in file_names_list){
  for(file in lis$file_vec){
    filename <- paste0(global_suffix, lis$prefix, file, lis$suffix)
    load(filename)

    if(lis$object_name == "adams"){
      seurat_obj <- adams
    } else if(lis$object_name == "habermann"){
      seurat_obj <- habermann
    } else if(lis$object_name == "regevEpi"){
      seurat_obj <- regevEpi
    } else if(lis$object_name == "sns"){
      seurat_obj <- sns
    } else {
      stop()
    }

    p <- length(Seurat::VariableFeatures(seurat_obj))
    n <- length(SeuratObject::Cells(seurat_obj))
    k <- length(unique(seurat_obj@meta.data[,lis$individual_var]))

    tab_mat <- table(seurat_obj@meta.data[,lis$individual_var],
                     seurat_obj@meta.data[,lis$cc_var])
    k0 <- length(which(tab_mat[,1] != 0))

    median_cell <- stats::median(table(seurat_obj@meta.data[,lis$individual_var]))

    mat <- SeuratObject::LayerData(seurat_obj,
                                   features = Seurat::VariableFeatures(seurat_obj),
                                   layer = "counts",
                                   assay = "RNA")
    sparsity <- length(mat@x)/prod(dim(mat))*100
    median_depth <- stats::median(seurat_obj$nCount_RNA)

    print(paste0(lis$prefix, file))
    print(paste0("p = ", p, ", n = ", n, ", k = ", k, ", k0 = ", k0))
    print(paste0("median_cell = ", median_cell, ", sparsity = ", round(sparsity,1),
                 ", median_depth = ", median_depth))
    print(head(tab_mat))
    print("======")
  }
}


